"""Step 14 independent pressure-coordinate flux solver. CGS inputs.
Ordinary mass, baryon integral, weighted baryon and weighted energy integrals
are separate states; no numerical rho derivative enters the primary solution.
No EOS extrapolation, A=0 shortcut, or fallback on failed integrations.
"""
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import PchipInterpolator
G=6.67430e-8; C=2.99792458e10; MSUN=1.98847e33; L=1e6
ROOT=Path(__file__).resolve().parents[2]
class EOS:
    def __init__(self,kind='compose',interp='linear'):
        if interp not in ('linear','pchip'):raise ValueError(interp)
        self.kind,self.interp=kind,interp
        if kind=='compose':
            path=ROOT/'paper3_NeutronStars/ns_v7'
            words=(path/'eos.nb').read_text().split();lo,hi=map(int,words[:2]);nb=np.array(words[2:],float)
            if len(nb)!=hi-lo+1:raise ValueError('density grid length')
            with (path/'eos.thermo').open() as f:
                mn=float(f.readline().split()[0]);rows=np.loadtxt(f)
            ids=rows[:,1].astype(int)
            if len(set(ids))!=len(ids) or set(ids)!=set(range(lo,hi+1)):raise ValueError('EOS index mapping')
            n=nb[ids-lo];self.P=rows[:,3]*n*1.602176634e33
            self.eps=n*mn*(1+rows[:,9])*1.602176634e33
            self.rho=n*1e39*1.66053906660e-24
        elif kind in ('dh','dh_bug'):
            a=np.loadtxt(ROOT/'paper3_NeutronStars/037/eos_sly4_dh2001_v2.dat')
            self.rho=a[:,0];self.P=a[:,1];self.eps=a[:,0]*C*C if kind=='dh' else a[:,2]
        else:raise ValueError(kind)
        order=np.argsort(self.P);self.P,self.eps,self.rho=[v[order] for v in (self.P,self.eps,self.rho)]
        if not all(np.all(np.isfinite(v)) and np.all(v>0) and np.all(np.diff(v)>0) for v in (self.P,self.eps,self.rho)):raise ValueError('nonmonotone EOS')
        self.lp=np.log(self.P);self.le=np.log(self.eps);self.lr=np.log(self.rho)
        self.pe=PchipInterpolator(self.lp,self.le,extrapolate=False);self.pr=PchipInterpolator(self.lp,self.lr,extrapolate=False)
        self.P_min=float(self.P[0]);self.P_max=float(self.P[-1])
    def values(self,P):
        if not np.isfinite(P) or P<self.P_min*(1-1e-13) or P>self.P_max*(1+1e-13):raise ValueError('EOS out of range')
        p=np.clip(np.log(P),self.lp[0],self.lp[-1])
        if self.interp=='linear':return np.exp(np.interp(p,self.lp,self.le)),np.exp(np.interp(p,self.lp,self.lr))
        return np.exp(self.pe(p)),np.exp(self.pr(p))
    def eps_of_P(self,P):return self.values(P)[0]
    def rho_of_P(self,P):return self.values(P)[1]
    def drho_dP(self,P):
        _,rho=self.values(P);p=np.log(P)
        if self.interp=='pchip':slope=self.pr.derivative()(p)
        else:
            i=np.clip(np.searchsorted(self.lp,p,side='right')-1,0,len(self.lp)-2)
            slope=(self.lr[i+1]-self.lr[i])/(self.lp[i+1]-self.lp[i])
        return rho*slope/P

def delta(rho,A,alpha):
    z=np.exp(alpha*np.log(rho/150.));return -3*A*z/(1+z)

def solve(eos,Pc,A=0.,alpha=.12,embedding='emb1',r0=100.,rstar=1e6,
          rtol=2e-9,atol=1e-17,max_step=.1,surface_factor=1.001,method='DOP853',legacy_seed=False):
    vals=[Pc,A,alpha,r0,rstar,rtol,atol,max_step,surface_factor]
    if not np.all(np.isfinite(vals)) or not (0<=A<.1 and alpha>0 and r0>0 and rstar>0 and rtol>0 and atol>0 and max_step>0 and surface_factor>=1):raise ValueError('invalid solver inputs')
    if embedding not in ('emb1','emb2'):raise ValueError('unknown embedding')
    Ps=eos.P_min*surface_factor
    if not Ps<Pc<=eos.P_max:raise ValueError('central/surface pressure outside EOS')
    ec,rc=eos.values(Pc);dc=delta(rc,A,alpha)
    Bgr=2*np.pi*G/C**4*(ec+Pc)*(ec/3+Pc)
    B=Bgr*3/(3+dc) if embedding=='emb1' else 2*np.pi*G/C**4*(ec+Pc)*(ec/(3+dc)+Pc)
    if legacy_seed:B=Bgr
    drop=B*r0*r0/Pc
    if not 0<drop<.01:raise ValueError('central expansion out of range')
    t0=-np.log1p(-drop);t1=np.log(Pc/Ps)
    if t0>=t1:raise ValueError('start outside surface')
    m=4*np.pi/3*r0**3*ec/C**2/MSUN;b0=4*np.pi/3*r0**3*rc/MSUN
    h=(r0/rstar)**dc;b=b0*h*3/(3+dc);e=m*h*3/(3+dc)
    if legacy_seed and embedding=='emb2':e=m*h
    y0=[r0/L,m,b0,b,e]
    def rhs(t,y):
        u,m,b0,b,e=y;r=u*L
        if not np.all(np.isfinite(y)) or min(y)<=0:raise RuntimeError('nonpositive/nonfinite state')
        f=1-2*G*MSUN*m/(r*C*C)
        if f<=0:raise RuntimeError('trapped surface')
        P=Pc*np.exp(-t);eps,rho=eos.values(P);h=np.exp(delta(rho,A,alpha)*np.log(r/rstar))
        active=m*(b/(h*b0)) if embedding=='emb1' else e/h
        press=4*np.pi*r**3*P/(C*C*MSUN)
        # In embedding 1 the pressure term is multiplied by the same ratio.
        if embedding=='emb1':press*=b/(h*b0)
        grad=-G*(eps+P)*MSUN*(active+press)/(C*C*r*r*f)
        dr=-P/grad
        if not np.isfinite(dr) or dr<=0:raise RuntimeError('invalid gradient')
        vol=4*np.pi*r*r*dr/MSUN
        return [dr/L,vol*eps/C**2,vol*rho,vol*h*rho,vol*h*eps/C**2]
    sol=solve_ivp(rhs,(t0,t1),y0,method=method,rtol=rtol,atol=atol,max_step=max_step)
    if not sol.success or sol.status!=0 or sol.t[-1]!=t1 or not np.all(np.isfinite(sol.y)):raise RuntimeError('integration did not reach pressure surface: '+sol.message)
    u,m,b0,b,e=sol.y[:,-1];r=u*L
    if min(sol.y.ravel())<=0 or 2*G*MSUN*m/(r*C*C)>=1:raise RuntimeError('invalid completed star')
    return {'M':float(m),'R':float(r/1e5),'Pc':float(Pc),'P_surface':float(Ps),'steps':len(sol.t),'nfev':sol.nfev}
