"""Independent radial-field cross-check of the pressure-coordinate flux solver.
Analytic derivative of the selected EOS interpolant or the legacy .001P stencil.
Stage evaluations below the surface use constant endpoint EOS solely to locate
an event; no accepted output below the declared EOS floor is returned.
"""
import numpy as np
from scipy.integrate import solve_ivp
from model import G,C,MSUN,L,delta
S=G*MSUN/L**2

def solve_direct(eos,Pc,A=0.,alpha=.12,embedding='emb1',r0=100.,rstar=1e6,
                 rtol=2e-10,atol=1e-17,max_step=2000.,surface_factor=1.001,derivative='exact',rmax=3e6):
    if embedding not in ('emb1','emb2') or derivative not in ('exact','fd'):raise ValueError('mode')
    if not all(np.isfinite([Pc,A,alpha,r0,rstar,rtol,atol,max_step,surface_factor,rmax])) or not (0<=A<.1 and alpha>0 and r0>0 and rstar>0 and rtol>0 and atol>0 and max_step>0 and surface_factor>=1 and rmax>r0):raise ValueError('invalid inputs')
    Ps=eos.P_min*surface_factor
    if not Ps<Pc<=eos.P_max:raise ValueError('EOS range')
    ec,rc=eos.values(Pc);dc=delta(rc,A,alpha)
    B=2*np.pi*G/C**4*(ec+Pc)*(ec/3+Pc)*3/(3+dc) if embedding=='emb1' else 2*np.pi*G/C**4*(ec+Pc)*(ec/(3+dc)+Pc)
    if B*r0*r0/Pc>=.01:raise ValueError('central series')
    m0=4*np.pi/3*r0**3*ec/C**2/MSUN;h0=(r0/rstar)**dc
    y0=[np.log(Pc-B*r0*r0),m0]
    y0+= [4*np.pi*G*rc*r0/(3+dc)/S,4*np.pi*G*rc*r0/3/S] if embedding=='emb1' else [m0*h0*3/(3+dc)]
    def rhs(u,y):
        r=u*L;P=np.exp(y[0]);eps,rho=eos.values(max(P,Ps));m=y[1]*MSUN
        if not np.all(np.isfinite(y)) or min(y[1:])<=0:raise RuntimeError('invalid state')
        f=1-2*G*m/(r*C*C)
        if f<=0:raise RuntimeError('trapped surface')
        d=delta(rho,A,alpha);h=(r/rstar)**d
        if embedding=='emb1':
            R=y[2]/y[3];dp=-G*(eps+P)*(m+4*np.pi*r**3*P/C**2)/(C*C*r*r*f)*R
        else:dp=-G*(eps+P)*(MSUN*y[2]/h+4*np.pi*r**3*P/C**2)/(C*C*r*r*f)
        out=[L*dp/P,L*4*np.pi*r*r*eps/C**2/MSUN]
        if embedding=='emb1':
            if P<Ps:drho=0.
            elif derivative=='fd':
                hi=min(P*1.001,eos.P_max);lo=max(P*.999,eos.P_min)
                drho=(eos.rho_of_P(hi)-eos.rho_of_P(lo))/(hi-lo)
            else:drho=eos.drho_dP(P)
            x=(rho/150)**alpha;dd=-3*A*alpha*x/(rho*(1+x)**2)*drho*dp
            src=4*np.pi*G*rho/S
            out += [L*(src-((2+d)/r+dd*np.log(r/rstar))*y[2]),L*(src-2*y[3]/r)]
        else:out += [out[1]*h]
        return out
    def surf(u,y):return y[0]-np.log(Ps)
    surf.terminal=True;surf.direction=-1
    sol=solve_ivp(rhs,(r0/L,rmax/L),y0,method='DOP853',rtol=rtol,atol=atol,max_step=max_step/L,events=surf)
    if not sol.success or sol.status!=1 or len(sol.t_events[0])!=1 or not np.all(np.isfinite(sol.y)):raise RuntimeError('no successful unique surface: '+sol.message)
    end=sol.y_events[0][0];r=sol.t_events[0][0]*L
    if abs(end[0]-np.log(Ps))>1e-10 or min(end[1:])<=0 or 2*G*MSUN*end[1]/(r*C*C)>=1:raise RuntimeError('bad surface')
    return {'M':float(end[1]),'R':float(r/1e5),'Pc':float(Pc),'P_surface':float(Ps),'steps':len(sol.t),'nfev':sol.nfev}
