"""Conditional GR-form quadrupole equation on checked full-law backgrounds.
Piecewise log-pressure integration resolves every EOS interpolation knot.
Schwarzschild extraction and cutoff matching are explicit prescriptions.
"""
from pathlib import Path
import sys
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from decimal import Decimal,localcontext
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'paper3_NeutronStars/ns_step14'))
from model import EOS,G,C,MSUN,L,delta,solve as background
from resource_control import require_controls
KG=G/C**4

def dedp(eos,P,idx=None):
 eps=eos.eps_of_P(P)
 if getattr(eos,'analytic',False):return eos.deps_dP(P)
 if eos.interp=='pchip':s=eos.pe.derivative()(np.log(P))
 else:
  i=int(np.clip(np.searchsorted(eos.lp,np.log(P),side='right')-1,0,len(eos.lp)-2)) if idx is None else idx
  s=(eos.le[i+1]-eos.le[i])/(eos.lp[i+1]-eos.lp[i])
 val=eps*s/P
 if not np.isfinite(val) or val<=0:raise RuntimeError('Invalid EOS energy derivative')
 return val

def love_k2(c,y):
 if not 0<c<.5 or not np.isfinite(y):raise ValueError('Invalid compactness/y')
 # Decimal prevents cancellation in the Newtonian benchmark and extraction.
 with localcontext() as ctx:
  ctx.prec=60;c=Decimal(str(c));y=Decimal(str(y));f=1-2*c;b=2+2*c*(y-1)-y
  num=Decimal(8)/5*c**5*f**2*b
  den=2*c*(6-3*y+3*c*(5*y-8))+4*c**3*(13-11*y+c*(3*y-2)+2*c*c*(1+y))+3*f*f*b*f.ln()
  k=float(num/den)
 if not np.isfinite(k) or k<=0:raise RuntimeError('Invalid Love number')
 return k

def solve(eos,Pc,A=0.,alpha=.12,embedding='emb1',rstar=1e6,r0=100.,rtol=2e-9,atol=1e-17,max_step=.08,surface_factor=1.001,form='y',metric='hydrostatic',segment=True):
 require_controls()
 if embedding not in ('emb1','emb2') or form not in ('y','H') or metric not in ('hydrostatic','einstein'):raise ValueError('Unknown prescription')
 if not all(np.isfinite([Pc,A,alpha,rstar,r0,rtol,atol,max_step,surface_factor])) or not (0<=A<.1 and alpha>0 and rstar>0 and r0>0 and rtol>0 and atol>0 and max_step>0 and surface_factor>=1):raise ValueError('Invalid inputs')
 Ps=eos.P_min*surface_factor
 if not Ps<Pc<=eos.P_max:raise ValueError('EOS domain')
 ec,rc=eos.values(Pc);dc=delta(rc,A,alpha)
 B=2*np.pi*G/C**4*(ec+Pc)*(ec/3+Pc)*3/(3+dc) if embedding=='emb1' else 2*np.pi*G/C**4*(ec+Pc)*(ec/(3+dc)+Pc)
 drop=B*r0*r0/Pc
 if not 0<drop<.01:raise ValueError('Central expansion')
 t0=-np.log1p(-drop);t1=np.log(Pc/Ps)
 if t0>=t1:raise ValueError('Start beyond surface')
 m=4*np.pi/3*r0**3*ec/C**2/MSUN;b0=4*np.pi/3*r0**3*rc/MSUN;h=(r0/rstar)**dc
 y2=-4*np.pi*KG/7*(ec/3+11*Pc+(ec+Pc)*dedp(eos,Pc))
 state=[r0/L,m,b0,b0*h*3/(3+dc),m*h*3/(3+dc)]
 state += [2+y2*r0*r0] if form=='y' else [1+y2*r0*r0/2,y2*r0*r0]
 def terms(t,z):
  r=z[0]*L;m=z[1];P=Pc*np.exp(-t);eps,rho=eos.values(P);mg=G*MSUN*m/C**2;f=1-2*mg/r
  if not np.all(np.isfinite(z)) or min(z[:5])<=0 or f<=0:raise RuntimeError('Invalid background state')
  w=np.exp(delta(rho,A,alpha)*np.log(r/rstar));ratio=z[3]/(w*z[2]);press=4*np.pi*r**3*P/(C*C*MSUN)
  active=(m+press)*ratio if embedding=='emb1' else z[4]/w+press
  g=G*MSUN*active/(C*C*r*r*f);dp=-(eps+P)*g;dr=-P/dp
  if not np.isfinite(dr) or dr<=0:raise RuntimeError('Invalid pressure gradient')
  return r,P,eps,rho,mg,f,w,g,dr,ratio,active
 nfev=0;steps=0;ends=[];z=np.array(state)
 knots=[] if not segment or getattr(eos,'analytic',False) else [np.log(Pc/p) for p in eos.P if Ps<p<Pc*np.exp(-t0)]
 edges=np.unique([t0,*knots,t1]);idx=None
 for left,right in zip(edges[:-1],edges[1:]):
  if not getattr(eos,'analytic',False) and eos.interp=='linear' and segment:
   idx=int(np.clip(np.searchsorted(eos.lp,np.log(Pc)-(left+right)/2)-1,0,len(eos.lp)-2))
  def rhs(t,z):
   r,P,eps,rho,mg,f,w,g,dr,ratio,active=terms(t,z);vol=4*np.pi*r*r*dr/MSUN
   out=[dr/L,vol*eps/C**2,vol*rho,vol*w*rho,vol*w*eps/C**2]
   pg=KG*P;eg=KG*eps
   nu=2*g if metric=='hydrostatic' else 2*(mg+4*np.pi*r**3*pg)/(r*r*f)
   F=(1+4*np.pi*r*r*(pg-eg))/f
   q=(4*np.pi*r*r*(5*eg+9*pg+(eg+pg)*dedp(eos,P,idx))-6)/f-(r*nu)**2
   if form=='y':out += [-(z[5]**2+F*z[5]+q)*dr/r]
   else:out += [z[6]*dr/r,(-(4+F)*z[6]-(4+2*F+q)*z[5])*dr/r]
   return out
  sol=solve_ivp(rhs,(left,right),z,method='DOP853',rtol=rtol,atol=[atol]*5+[max(atol,rtol*.01)]*(len(z)-5),max_step=max_step)
  if not sol.success or sol.status!=0 or sol.t[-1]!=right or not np.all(np.isfinite(sol.y)):raise RuntimeError('Incomplete tidal segment: '+sol.message)
  z=sol.y[:,-1];nfev+=sol.nfev;steps+=len(sol.t);ends.append({'t':float(right),'state':z.tolist()})
 r,P,eps,rho,mg,f,w,g,dr,ratio,active=terms(t1,z)
 y=float(z[5] if form=='y' else 2+z[6]/z[5]);compact=mg/r
 # Distinct cutoff extraction conventions; none asserts a physical shell solution.
 jumps={'GR_zero_pressure':4*np.pi*r**3*KG*eps/mg,
        'GR_finite_pressure':4*np.pi*r**3*KG*eps/(mg+4*np.pi*r**3*KG*P),
        'frozen_hydrostatic':4*np.pi*r*KG*eps/(f*g),'none':0.}
 extraction={}
 for name,jump in jumps.items():
  yy=y-jump;k=love_k2(compact,yy);extraction[name]={'jump':jump,'y_out':yy,'k2':k,'Lambda':2/3*k/compact**5}
 primary=extraction['GR_finite_pressure']
 return {'M':float(z[1]),'R':float(r/1e5),'Pc':float(Pc),'P_surface':float(Ps),'C':float(compact),'y_in':y,**primary,'extractions':extraction,'nu_prime_in':2*g,'nu_prime_Schwarzschild':2*mg/(r*r*f),'metric_derivative_relative_mismatch':g/(mg/(r*r*f))-1,'surface_operator_ratio':float(ratio),'surface_pressure_to_energy':float(P/eps),'nfev':nfev,'steps':steps,'segments':ends}

def fixed(eos,target,A,alpha,embedding='emb1',**settings):
 require_controls();evs=[];cache={}
 # Root the checked five-state background, then verify augmented mass/radius.
 bg={k:v for k,v in settings.items() if k in ('rstar','r0','rtol','atol','max_step','surface_factor')}
 bg.update(rtol=min(bg.get('rtol',2e-9),2e-11),atol=min(bg.get('atol',1e-17),1e-19),max_step=min(bg.get('max_step',.08),.04))
 def ev(lp):
  if lp not in cache:
   cache[lp]=background(eos,float(np.exp(lp)),A,alpha,embedding,**bg);evs.append(cache[lp])
  return cache[lp]
 lo,hi=np.log(2e34),np.log(min(7e35,eos.P_max))
 if not ev(lo)['M']<target<ev(hi)['M']:raise RuntimeError('Target mass not bracketed')
 root=brentq(lambda lp:ev(lp)['M']-target,lo,hi,xtol=2e-10)
 bgstar=ev(root);star=solve(eos,float(np.exp(root)),A,alpha,embedding,**settings)
 if abs(star['M']-target)>5e-7 or abs(star['M']/bgstar['M']-1)>5e-7 or abs(star['R']/bgstar['R']-1)>5e-7:raise RuntimeError('Augmented background mismatch')
 return {'target_mass':target,'A':A,'alpha':alpha,'embedding':embedding,'settings':settings,'star':star,'background':bgstar,'root_evaluations':evs}

def pair(gr,modified):
 a=gr['star'];b=modified['star'];kr=b['k2']/a['k2'];rr=b['R']/a['R'];mr=a['M']/b['M'];lr=b['Lambda']/a['Lambda']
 exact=kr*rr**5*mr**5
 if abs(exact-lr)>1e-12:raise RuntimeError('Decomposition mismatch')
 return {'delta_Lambda':lr-1,'delta_k2':kr-1,'delta_R':rr-1,'radius_only_exact':rr**5-1,'radius_only_linear':5*(rr-1),'multiplicative_cross_term':(kr-1)*(rr**5-1),'mass_residual_factor':mr**5,'identity_residual':exact-lr}
