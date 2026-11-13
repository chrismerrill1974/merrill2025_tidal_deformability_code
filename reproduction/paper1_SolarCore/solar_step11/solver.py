"""Step 11 solar-a response: ODE quadratures and dense profile evaluation.
Same expanded operator and explicit surface regularization as solarcore_v3.
"""
from dataclasses import dataclass, asdict
from pathlib import Path
import sys
import math
import numpy as np
from scipy.integrate import solve_ivp
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from numerical_guardrails import positive

BETA=np.array([4.,2.,10.,20.])
@dataclass(frozen=True)
class Settings:
    rtol: float=2e-12
    atol: float=2e-14
    max_step: float=.04
    xi0: float=1e-6
    floor: float=1e-8
    calibration_tol: float=1e-12
    max_iter: int=60
    reference: str='a'
    si: bool=False
    profile_n: int=801


def integrate(A,alpha,rho,rho0,L,s):
    positive(rho=rho,rho0=rho0,alpha=alpha,rtol=s.rtol,atol=s.atol,max_step=s.max_step,xi0=s.xi0,floor=s.floor)
    if not 0<=A<1 or not math.isfinite(L):raise ValueError('Invalid model')
    def delta(th):
        x=(rho*max(th,1e-30)**3/rho0)**alpha
        return -3*A*x/(1+x), -9*A*alpha*x/(max(th,1e-30)*(1+x)**2)
    def rhs(x,y):
        th=max(y[0],0.);d,dd=delta(th)
        coeff=(2+d)/x+(dd*y[1]*(math.log(x)-L) if th>s.floor else 0.)
        return [y[1],-th**3-coeff*y[1],x*x*th**3,*[x*x*th**g for g in (10,8,16,26)]]
    d,_=delta(1.);x=s.xi0
    y=[1-x*x/(2*(3+d)),-x/(3+d),*[x**3/3]*5]
    def surface(x,y):return y[0]
    surface.terminal=True;surface.direction=-1
    sol=solve_ivp(rhs,(x,20),y,method='DOP853',rtol=s.rtol,atol=s.atol,max_step=s.max_step,events=surface,dense_output=True)
    if not sol.success:raise RuntimeError('ODE solver failed: '+sol.message)
    if sol.status!=1 or len(sol.t_events[0])!=1 or not np.isfinite(sol.y).all():raise RuntimeError('Invalid surface')
    if abs(sol.y[0,-1])>1e-8 or sol.y[1,-1]>=0:raise RuntimeError('Invalid endpoint')
    return sol


def calibrate(A,alpha,s=Settings()):
    if s.reference not in ('a','Rsun'):raise ValueError('Unknown reference')
    M,R,rho0=(1.98847e30,6.957e8,150e3) if s.si else (1.98847e33,6.957e10,150.)
    rho=80*rho0/150;L=math.log(6.897) if s.reference=='Rsun' else 0.
    for iteration in range(s.max_iter):
        sol=integrate(A,alpha,rho,rho0,L,s);xi=sol.t[-1];a=R/xi
        I=sol.y[2,-1];positive(I=I)
        new=M/(4*np.pi*a**3*I);ln=math.log(xi) if s.reference=='Rsun' else 0.
        residual=max(abs(new/rho-1),abs(ln-L))
        if residual<s.calibration_tol:
            return dict(sol=sol,xi=xi,a=a,rho=rho,T=a*a*rho,kernels=sol.y[3:,-1],iterations=iteration+1,residual=residual)
        rho,L=new,ln
    raise RuntimeError(f'Calibration failed A={A}, alpha={alpha}, residual={residual}')


def response(A,alpha,s=Settings(),std=None):
    if std is None:std=calibrate(0,.1,s)
    mod=std if A==0 else calibrate(A,alpha,s)
    ratios=(mod['rho']/std['rho'])**2*(mod['a']/std['a'])**3*(mod['T']/std['T'])**BETA*mod['kernels']/std['kernels']
    x=np.linspace(.01,.98,s.profile_n)
    t=mod['sol'].sol(x*mod['xi'])[0];t0=std['sol'].sol(x*std['xi'])[0]
    if np.any(t<=0) or np.any(t0<=0):raise RuntimeError('Incomplete sound profile')
    h=float(np.max(np.abs(np.sqrt(mod['T']/std['T']*t/t0)-1)))
    row=dict(A=float(A),alpha=float(alpha),xi1=mod['xi'],rho_c=mod['rho'],Tc_ratio=mod['T']/std['T'],h=h,iterations=mod['iterations'],calibration_residual=mod['residual'])
    row.update(dict(zip(('R_pp','R_pep','R_be7','R_b8'),map(float,ratios))))
    if not all(math.isfinite(v) for v in row.values()):raise RuntimeError('Nonfinite response')
    return row
