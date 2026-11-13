"""Audit existing NS RHSs without rewriting production or running an EOS survey.
Remove only the A==0 dispatch in an in-memory AST; capture actual nested RHS/IC.
Retain source hashes and measured discrepancies, including expected failures.
"""
from pathlib import Path
import ast, hashlib, importlib.util, json, platform, sys
import numpy as np
import scipy
from scipy.integrate import quad
ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'2026_june_14_make_it_consistent/paper3_keystone/tov_operator.py'
spec=importlib.util.spec_from_file_location('step13_operator',SOURCE)
op=importlib.util.module_from_spec(spec);spec.loader.exec_module(op)
tree=ast.parse(SOURCE.read_text());fn=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='solve_star_op')
dispatch=fn.body[1]
assert isinstance(dispatch,ast.If) and ast.unparse(dispatch.test)=='A == 0.0'
fn.body.pop(1)
ns=op.__dict__.copy();exec(compile(ast.fix_missing_locations(ast.Module(body=[fn],type_ignores=[])),str(SOURCE),'exec'),ns)
solve=ns['solve_star_op'];G,C=op.G,op.C
class EOS:
    P_min=1e18
    def __init__(self,c=C,K=1e5,energy_ratio=1.):self.c,self.K,self.energy_ratio=c,K,energy_ratio
    def rho_of_P(self,P):return np.sqrt(P/self.K)
    def eps_of_P(self,P):return self.energy_ratio*self.rho_of_P(P)*self.c**2
class Captured(Exception):pass
def capture(eos,Pc,A,alpha,emb,rstar=1e6):
    data={}
    def fake(rhs,span,y0,**kw):
        data.update(rhs=rhs,span=span,y0=np.array(y0));raise Captured
    original=ns['solve_ivp'];ns['solve_ivp']=fake
    try:
        try:solve(eos,Pc,A,alpha,embedding=emb,rstar_cm=rstar)
        except Captured:pass
    finally:ns['solve_ivp']=original
    assert data
    return data

def run():
    out={'scope':'RHS/central-series audit and synthetic EOS limit checks; no SLy inference', 'runtime':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__},'checks':{}}
    eos=EOS();Pc=1e33;r=3e5;rho=eos.rho_of_P(Pc);eps=eos.eps_of_P(Pc);m=4*np.pi*r**3*rho/3
    aN=G*m/r**2
    expected=-G*(eps+Pc)*(m+4*np.pi*r**3*Pc/C**2)/(C**2*r**2*(1-2*G*m/(r*C**2)))
    gr={}
    for emb,state in [('emb1',[Pc,m,aN,aN]),('emb2',[Pc,m,m])]:
        d=capture(eos,Pc,0.,.12,emb);rhs=d['rhs'](r,state)
        gr[emb]={'rhs_relative_error':abs(rhs[0]/expected-1)}
        baseline=op.solve_star_v7(eos,Pc,rtol=1e-10,max_step=2000.)
        actual=solve(eos,Pc,0.,.12,embedding=emb,rtol=1e-10,max_step=2000.)
        gr[emb]['synthetic_M_R_relative_errors']=list(np.abs(np.array(actual)/baseline-1))
        assert gr[emb]['rhs_relative_error']<1e-14
        assert max(gr[emb]['synthetic_M_R_relative_errors'])<1e-6
    out['checks']['zero_A_modified_paths']=gr
    central={};A=.0087;alpha=.122;dc=op.delta(rho,A,alpha)
    for emb in ['emb1','emb2']:
        d=capture(eos,Pc,A,alpha,emb);r0=d['span'][0]
        B=2*np.pi*G/C**4*(eps+Pc)*(eps/3+Pc)*3/(3+dc) if emb=='emb1' else 2*np.pi*G/C**4*(eps+Pc)*(eps/(3+dc)+Pc)
        coded=(Pc-d['y0'][0])/r0**2
        row={'delta_c':dc,'required_pressure_quadratic_coefficient':B,'coded_pressure_quadratic_coefficient':coded,'relative_coefficient_error':coded/B-1}
        assert abs(row['relative_coefficient_error'])>1e-3
        if emb=='emb2':
            exact=4*np.pi*eps/C**2*r0**3*(r0/1e6)**dc/(3+dc)
            row['weighted_mass_seed_relative_error']=d['y0'][2]/exact-1
            assert abs(row['weighted_mass_seed_relative_error']-dc/3)<1e-14
        central[emb]=row
    out['checks']['known_central_series_failures']=central
    # Smooth prescribed density profile; auxiliary fields from independent quadrature.
    L=2e5;r=3e5;rhoc=1e10;A=.0087;alpha=.02;rstar=1e6
    rh=lambda s:rhoc/(1+(s/L)**2)
    w=lambda s:(s/rstar)**op.delta(rh(s),A,alpha)
    mass=4*np.pi*quad(lambda s:s*s*rh(s),0,r,epsabs=0,epsrel=1e-11)[0]
    weighted=4*np.pi*quad(lambda s:s*s*w(s)*rh(s),0,r,epsabs=0,epsrel=1e-11)[0]
    g=G*weighted/(r*r*w(r));gN=G*mass/r**2
    weak={}
    for emb in ['emb1','emb2']:
        errors=[]
        for multiplier in [1.,10.,100.]:
            c=C*multiplier;ns['C']=c;e=EOS(c=c);P=e.K*rh(r)**2
            d=capture(e,P,A,alpha,emb,rstar)
            state=[P,mass,g,gN] if emb=='emb1' else [P,mass,weighted]
            grad=d['rhs'](r,state)[0];errors.append(abs(grad/(-rh(r)*g)-1))
        assert errors[2]<errors[1]<errors[0] and errors[-1]<1e-9
        weak[emb]={'relative_errors_for_c_multipliers_1_10_100':errors}
    ns['C']=C
    out['checks']['cold_Newtonian_limit']=weak
    # Different physical reference lengths versus a mere unit conversion.
    refs={}
    for constant in [False,True]:
        vals=[]
        for ref in [1e5,1e6,1e7]:
            ww=lambda s:(s/ref)**(-.02 if constant else op.delta(rh(s),A,alpha))
            integ=4*np.pi*quad(lambda s:s*s*ww(s)*rh(s),0,r,epsabs=0,epsrel=1e-11)[0]
            vals.append(G*integ/(r*r*ww(r)))
        refs['constant_delta' if constant else 'variable_delta']={'fields':vals,'relative_span':max(vals)/min(vals)-1}
    assert refs['constant_delta']['relative_span']<1e-12
    assert refs['variable_delta']['relative_span']>1e-5
    out['checks']['physical_reference_length']=refs
    # Integral form differentiates to the expanded operator including delta-prime log.
    def field(s):
        integ=4*np.pi*quad(lambda t:t*t*w(t)*rh(t),0,s,epsabs=0,epsrel=1e-11)[0]
        return G*integ/(s*s*w(s))
    h=r*1e-5;gp=(field(r+h)-field(r-h))/(2*h)
    rhop=-2*rhoc*r/L**2/(1+(r/L)**2)**2
    dp=op.ddelta_drho(rh(r),A,alpha)*rhop
    residual=gp+((2+op.delta(rh(r),A,alpha))/r+dp*np.log(r/rstar))*g-4*np.pi*G*rh(r)
    rel=abs(residual/(4*np.pi*G*rh(r)));assert rel<1e-8
    out['checks']['quadrature_operator_relative_residual']=rel
    # A weak metric alone is insufficient if energy and rest-mass densities differ.
    ns['C']=C*100;e=EOS(c=C*100,energy_ratio=1.2);P=e.K*rh(r)**2
    d=capture(e,P,A,alpha,'emb1',rstar)
    grad=d['rhs'](r,[P,1.2*mass,g,gN])[0]
    ratio=grad/(-rh(r)*g);assert abs(ratio-1.44)<1e-8
    out['checks']['noncold_source_counterexample']={'epsilon_over_rhob_c2':1.2,'gradient_over_minus_rhob_gop':ratio}
    ns['C']=C
    # CGS->SI: length *.01, mass *.001, pressure/energy density *.1, rho *1000.
    unit={}
    for emb in ['emb1','emb2']:
        P=1e33;r=3e5;e=EOS();rho=e.rho_of_P(P);m=4*np.pi*r**3*rho/3
        dc=op.delta(rho,A,alpha);w0=(r/1e6)**dc;aN=G*m/r**2
        state=np.array([P,m,aN*1.01,aN]) if emb=='emb1' else np.array([P,m,m*w0*1.01])
        cgs=np.array(capture(e,P,A,alpha,emb)['rhs'](r,state))
        class SI:
            P_min=e.P_min*.1
            def rho_of_P(self,p):return e.rho_of_P(p/.1)*1000
            def eps_of_P(self,p):return e.eps_of_P(p/.1)*.1
        # Helpers retain their production globals, so transform their dimensional constants too.
        ns['G']=G*1e-3;ns['C']=C*.01;op.RHO0_RDT=150000.
        factors=np.array([.1,.001,.01,.01]) if emb=='emb1' else np.array([.1,.001,.001])
        si=np.array(capture(SI(),P*.1,A,alpha,emb,1e4)['rhs'](r*.01,state*factors))
        err=float(np.max(np.abs(si/(cgs*factors/.01)-1)))
        assert err<1e-10
        unit[emb]=err
        ns['G']=G;ns['C']=C;op.RHO0_RDT=150.
    out['checks']['CGS_SI_RHS_max_relative_errors']=unit
    out['status']='audit_complete_with_documented_central_initialization_failures'
    sources=[SOURCE,ROOT/'paper3_NeutronStars/ns_v7/tov.py',ROOT/'paper3_NeutronStars/ns_v7/eos.py',Path(__file__)]
    out['source_sha256']={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
    dest=ROOT/'2026_june_14_make_it_consistent/step13_evidence/audit.json';dest.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':run()
