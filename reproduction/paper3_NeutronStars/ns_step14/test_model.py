"""Meaningful failure and formulation checks for the new versioned NS path."""
import unittest
from unittest.mock import patch
from types import SimpleNamespace
import numpy as np
import model
from model import EOS,solve
from direct import solve_direct
class Checks(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.e=EOS()
    def test_rejects_invalid_parameters(self):
        for kw in [dict(A=-.01),dict(A=np.nan),dict(alpha=0),dict(r0=0),dict(r0=1e6),dict(rstar=0),dict(rtol=0),dict(atol=0),dict(max_step=0),dict(surface_factor=.5),dict(embedding='unknown')]:
            with self.subTest(kw=kw),self.assertRaises(ValueError):solve(self.e,2e35,**kw)
    def test_rejects_EOS_extrapolation(self):
        for p in [self.e.P_min*.9,self.e.P_max*1.01,np.nan]:
            with self.subTest(p=p),self.assertRaises(ValueError):self.e.values(p)
        with self.assertRaises(ValueError):solve(self.e,self.e.P_max*1.01)
    def test_solver_failure_is_not_output(self):
        with patch.object(model,'solve_ivp',return_value=SimpleNamespace(success=False,message='injected step failure')):
            with self.assertRaisesRegex(RuntimeError,'injected'):solve(self.e,2e35)
    def test_partial_endpoint_is_not_output(self):
        fake=SimpleNamespace(success=True,status=0,t=np.array([0.,1.]),y=np.ones((5,2)),message='partial')
        with patch.object(model,'solve_ivp',return_value=fake):
            with self.assertRaises(RuntimeError):solve(self.e,2e35)
    def test_nonfinite_endpoint_is_not_output(self):
        def fake(rhs,span,y0,**kw):return SimpleNamespace(success=True,status=0,t=np.array(span),y=np.full((5,2),np.nan),message='nan')
        with patch.object(model,'solve_ivp',side_effect=fake):
            with self.assertRaises(RuntimeError):solve(self.e,2e35)
    def test_wrong_success_status_is_not_output(self):
        with patch.object(model,'solve_ivp',return_value=SimpleNamespace(success=True,status=1,message='unexpected termination')):
            with self.assertRaises(RuntimeError):solve(self.e,2e35)
    def test_missing_radial_surface_is_failure(self):
        with self.assertRaisesRegex(RuntimeError,'surface'):solve_direct(self.e,2e35,rmax=110.)
    def test_correct_central_series_is_used(self):
        # Inspect initial states supplied to the integrator; pressure start and weighted masses.
        for emb in ['emb1','emb2']:
            captured={}
            def fake(rhs,span,y0,**kw):
                captured.update(span=span,y0=y0);raise RuntimeError('captured')
            with patch.object(model,'solve_ivp',side_effect=fake),self.assertRaisesRegex(RuntimeError,'captured'):
                solve(self.e,2e35,.0087,.12,emb)
            ec,rho=self.e.values(2e35);d=model.delta(rho,.0087,.12);r0=100.
            B=2*np.pi*model.G/model.C**4*(ec+2e35)*(ec/3+2e35)*3/(3+d) if emb=='emb1' else 2*np.pi*model.G/model.C**4*(ec+2e35)*(ec/(3+d)+2e35)
            self.assertAlmostEqual(captured['span'][0]/(-np.log1p(-B*r0*r0/2e35)),1.,places=14)
            h=(r0/1e6)**d
            self.assertAlmostEqual(captured['y0'][4]/captured['y0'][1]/h,3/(3+d),places=14)
    def test_GR_both_modified_paths(self):
        a=solve(self.e,2e35,embedding='emb1');b=solve(self.e,2e35,embedding='emb2')
        for k in ['M','R']:self.assertLess(abs(a[k]/b[k]-1),1e-9)
    def test_uniform_density_GR_analytic_solution(self):
        rho=1e15;eps=rho*model.C**2;R=8e5
        k=8*np.pi*model.G*rho/(3*model.C**2);s=np.sqrt(1-k*R*R)
        pc=eps*(1-s)/(3*s-1)
        class Uniform:
            P_min=pc*1e-8;P_max=pc*2
            def values(self,P):return eps,rho
        e=Uniform();q=e.P_min/eps;rs=np.sqrt((1-(s*(1+3*q)/(1+q))**2)/k)
        mass=4*np.pi*rho*rs**3/3/model.MSUN
        for emb in ['emb1','emb2']:
            a=solve(e,pc,embedding=emb,surface_factor=1.,rtol=2e-11,atol=1e-19,max_step=.04)
            self.assertLess(abs(a['R']/(rs/1e5)-1),1e-8)
            self.assertLess(abs(a['M']/mass-1),1e-8)
    def test_direct_nonzero_matches_flux(self):
        for emb in ['emb1','emb2']:
            a=solve(self.e,2e35,.0087,.12,emb);b=solve_direct(self.e,2e35,.0087,.12,emb)
            for k in ['M','R']:self.assertLess(abs(a[k]/b[k]-1),5e-7)
if __name__=='__main__':unittest.main()
