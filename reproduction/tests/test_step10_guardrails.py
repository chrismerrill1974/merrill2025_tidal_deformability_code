"""Run separately with RDT_TEST_KIND=solar or wd to isolate legacy module names."""
import csv
import importlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
KIND=os.environ['RDT_TEST_KIND']
FOLDER=ROOT/('paper1_SolarCore/solarcore_v3' if KIND=='solar' else 'paper2_WhiteDwarfs/wd_v5')
sys.path.insert(0,str(FOLDER));sys.path.insert(0,str(ROOT))
import model
import validate
import numerical_guardrails as guards
producer=importlib.import_module('scan' if KIND=='solar' else 'predict')


class FailureHandling(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.path=Path(self.temp.name)/'data.csv'
        self.params=model.Params(.0042,.02)

    def art(self):
        if KIND=='solar': return producer.artifact(self.path,[0,.0042],[.02])
        return producer.artifact(self.path,[(.02,.0042)],[1e5])

    def row(self):
        if KIND=='solar':
            ratios={b:1. for b in model.BETA}
            return dict(A=0.,alpha=.02,xi1=6.89,rho_c=76.,Tc_ratio=1.,
                        R_pp=1.,R_pep=1.,R_be7=1.,R_b8=1.,
                        chi2=model.chi2_fluxes(ratios),max_abs_dcs=0.)
        return dict(alpha=.02,A=.0042,rho_c=1e5,M_std=.15,R_std_km=16000.,
                    dlnR=-.002,dlnM=-.003,dlnR_fixedM=-.003)

    def run_producer(self):
        if KIND=='solar':return producer.main(self.path,[0,.0042],[.02])
        return producer.main(self.path,[(.02,.0042)],[1e5])

    def test_solver_failure_even_with_event(self):
        sol=SimpleNamespace(success=False,message='injected failure',status=1,
                            t_events=[np.array([3.])])
        with patch.object(model,'solve_ivp',return_value=sol),self.assertRaisesRegex(RuntimeError,'ODE solver failed'):
            model.lane_emden(self.params,1e5,150.)

    def test_missing_surface(self):
        sol=SimpleNamespace(success=True,status=0,t_events=[np.array([])])
        with patch.object(model,'solve_ivp',return_value=sol),self.assertRaisesRegex(RuntimeError,'surface'):
            model.lane_emden(self.params,1e5,150.)

    def test_nonfinite_solution(self):
        sol=SimpleNamespace(success=True,status=1,t_events=[np.array([3.])],
            t=np.array([1e-6,3.]),y=np.array([[1.,0.],[0.,np.nan]]),y_events=[np.array([[0.,-.2]])])
        with patch.object(model,'solve_ivp',return_value=sol),self.assertRaisesRegex(RuntimeError,'Nonfinite'):
            model.lane_emden(self.params,1e5,150.)

    def test_invalid_inputs(self):
        with self.assertRaises(ValueError):model.lane_emden(model.Params(-.1,.02),1e5,150.)
        with self.assertRaises(ValueError):model.lane_emden(self.params,float('nan'),150.)

    def test_gate_false_or_exception_raises(self):
        with self.assertRaises(RuntimeError):guards.run_gates([lambda:False])
        def broken():raise ValueError('injected')
        with self.assertRaises(RuntimeError):guards.run_gates([broken])

    def test_validator_cli_failure_nonzero(self):
        function='calibrate' if KIND=='solar' else 'build_wd'
        code=f'import model,runpy; model.{function}=lambda *a,**k: (_ for _ in ()).throw(RuntimeError("injected")); runpy.run_path("validate.py",run_name="__main__")'
        run=subprocess.run([sys.executable,'-c',code],cwd=FOLDER,capture_output=True,text=True)
        self.assertNotEqual(run.returncode,0);self.assertIn('Validation failed',run.stderr)

    def test_producer_gate_blocks_output_and_records_failure(self):
        solve='calibrate' if KIND=='solar' else 'mr_pair'
        with patch.object(validate,'run_all',return_value=False),patch.object(producer,solve) as mocked:
            with self.assertRaises(RuntimeError):self.run_producer()
            mocked.assert_not_called()
        art=self.art();art.load()
        self.assertEqual(art.rows,[]);self.assertEqual(art.state,'failed')
        self.assertIn('validation',art.log.read_text())
        with self.assertRaises(ValueError):art.load(complete=True)

    def test_point_failure_stays_visible_and_not_complete(self):
        with patch.object(validate,'run_all',return_value=True):
            if KIND=='solar':
                base=SimpleNamespace(xi1=6.89,rho_c=76.)
                with patch.object(producer,'calibrate',side_effect=[base,RuntimeError('point failed')]):
                    with self.assertRaises(RuntimeError):self.run_producer()
            else:
                with patch.object(producer,'mr_pair',side_effect=RuntimeError('point failed')):
                    with self.assertRaises(RuntimeError):self.run_producer()
        art=self.art();art.load();self.assertEqual(art.state,'failed')
        self.assertIn('0.0042',art.log.read_text())
        with self.assertRaises(ValueError):art.load(complete=True)

    def test_legacy_csv_refused_without_mutation(self):
        self.path.write_text('historical data\n');before=self.path.read_bytes()
        with self.assertRaisesRegex(ValueError,'Unmanifested'):self.art().load()
        self.assertEqual(self.path.read_bytes(),before)
        self.assertFalse(self.art().meta.exists())

    def test_tampered_csv_rejected(self):
        art=self.art();art.append(self.row());self.path.write_text(self.path.read_text()+'broken\n')
        with self.assertRaisesRegex(ValueError,'checksum'):self.art().load()

    def test_incompatible_settings_refused(self):
        art=self.art();art.append(self.row())
        other=self.art();other.contract['settings']['theta_floor']=.9
        with self.assertRaisesRegex(ValueError,'provenance'):other.load()

    def test_nonfinite_row_not_appended(self):
        art=self.art();row=self.row();row['rho_c']=float('nan')
        with self.assertRaises(ValueError):art.append(row)
        self.assertFalse(self.path.exists())

    def test_row_algebra_checked_even_with_new_checksum(self):
        art=self.art();row=self.row()
        row['chi2' if KIND=='solar' else 'dlnR_fixedM']+=.1
        art.rows=[row];art.save('partial')
        with self.assertRaises(ValueError):self.art().load()

    def test_duplicate_and_outside_grid_rejected(self):
        art=self.art();row=self.row();art.append(row)
        with self.assertRaises(ValueError):art.append(row)
        row['alpha']=.03
        with self.assertRaises(ValueError):art.append(row)
        art.rows=[self.row(),self.row()];art.save('partial')
        with self.assertRaisesRegex(ValueError,'Duplicate'):self.art().load()

    def test_partial_cannot_be_analysis_input(self):
        art=self.art();art.append(self.row())
        with self.assertRaises(ValueError):self.art().load(complete=True)

    def test_complete_roundtrip_and_resume(self):
        art=self.art();art.append(self.row())
        if KIND=='solar':
            r=self.row();r['A']=.0042;art.append(r)
        art.save('complete');self.assertEqual(len(self.art().load(complete=True)),len(art.expected))
        solve='calibrate' if KIND=='solar' else 'mr_pair'
        with patch.object(validate,'run_all',return_value=True),patch.object(producer,solve) as solver:
            if KIND=='solar':solver.return_value=SimpleNamespace(xi1=6.89,rho_c=76.)
            self.run_producer()
            self.assertEqual(solver.call_count,1 if KIND=='solar' else 0)

    @unittest.skipUnless(KIND=='solar','solar calibration only')
    def test_calibration_nonconvergence_and_reference_validation(self):
        with self.assertRaisesRegex(RuntimeError,'did not converge'):
            model.calibrate(self.params,max_iter=1)
        with self.assertRaises(ValueError):model.calibrate(self.params,r0_mode='typo')
        with self.assertRaises(ValueError):model.calibrate(self.params,max_iter=0)

    @unittest.skipUnless(KIND=='solar','solar reference iteration only')
    def test_reference_iteration_must_converge_even_with_fixed_density(self):
        M,R,_=model.units()
        calls=[]
        def profile(*args,**kwargs):
            endpoint=6. if len(calls)%2==0 else 7.
            calls.append(endpoint)
            return np.array([1e-6,endpoint/2,endpoint]), np.array([1.,.5,0.]), np.array([0.,-.1,-.1])
        def integral(*args,**kwargs):return M/(4*np.pi*(R/calls[-1])**3*80.)
        with patch.object(model,'lane_emden',side_effect=profile),patch.object(model,'simpson',side_effect=integral):
            with self.assertRaisesRegex(RuntimeError,'did not converge'):
                model.calibrate(self.params,r0_mode='Rsun',max_iter=4)
        self.assertEqual(len(calls),4)

    @unittest.skipUnless(KIND=='solar','solar analysis only')
    def test_analysis_refuses_archive_before_figures(self):
        import analyze
        self.path.write_text('historical\n');out=Path(self.temp.name)/'figures'
        with self.assertRaises(ValueError):analyze.main(self.path,str(out))
        self.assertFalse(out.exists())


if __name__=='__main__':unittest.main(verbosity=2)
