import sys,unittest,json,tempfile,subprocess
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'paper1_SolarCore/solar_step11'))
from inference import score,intervals,C,ROWS
from solver import response,Settings,calibrate
from dataclasses import replace

class InferenceTests(unittest.TestCase):
    def test_standalone_summary_import(self):
        folder=str(Path(__file__).resolve().parents[1]/'paper1_SolarCore/solar_step11')
        with tempfile.TemporaryDirectory() as d:
            result=subprocess.run([sys.executable,'-c',f'import sys; sys.path.insert(0,{folder!r}); import summarize'],cwd=d,capture_output=True,text=True)
        self.assertEqual(result.returncode,0,result.stderr)

    def test_incomplete_region_refused(self):
        from summarize import load_region
        with tempfile.TemporaryDirectory() as d:
            (Path(d)/'manifest.json').write_text(json.dumps({'status':'failed'}))
            with self.assertRaisesRegex(ValueError,'incomplete/failed'):load_region(d)
    def test_corrupt_region_refused(self):
        from summarize import load_region
        with tempfile.TemporaryDirectory() as d:
            (Path(d)/'region.json').write_text('{}')
            (Path(d)/'manifest.json').write_text(json.dumps({'status':'complete','sha256':{'region.json':'bad'}}))
            with self.assertRaisesRegex(ValueError,'Checksum'):load_region(d)

    def test_disconnected(self):
        x=np.linspace(0,4,81);f=lambda a:10*(a-1)**2*(a-3)**2
        r=intervals(x,[f(a) for a in x],threshold=.1,evaluate=f)
        self.assertEqual(len(r['intervals']),2)
        self.assertGreater(r['upper_endpoint'],3)
        self.assertLess(r['maximum_root_score_residual'],1e-7)
    def test_grid_edge(self):
        x=np.linspace(0,1,11);r=intervals(x,x*x,threshold=2.71)
        json.dumps(r)
        self.assertIsNone(r['upper_endpoint']);self.assertTrue(r['intervals'][-1]['upper_grid_limited'])
    def test_boundary_minimum(self):
        x=np.linspace(0,.04,81);f=lambda a:1e5*(a-.0001)**2
        r=intervals(x,[f(a) for a in x],evaluate=f)
        self.assertAlmostEqual(r['minimum_A'],.0001,places=8)
    def test_source_score_and_legacy(self):
        row=dict(R_pp=1.,R_pep=1.,R_be7=1.,R_b8=1.,h=0.)
        expected=sum((1-C['y'][b])**2/C['variance'][b] for b in C['y'])
        self.assertAlmostEqual(score(row,'source_precision_diagonal_proxy_v1'),expected)
        self.assertAlmostEqual(score(row,'legacy_diagonal_proxy'),.7791452991452995)
    def test_unit_invariance(self):
        a=response(.0087,.15);b=response(.0087,.15,replace(Settings(),si=True))
        for k in ['R_pp','R_pep','R_be7','R_b8','h','xi1']:self.assertLess(abs(a[k]-b[k]),2e-9)
    def test_nonconvergence(self):
        with self.assertRaisesRegex(RuntimeError,'Calibration failed'):calibrate(.0042,.02,replace(Settings(),max_iter=1))
    def test_standard_polytrope(self):
        standard=calibrate(0,.1)
        self.assertAlmostEqual(standard['xi'],6.89684861937696,places=8)
        self.assertAlmostEqual(standard['sol'].y[2,-1],2.01823595,places=7)
        self.assertLess(standard['residual'],1e-12)
    def test_null(self):
        a=response(0,.02)
        self.assertEqual(a['h'],0.)
        for k in ['R_pp','R_pep','R_be7','R_b8']:self.assertEqual(a[k],1.)
if __name__=='__main__':unittest.main(verbosity=2)
