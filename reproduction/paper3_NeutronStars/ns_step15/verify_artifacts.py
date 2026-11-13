"""Verify retained step 15 artifact hashes and required completion evidence."""
from pathlib import Path
import json
from resource_control import require_controls

def main():
    require_controls()
    from scan import ROOT,HERE,sha,inputs
    inputs()
    m=json.loads((ROOT/'2026_june_14_make_it_consistent/STEP15_ARTIFACTS.json').read_text())
    for p,h in m['sha256'].items():
        if sha(ROOT/p)!=h:raise ValueError('Artifact mismatch: '+p)
    for folder in ['resumed95_results','refinement','analysis','checks','figures']:
        d=HERE/folder;j=json.loads((d/'manifest.json').read_text())
        assert j['status']=='complete',folder
        for p,h in j['sha256'].items():assert sha(d/p)==h,(folder,p)
    for name in ['probe95','resume95','refinement','analysis','checks','figures']:
        r=json.loads((HERE/'resource_evidence'/f'{name}.json').read_text())
        assert r['status']=='complete' and r['nice']==19 and r['worker_processes']==1,name
        assert r['cpu_fraction_target']==.95 and all(v=='1' for v in r['native_thread_environment'].values()),name
    print(json.dumps({'status':'pass','artifact_hashes':len(m['sha256']),'complete_scientific_stages':5,'resource_logs_checked':6,'prior_solar_and_NS_sources':'unchanged'}))
if __name__=='__main__':main()
