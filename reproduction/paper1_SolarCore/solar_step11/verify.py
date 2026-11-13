"""Bind each reported minimum/crossing to its saved response and score inputs."""
from pathlib import Path
import json
from inference import score,VARIANTS
from summarize import load_region
HERE=Path(__file__).resolve().parent
if __name__=='__main__':
    output={}
    for name in ['solar_a','solar_R']:
        folder=HERE/'results'/name;region=load_region(folder)
        rows=[json.loads(l) for l in (folder/'evaluations.jsonl').read_text().splitlines()]
        cache={(r['A'],r['alpha']):r for r in rows}
        m=json.loads((folder/'manifest.json').read_text());grid=m['contract']['settings']['A']
        errors=[];min_errors=[];membership_errors=[]
        for sl in region['slices']:
            alpha=sl['alpha']
            for variant,v in sl['variants'].items():
                minimum=score(cache[(v['minimum_A'],alpha)],variant)
                min_errors.append(abs(minimum-v['minimum_Q']))
                for component in v['intervals']:
                    for end in [component['lower'],component['upper']]:
                        if end in [grid[0],grid[-1]]:continue
                        q=score(cache[(end,alpha)],variant)
                        errors.append(abs(q-v['minimum_Q']-2.71))
                for A in grid:
                    q=score(cache[(A,alpha)],variant)-v['minimum_Q']
                    accepted=any(c['lower']<=A<=c['upper'] for c in v['intervals'])
                    if (q<=2.71)!=accepted:membership_errors.append([A,alpha,variant,q])
        assert max(min_errors,default=0)<1e-12
        assert max(errors,default=0)<2e-6
        assert not membership_errors,membership_errors
        output[name]={'slices':len(region['slices']),'evaluations':len(rows),'declared_grid_points':len(grid)*len(region['slices']),'crossings_checked':len(errors),'max_crossing_score_residual':max(errors,default=0),'max_minimum_score_difference':max(min_errors,default=0),'grid_membership_mismatches':0,'max_calibration_residual':max(r['calibration_residual'] for r in rows),'max_calibration_iterations':max(r['iterations'] for r in rows)}
    (HERE/'results/verification.json').write_text(json.dumps(output,indent=2)+'\n');print(json.dumps(output,indent=2))
