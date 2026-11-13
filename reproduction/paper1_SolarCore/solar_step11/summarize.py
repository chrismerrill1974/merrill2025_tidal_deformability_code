"""Read only complete, checksum-verified artifacts; summarize numerical scope."""
from pathlib import Path
import json,sys
import numpy as np
from inference import VARIANTS,CONTRACT,INPUT,ROOT
sys.path.insert(0,str(ROOT))
from numerical_guardrails import sha
HERE=Path(__file__).resolve().parent

def load_region(folder):
    p=Path(folder);m=json.loads((p/'manifest.json').read_text())
    if m['status']!='complete':raise ValueError('Refuse incomplete/failed region')
    for f,h in m['sha256'].items():
        if sha(p/f)!=h:raise ValueError('Checksum mismatch')
    for f,h in m['contract']['source_sha256'].items():
        if sha(HERE/f)!=h:raise ValueError('Incompatible source')
    settings=m['contract']['settings']
    if sha(CONTRACT)!=settings['score_contract_sha256'] or sha(INPUT)!=settings['input_sha256']:raise ValueError('Incompatible input contract')
    rows=[json.loads(line) for line in (p/'evaluations.jsonl').read_text().splitlines()]
    keys={(r['A'],r['alpha']) for r in rows}
    if len(keys)!=len(rows) or len(rows)!=m['evaluations']:raise ValueError('Duplicate/count mismatch')
    if not {(a,l) for a in settings['A'] for l in settings['alpha']}<=keys:raise ValueError('Incomplete declared grid')
    if any(r['calibration_residual']>=settings['solver']['calibration_tol'] or not all(np.isfinite(v) for v in r.values()) for r in rows):raise ValueError('Invalid response row')
    region=json.loads((p/'region.json').read_text())
    if [s['alpha'] for s in region['slices']]!=settings['alpha']:raise ValueError('Incomplete slices')
    return region

def summarize():
    a=load_region(HERE/'results/solar_a');r=load_region(HERE/'results/solar_R')
    audit=json.loads((HERE/'results/numerical_audit/audit.json').read_text())
    report={'envelopes':a['envelopes'],'reference_R_envelopes':r['envelopes'],'mesh_errors':{},'numerical_endpoint_checks':{},'response_differences':{},'sampled_joint_minima':{}}
    for variant in VARIANTS:
        vals=[sl['variants'][variant] for sl in a['slices']]
        best=min(a['slices'],key=lambda sl:sl['variants'][variant]['minimum_Q'])
        report['sampled_joint_minima'][variant]={'alpha':best['alpha'],'A':best['variants'][variant]['minimum_A'],'Q':best['variants'][variant]['minimum_Q'],'interpretation':'minimum over refined A slices at sampled alpha; no significance claim'}
        report['mesh_errors'][variant]={'coarse_endpoint_max_abs_difference':max([abs(v['upper_endpoint']-v['coarse_mesh_endpoint']) for v in vals if v['upper_endpoint'] is not None and v['coarse_mesh_endpoint'] is not None],default=None),'fine_endpoint_max_abs_difference':max([abs(v['upper_endpoint']-v['fine_mesh_endpoint']) for v in vals if v['upper_endpoint'] is not None and v['fine_mesh_endpoint'] is not None],default=None),'component_count_disagreements':sum(len(set(v['mesh_component_counts']))>1 for v in vals),'max_root_score_residual':max(v['maximum_root_score_residual'] for v in vals),'minima_A_range':[min(v['minimum_A'] for v in vals),max(v['minimum_A'] for v in vals)]}
    for key in ['tight','start','profile','floor_1e-6','floor_1e-10','si']:
        vals=[];failed=[]
        for row in audit['endpoint_checks']:
            t=row['settings'][key];b=row['settings']['baseline']
            if t['status']=='complete' and b['status']=='complete':vals.append(abs(t['endpoint']-b['endpoint']))
            else:failed.append({'alpha':row['alpha'],'test':t,'baseline':b})
        report['numerical_endpoint_checks'][key]={'successful_pairs':len(vals),'maximum_abs_endpoint_shift':max(vals,default=None),'failures':failed}
    for key in ['tight','start','profile','floor_1e-6','floor_1e-10','si','loose']:
        vals=[];failed=[]
        for row in audit['points']:
            b=row['settings']['baseline'];t=row['settings'][key]
            if b['status']=='complete' and t['status']=='complete':vals.append({k:abs(t['response'][k]-b['response'][k]) for k in ['R_pp','R_pep','R_be7','R_b8','h']})
            else:failed.append({'A':row['A'],'alpha':row['alpha'],'test':t})
        report['response_differences'][key]={'successful_pairs':len(vals),'max_abs':{k:max(v[k] for v in vals) for k in vals[0]} if vals else {},'failures':failed}
    (HERE/'results/summary.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))
if __name__=='__main__':summarize()
