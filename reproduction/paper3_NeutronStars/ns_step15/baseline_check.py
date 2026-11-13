"""Retain the failed default-precision GR gate and its tighter resolution."""
from pathlib import Path
import json
from scan import metrics,sha
HERE=Path(__file__).resolve().parent
if __name__=='__main__':
    result={}
    for name,settings in [('default',dict(rtol=2e-9,atol=1e-17,max_step=.1)),('tight',dict(rtol=2e-11,atol=1e-19,max_step=.04))]:
        a=metrics(dict(A=0.,alpha=.1,embedding='emb1',rstar_cm=1e6,settings=settings));b=metrics(dict(A=0.,alpha=.5,embedding='emb2',rstar_cm=5e5,settings=settings))
        diff={k:b[k][v]/a[k][v]-1 for k,v in [('domain_peak','M'),('R14','R')]}
        result[name]={'primary':a,'crosscheck':b,'relative_differences':diff,'passes_1e8_gate':max(abs(x) for x in diff.values())<=1e-8}
    if not result['tight']['passes_1e8_gate']:raise RuntimeError('Tight GR gate still failed')
    result['script_sha256']=sha(__file__)
    (HERE/'diagnostics/baseline_resolution.json').write_text(json.dumps(result,indent=2)+'\n')
    print({k:v['relative_differences'] for k,v in result.items() if isinstance(v,dict)})
