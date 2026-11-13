"""Retry endpoint comparisons locally, preserving earlier wide-bracket failures."""
from pathlib import Path
import json
from dataclasses import replace,asdict
from scipy.optimize import minimize_scalar,brentq
from solver import Settings,response,calibrate
from inference import score
from numerical_guardrails import fingerprint,sha
HERE=Path(__file__).resolve().parent
if __name__=='__main__':
    source=HERE/'results/numerical_audit/audit.json';out=HERE/'results/local_endpoint_checks.json'
    if out.exists():raise ValueError('Use a new artifact path')
    settings={'tight':replace(Settings(),rtol=2e-13,atol=2e-15,max_step=.02),'floor_1e-10':replace(Settings(),floor=1e-10)}
    record={'contract':fingerprint(HERE,['solver.py','inference.py','local_checks.py'],{'settings':{k:asdict(v) for k,v in settings.items()},'input_sha256':sha(source),'relative_bracket_half_width':.02}),'checks':[]}
    for row in json.loads(source.read_text())['endpoint_checks']:
        alpha=row['alpha'];base=row['settings']['baseline']['endpoint']
        item={'alpha':alpha,'baseline':base,'comparisons':{}}
        for label,s in settings.items():
            last=None
            try:
                std=calibrate(0,.1,s)
                def f(A):
                    global last
                    last=float(A)
                    return score(response(A,alpha,s,std),'source_precision_diagonal_proxy_v1')
                opt=minimize_scalar(f,bounds=(0,.001),method='bounded',options={'xatol':2e-10})
                if not opt.success:raise RuntimeError('Minimum refinement failed')
                qm=min(f(0),float(opt.fun))
                endpoint=float(brentq(lambda A:f(A)-qm-2.71,base*.98,base*1.02,xtol=2e-10))
                item['comparisons'][label]={'status':'complete','endpoint':endpoint,'shift':endpoint-base}
            except Exception as exc:item['comparisons'][label]={'status':'failed','A':last,'error':str(exc)}
        record['checks'].append(item)
        print(alpha,item['comparisons'],flush=True)
    out.write_text(json.dumps(record,indent=2)+'\n')
