"""Numerical comparisons at the locally refined primary alpha peak."""
from pathlib import Path
import json
from dataclasses import replace,asdict
from scipy.optimize import minimize_scalar,brentq
from solver import Settings,calibrate,response
from inference import score
from numerical_guardrails import fingerprint,sha
HERE=Path(__file__).resolve().parent
if __name__=='__main__':
    source=HERE/'results/alpha_peak_checks.json';out=HERE/'results/peak_numerics.json'
    if out.exists():raise ValueError('New artifact required')
    p=json.loads(source.read_text())['checks']['source_precision_diagonal_proxy_v1']
    assert p['status']=='complete'
    alpha=p['alpha'];base=p['endpoint']
    settings={'tight':replace(Settings(),rtol=2e-13,atol=2e-15,max_step=.02),'profile':replace(Settings(),profile_n=3201),'floor_1e-6':replace(Settings(),floor=1e-6),'floor_1e-10':replace(Settings(),floor=1e-10)}
    result={'contract':fingerprint(HERE,['solver.py','inference.py','peak_numerics.py'],{'settings':{k:asdict(v) for k,v in settings.items()},'input_sha256':sha(source)}),'alpha':alpha,'baseline_endpoint':base,'checks':{}}
    for name,s in settings.items():
        try:
            std=calibrate(0,.1,s)
            f=lambda A:score(response(A,alpha,s,std),'source_precision_diagonal_proxy_v1')
            opt=minimize_scalar(f,bounds=(0,.004),method='bounded',options={'xatol':2e-10})
            if not opt.success:raise RuntimeError('Minimum failed')
            qmin=min(f(0),float(opt.fun))
            end=float(brentq(lambda A:f(A)-qmin-2.71,base*.98,base*1.02,xtol=2e-10))
            result['checks'][name]={'status':'complete','endpoint':end,'shift':end-base}
        except Exception as e:result['checks'][name]={'status':'failed','error':str(e)}
    out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result['checks'],indent=2))
