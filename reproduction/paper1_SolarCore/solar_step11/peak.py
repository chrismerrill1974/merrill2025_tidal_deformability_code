"""Local alpha refinement around each observed sampled-envelope maximum.
This does not certify a global continuous-alpha bound.
"""
from pathlib import Path
import json
import numpy as np
from scipy.optimize import minimize_scalar,brentq
from solver import Settings,response,calibrate
from inference import score,VARIANTS
from summarize import load_region
from numerical_guardrails import fingerprint,sha
HERE=Path(__file__).resolve().parent
if __name__=='__main__':
    source=HERE/'results/solar_a';region=load_region(source);out=HERE/'results/alpha_peak_checks.json'
    if out.exists():raise ValueError('New output required')
    std=calibrate(0,.1);cache={};evaluations=[]
    def get(A,alpha):
        key=(float(A),float(alpha))
        if key not in cache:
            cache[key]=response(*key,std=std);evaluations.append(cache[key])
        return cache[key]
    record={'contract':fingerprint(HERE,['solver.py','inference.py','peak.py'],{'solver':Settings().__dict__,'region_sha256':sha(source/'region.json'),'alpha_xtol':2e-6}),'checks':{}}
    for variant in VARIANTS:
        e=region['envelopes'][variant]
        if e['grid_limited_slices']:
            record['checks'][variant]={'status':'grid_limited','reason':'No finite overall upper envelope'};continue
        _,alpha=e['sampled_maximum'];alphas=[r['alpha'] for r in region['slices']];i=alphas.index(alpha)
        if i in (0,len(alphas)-1):
            record['checks'][variant]={'status':'boundary_sample_maximum','sampled':e['sampled_maximum']};continue
        try:
            def endpoint(al):
                f=lambda A:score(get(A,al),variant)
                opt=minimize_scalar(f,bounds=(0,.004),method='bounded',options={'xatol':2e-10})
                if not opt.success:raise RuntimeError('Minimum failed')
                qm=min(f(0),float(opt.fun))
                return float(brentq(lambda A:f(A)-qm-2.71,.002,.03,xtol=2e-10))
            fit=minimize_scalar(lambda al:-endpoint(al),bounds=(alphas[i-1],alphas[i+1]),method='bounded',options={'xatol':2e-6})
            if not fit.success:raise RuntimeError('Alpha maximum refinement failed')
            record['checks'][variant]={'status':'complete','alpha':float(fit.x),'endpoint':float(-fit.fun),'sampled':e['sampled_maximum'],'increase_over_sample':float(-fit.fun-e['sampled_maximum'][0]),'bracket':[alphas[i-1],alphas[i+1]]}
        except Exception as exc:
            record['checks'][variant]={'status':'failed','error':str(exc)}
        print(variant,record['checks'][variant],flush=True)
    record['evaluations']=evaluations;out.write_text(json.dumps(record,indent=2)+'\n')
