"""Run a new, fail-closed solar-region artifact; no historical cache admission."""
from pathlib import Path
import sys,json,argparse,hashlib,traceback
from dataclasses import asdict
import numpy as np
from solver import Settings,response,calibrate
from inference import score,intervals,VARIANTS,CONTRACT,INPUT
from numerical_guardrails import fingerprint,sha
HERE=Path(__file__).resolve().parent

def run(output,reference='a',na=81,nalpha=65):
    out=Path(output);out.mkdir(parents=True,exist_ok=False)
    s=Settings(reference=reference)
    contract=fingerprint(HERE,['solver.py','inference.py','run.py'],{'solver':asdict(s),'A':np.linspace(0,.04,na).tolist(),'alpha':np.geomspace(.02,.5,nalpha).tolist(),'score_contract_sha256':sha(CONTRACT),'input_sha256':sha(INPUT),'threshold':2.71,'crossing_xtol':2e-10})
    manifest={'status':'running','contract':contract}
    def save_manifest(): (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    save_manifest();cache={};point=None
    try:
        std=calibrate(0,.1,s)
        def get(A,alpha):
            nonlocal point
            key=(float(A),float(alpha));point={'A':key[0],'alpha':key[1]}
            if key not in cache:
                row=response(*key,s,std);cache[key]=row
                with (out/'evaluations.jsonl').open('a') as f:f.write(json.dumps(row)+'\n')
            return cache[key]
        aa=np.linspace(0,.04,na);alphas=np.geomspace(.02,.5,nalpha)
        result={'score_contract':'source_precision_diagonal_proxy_v1','coverage_claimed':False,'interpretation':'sampled nominal fixed-alpha slice sets; toy response and deterministic helioseismic penalty','reference':reference,'slices':[],'envelopes':{}}
        for alpha in alphas:
            rows=[get(A,alpha) for A in aa]
            slice_result={'alpha':float(alpha),'variants':{}}
            for variant in VARIANTS:
                values=[score(r,variant) for r in rows]
                refined=intervals(aa,values,evaluate=lambda A:score(get(A,alpha),variant))
                coarse=intervals(aa[::2],values[::2]);fine=intervals(aa,values)
                refined['coarse_mesh_endpoint']=coarse['upper_endpoint'];refined['fine_mesh_endpoint']=fine['upper_endpoint']
                refined['mesh_component_counts']=[len(coarse['intervals']),len(fine['intervals']),len(refined['intervals'])]
                slice_result['variants'][variant]=refined
            result['slices'].append(slice_result)
            (out/'region.partial.json').write_text(json.dumps(result,indent=2)+'\n')
            print(reference,'alpha',alpha,'endpoint',slice_result['variants'][VARIANTS[0]]['upper_endpoint'], 'evaluations',len(cache),flush=True)
        for variant in VARIANTS:
            values=[(v['variants'][variant]['upper_endpoint'],v['alpha']) for v in result['slices']]
            finite=[v for v in values if v[0] is not None]
            result['envelopes'][variant]={'sampled_maximum':max(finite) if finite else None,'grid_limited_slices':sum(v[0] is None for v in values),'coarser_alpha_sample_maximum':max([v for v in values[::2] if v[0] is not None],default=None)}
        (out/'region.json').write_text(json.dumps(result,indent=2)+'\n');(out/'region.partial.json').unlink()
        manifest.update(status='complete',evaluations=len(cache),sha256={f:sha(out/f) for f in ['evaluations.jsonl','region.json']});save_manifest()
    except Exception as exc:
        manifest.update(status='failed',failed_point=point,error=str(exc));save_manifest()
        (out/'failure.txt').write_text(traceback.format_exc());raise

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);p.add_argument('--reference',default='a',choices=['a','Rsun']);p.add_argument('--na',type=int,default=81);p.add_argument('--nalpha',type=int,default=65)
    a=p.parse_args();run(a.output,a.reference,a.na,a.nalpha)
