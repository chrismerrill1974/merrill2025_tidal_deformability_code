"""Two local alpha refinements of observed NS deficit extrema.
New solar endpoints are recomputed from the source-precision score, not obtained
by interpolating an accepted contour. Finite refinement is not a global proof.
"""
from pathlib import Path
import argparse,concurrent.futures,json,sys,traceback
import numpy as np
from scan import ROOT,HERE,sha,metrics,REFERENCES,EMBEDDINGS,VARIANT,inputs
sys.path.insert(0,str(ROOT/'paper1_SolarCore/solar_step11'))
from solver import Settings,calibrate,response
from inference import intervals,score

def load(folder):
    m=json.loads((folder/'manifest.json').read_text())
    if m['status']!='complete':raise ValueError('Incomplete parent scan')
    for p,h in m['contract']['source_sha256'].items():
        if sha(ROOT/p)!=h:raise ValueError('source mismatch '+p)
    for p,h in m['sha256'].items():
        if sha(folder/p)!=h:raise ValueError('artifact mismatch '+p)
    rows=[json.loads(v) for v in (folder/'evaluations.jsonl').read_text().splitlines()]
    return m,rows,json.loads((folder/'baseline.json').read_text())['primary']

def run(parent,out,workers):
    out.mkdir(parents=True,exist_ok=False);inputs();pm,rows,baseline=load(parent)
    sources=[HERE/'refine.py',HERE/'scan.py',ROOT/'paper1_SolarCore/solar_step11/solver.py',ROOT/'paper1_SolarCore/solar_step11/inference.py',ROOT/'paper3_NeutronStars/ns_step14/model.py']
    manifest={'status':'running','parent_manifest_sha256':sha(parent/'manifest.json'),'source_sha256':{str(p.relative_to(ROOT)):sha(p) for p in sources},'rounds':2,'solar_settings':Settings().__dict__,'solar_A_mesh':np.linspace(0,.04,17).tolist(),'alpha_refinement':'geometric midpoints on each side of observed minima','workers':workers}
    def save():(out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    save();std=calibrate(0,.1);solar_rows=[];endpoints=[];newrows=[];rounds=[]
    journal=(out/'evaluations.jsonl').open('w')
    try:
        for iteration in [1,2]:
            alphas=sorted(set(r['point']['alpha'] for r in rows));candidates=set();fractions={1.}
            before={}
            for ref in REFERENCES:
                for emb in EMBEDDINGS:
                    vals=[r for r in rows if r['point']['rstar_cm']==ref and r['point']['embedding']==emb]
                    for q in ['delta_mass','delta_R14']:
                        best=min(vals,key=lambda r:r[q]);a=best['point']['alpha'];i=alphas.index(a)
                        before[f'{ref:g}/{emb}/{q}']={'alpha':a,'A':best['point']['A'],'value':best[q]}
                        fractions.add(best['point']['fraction'])
                        if i>0:candidates.add(float(np.sqrt(a*alphas[i-1])))
                        if i<len(alphas)-1:candidates.add(float(np.sqrt(a*alphas[i+1])))
            candidates=sorted(candidates-set(alphas));jobs=[]
            for alpha in candidates:
                cache={}
                def get(A):
                    A=float(A)
                    if A not in cache:
                        v=response(A,alpha,std=std);cache[A]=score(v,VARIANT);solar_rows.append({'round':iteration,'A':A,'alpha':alpha,'response':v,'score':cache[A]})
                    return cache[A]
                mesh=np.linspace(0,.04,17);v=intervals(mesh,[get(A) for A in mesh],evaluate=get)
                if v['upper_endpoint'] is None or len(v['intervals'])!=1 or v['intervals'][0]['lower']!=0:raise RuntimeError('refined solar topology changed')
                ep={'alpha':alpha,'boundary_A':v['upper_endpoint'],'origin':f'NS_alpha_refinement_{iteration}','solar_slice':v};endpoints.append(ep)
                for ref in REFERENCES:
                    for emb in EMBEDDINGS:
                        for fraction in sorted(fractions):
                            jobs.append({k:ep[k] for k in ['alpha','boundary_A','origin']}|{'A':v['upper_endpoint']*fraction,'fraction':fraction,'rstar_cm':ref,'embedding':emb})
                print('solar endpoint',iteration,alpha,v['upper_endpoint'],flush=True)
                (out/'solar_evaluations.partial.json').write_text(json.dumps(solar_rows)+'\n')
            with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as pool:
                for p,value in zip(jobs,pool.map(metrics,jobs,chunksize=1)):
                    value['delta_mass']=value['domain_peak']['M']/baseline['domain_peak']['M']-1
                    value['delta_R14']=value['R14']['R']/baseline['R14']['R']-1
                    value['reused_GR']=False;newrows.append(value);rows.append(value)
                    journal.write(json.dumps(value)+'\n');journal.flush()
            rounds.append({'round':iteration,'before':before,'new_alphas':candidates,'fractions':sorted(fractions),'new_rows':len(jobs)})
            print('round',iteration,'complete',len(newrows),'new NS rows',flush=True)
        journal.close()
        for name,value in [('solar_endpoints.json',endpoints),('solar_evaluations.json',solar_rows),('rounds.json',rounds)]: (out/name).write_text(json.dumps(value,indent=2)+'\n')
        (out/'solar_evaluations.partial.json').unlink(missing_ok=True)
        manifest.update(status='complete',rows=len(newrows),solar_evaluations=len(solar_rows),sha256={name:sha(out/name) for name in ['evaluations.jsonl','solar_endpoints.json','solar_evaluations.json','rounds.json']});save()
    except Exception as exc:
        journal.close();manifest.update(status='failed',error=repr(exc),traceback=traceback.format_exc());save();raise
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--scan',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--workers',type=int,default=4);a=p.parse_args();run(a.scan,a.output,a.workers)
