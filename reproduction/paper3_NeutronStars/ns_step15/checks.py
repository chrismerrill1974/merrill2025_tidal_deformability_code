"""Selected extrema controls, independent field roots and saturation diagnostics."""
from pathlib import Path
import argparse,json,sys,traceback
from unittest.mock import patch
import numpy as np
from scipy.optimize import brentq
from scan import ROOT,HERE,REFERENCES,EMBEDDINGS,sha,metrics
from refine import load
from analyze import refined
import model
from model import EOS
from direct import solve_direct
BASELINES={}
def task(job):
    point=job['point'];mode=job['mode'];p=dict(point);settings={'rtol':2e-11,'atol':1e-19,'max_step':.04};interp='linear'
    if mode=='start10':settings['r0']=10.
    if mode=='start1000':settings['r0']=1000.
    if mode=='pchip':interp='pchip'
    if mode=='saturated':settings={'rtol':2e-9,'atol':1e-17,'max_step':.1}
    p['settings']=settings;p['pressure_grid_n']=18
    if mode=='direct':
        eos=EOS();evs=[]
        kw={'rtol':2e-13,'atol':1e-21,'max_step':500.,'rstar':point['rstar_cm']}
        def fixed(A):
            def ev(lp):
                v=solve_direct(eos,float(np.exp(lp)),A,point['alpha'],point['embedding'],**kw);evs.append({'A':A,'result':v});return v
            root=brentq(lambda lp:ev(lp)['M']-1.4,np.log(1e35),np.log(2e35),xtol=2e-10)
            return ev(root)
        gr=fixed(0.);ns=fixed(point['A'])
        return {'job':job,'settings':kw,'GR':gr,'modified':ns,'delta_R14':ns['R']/gr['R']-1,'evaluations':evs}
    key=(tuple(sorted(settings.items())),interp)
    import scan
    # Immutable legacy files; patch only this worker's in-memory callable.
    with patch.object(scan,'EOS',lambda:EOS(interp=interp)):
        if key not in BASELINES:
            BASELINES[key]=metrics(dict(A=0.,alpha=.1,rstar_cm=1e6,embedding='emb1',settings=settings,pressure_grid_n=18))
        baseline=BASELINES[key]
        if mode=='saturated':
            with patch.object(model,'delta',lambda rho,A,alpha:-3*A):value=metrics(p)
        else:value=metrics(p)
    return {'job':job,'interpolation':interp,'law':'constant_delta_minus_3A' if mode=='saturated' else 'full_density_law','baseline':baseline,'modified':value,'delta_mass':value['domain_peak']['M']/baseline['domain_peak']['M']-1,'delta_R14':value['R14']['R']/baseline['R14']['R']-1}

def run(parent,refinement,out,workers):
    from resource_control import require_controls
    require_controls()
    if workers!=1:raise ValueError("Exactly one worker is permitted")
    out.mkdir(parents=True,exist_ok=False);_,rows,_=load(parent);rows+=refined(refinement,parent)
    jobs=[];seen=set()
    def add(point,mode,reason):
        key=(point['alpha'],point['A'],point['rstar_cm'],point['embedding'],mode)
        if key not in seen:jobs.append({'point':point,'mode':mode,'reason':reason});seen.add(key)
    for ref in REFERENCES:
        for emb in EMBEDDINGS:
            vals=[r for r in rows if r['point']['rstar_cm']==ref and r['point']['embedding']==emb]
            for q in ['delta_mass','delta_R14']:add(min(vals,key=lambda r:r[q])['point'],'tight','refined_extremum_'+q)
    for emb in EMBEDDINGS:
        for q in ['delta_mass','delta_R14']:
            best=min((r for r in rows if r['point']['embedding']==emb),key=lambda r:r[q])
            for mode in ['start10','start1000','pchip','direct']:add(best['point'],mode,'strongest_across_references_'+q)
        for origin,alpha in [('low_alpha',.02),('high_alpha',.5),('solar_peak',None)]:
            group=[r for r in rows if r['point']['embedding']==emb and r['point']['rstar_cm']==1e6 and r['point']['fraction']==1 and (r['point']['alpha']==alpha if alpha else r['point']['origin']=='solar_amplitude_peak')]
            if len(group)!=1:raise ValueError('control selection')
            for mode in ['tight','saturated']:add(group[0]['point'],mode,origin)
    files=[HERE/n for n in ['checks.py','scan.py','refine.py','analyze.py','resource_control.py']]+[ROOT/'paper3_NeutronStars/ns_step14'/n for n in ['model.py','direct.py']]
    manifest={'status':'running','parent_manifest_sha256':sha(parent/'manifest.json'),'refinement_manifest_sha256':sha(refinement/'manifest.json'),'source_sha256':{str(p.relative_to(ROOT)):sha(p) for p in files},'jobs':jobs,'workers':workers}
    def save():(out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    save();values=[];journal=(out/'evaluations.jsonl').open('w')
    try:
        for value in map(task,jobs):
            p=value['job']['point'];original=next(r for r in rows if r['point']==p)
            value['original_delta_mass']=original['delta_mass'];value['original_delta_R14']=original['delta_R14']
            value['difference_R14']=value['delta_R14']-original['delta_R14']
            if 'delta_mass' in value:value['difference_mass']=value['delta_mass']-original['delta_mass']
            if not np.isfinite(value['delta_R14']):raise RuntimeError('nonfinite check')
            values.append(value);journal.write(json.dumps(value)+'\n');journal.flush()
            print(len(values),'of',len(jobs),value['job']['mode'],value['difference_R14'],flush=True)
        journal.close()
        summary={}
        for mode in sorted(set(j['mode'] for j in jobs)):
            group=[v for v in values if v['job']['mode']==mode]
            summary[mode]={'cases':len(group),'max_abs_R14_change':max(abs(v['difference_R14']) for v in group),'max_abs_mass_change':max((abs(v['difference_mass']) for v in group if 'difference_mass' in v),default=None)}
        if summary['tight']['max_abs_R14_change']>1e-6 or summary['tight']['max_abs_mass_change']>1e-6:raise RuntimeError('screening numerical changes exceed admission threshold')
        if summary['direct']['max_abs_R14_change']>1e-6:raise RuntimeError('independent field mismatch exceeds admission threshold')
        (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
        manifest.update(status='complete',rows=len(values),sha256={n:sha(out/n) for n in ['evaluations.jsonl','summary.json']});save()
    except Exception as exc:
        journal.close();manifest.update(status='failed',error=repr(exc),traceback=traceback.format_exc());save();raise
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--scan',type=Path,required=True);p.add_argument('--refinement',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--workers',type=int,default=1,choices=[1]);a=p.parse_args();run(a.scan,a.refinement,a.output,a.workers)
