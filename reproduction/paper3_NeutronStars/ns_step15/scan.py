"""Step 15: finite-grid propagation of the complete selected solar slice set.
Both conditional NS prescriptions; fixed physical reference lengths. All writes
are parent-owned; workers only solve stars. No extrapolation or failed-row reuse.
"""
from pathlib import Path
import argparse,hashlib,json,os,platform,sys,traceback
import numpy as np
import scipy
from scipy.optimize import brentq,minimize_scalar
ROOT=Path(__file__).resolve().parents[2];HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'paper3_NeutronStars/ns_step14'))
from model import EOS,solve
VARIANT='source_precision_diagonal_proxy_v1'
FRACTIONS=[0.,.25,.5,.75,1.];REFERENCES=[5e5,1e6,2e6];EMBEDDINGS=['emb1','emb2']
SETTINGS={'rtol':2e-9,'atol':1e-17,'max_step':.1}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def inputs():
    old=json.loads((ROOT/'2026_june_14_make_it_consistent/STEP11_ARTIFACTS.json').read_text())
    for p,h in old['sha256'].items():
        if sha(ROOT/p)!=h:raise ValueError('Solar source/artifact mismatch '+p)
    m=json.loads((ROOT/'paper3_NeutronStars/ns_step14/results/manifest.json').read_text())
    if m['status']!='complete':raise ValueError('Incomplete NS audit')
    for p,h in m['contract']['source_sha256'].items():
        if sha(ROOT/p)!=h:raise ValueError('NS audit source mismatch '+p)
    region=json.loads((ROOT/'paper1_SolarCore/solar_step11/results/solar_a/region.json').read_text());points=[]
    for s in region['slices']:
        v=s['variants'][VARIANT]
        if len(v['intervals'])!=1 or v['intervals'][0]['lower']!=0 or v['upper_endpoint'] is None:raise ValueError('Not a closed connected solar slice')
        points.append({'alpha':s['alpha'],'boundary_A':v['upper_endpoint'],'origin':'solar_slice'})
    peak=json.loads((ROOT/'paper1_SolarCore/solar_step11/results/alpha_peak_checks.json').read_text())['checks'][VARIANT]
    if peak['status']!='complete':raise ValueError('Solar local refinement incomplete')
    points.append({'alpha':peak['alpha'],'boundary_A':peak['endpoint'],'origin':'solar_amplitude_peak'})
    if len(points)!=66:raise ValueError('Unexpected solar grid')
    return points,old,m

def metrics(point):
    """Independent calculations per parameter tuple; every stellar evaluation retained."""
    from resource_control import require_controls
    require_controls()
    eos=EOS();A=point['A'];alpha=point['alpha'];emb=point['embedding'];ref=point['rstar_cm']
    settings=dict(point.get('settings',SETTINGS),rstar=ref);cache={};evs=[]
    def ev(lp):
        lp=float(lp)
        if lp not in cache:
            result=solve(eos,float(np.exp(lp)),A,alpha,emb,**settings);cache[lp]=result;evs.append(result)
        return cache[lp]
    n=point.get('pressure_grid_n',9)
    grid=np.unique(np.r_[np.linspace(np.log(5e34),np.log(eos.P_max),n),np.log(eos.P_max*np.array([.98,.99,.995,.999]))])
    curve=[ev(t) for t in grid];candidates=[(curve[-1],'upper_EOS_boundary')]
    for i in range(1,len(grid)-1):
        if curve[i]['M']>=max(curve[i-1]['M'],curve[i+1]['M']):
            fit=minimize_scalar(lambda lp:-ev(lp)['M'],bounds=(grid[i-1],grid[i+1]),method='bounded',options={'xatol':2e-6})
            if not fit.success:raise RuntimeError('pressure peak refinement failed')
            candidates += [(ev(fit.x),'interior_refined'),(curve[i],'interior_sample')]
    peak,status=max(candidates,key=lambda pair:pair[0]['M'])
    if max(v['M'] for v in curve)>peak['M']+1e-9:raise RuntimeError('peak candidate coverage failed')
    if curve[0]['M']>=peak['M']:raise RuntimeError('lower pressure edge maximum')
    lo,hi=np.log(1e35),np.log(2e35)
    if not ev(lo)['M']<1.4<ev(hi)['M']:
        lo,hi=np.log(2e34),np.log(peak['Pc'])
    if not ev(lo)['M']<1.4<ev(hi)['M']:raise RuntimeError('mass 1.4 not bracketed')
    root=brentq(lambda p:ev(p)['M']-1.4,lo,hi,xtol=2e-9)
    r14=ev(root)
    if abs(r14['M']-1.4)>5e-7:raise RuntimeError('fixed mass root residual')
    return {'point':point,'domain_peak':peak,'peak_status':status,'R14':r14,'pressure_grid':curve,'evaluations':evs}

def summarize(rows):
    worst={}
    for ref in REFERENCES:
        for emb in EMBEDDINGS:
            vals=[r for r in rows if r['point']['rstar_cm']==ref and r['point']['embedding']==emb]
            worst[f'{ref:g}/{emb}']={q:min(vals,key=lambda r:r[q]) for q in ['delta_mass','delta_R14']}
    return worst

def run(out,workers):
    from resource_control import require_controls
    require_controls()
    if workers!=1:raise ValueError("Exactly one worker is permitted")
    out.mkdir(parents=True,exist_ok=False);points,solar,ns=inputs()
    files=[HERE/'scan.py',ROOT/'paper3_NeutronStars/ns_step14/model.py',ROOT/'paper3_NeutronStars/ns_v7/eos.nb',ROOT/'paper3_NeutronStars/ns_v7/eos.thermo',ROOT/'paper1_SolarCore/solar_step11/results/solar_a/region.json',ROOT/'paper1_SolarCore/solar_step11/results/alpha_peak_checks.json']
    manifest={'status':'running','contract':{'source_sha256':{str(p.relative_to(ROOT)):sha(p) for p in files},'solar_artifacts':solar['sha256'],'ns_audit_sources':ns['contract']['source_sha256'],'variant':VARIANT,'points':points,'fractions':FRACTIONS,'reference_lengths_cm':REFERENCES,'embeddings':EMBEDDINGS,'solver_settings':SETTINGS,'pressure_grid_n':9,'scope':'sampled conditional results, no certified continuum or true maximum-mass bound'},'runtime':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,'workers':workers}}
    def save(): (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    save();rows=[];jobs=[];journal=(out/'evaluations.jsonl').open('w')
    try:
        baseline=metrics({'A':0.,'alpha':.1,'embedding':'emb1','rstar_cm':1e6,'origin':'GR_reused','settings':dict(rtol=2e-11,atol=1e-19,max_step=.04)})
        # GR has no alpha/reference dependence. Confirm second modified path, no fallback.
        other=metrics({'A':0.,'alpha':.5,'embedding':'emb2','rstar_cm':5e5,'origin':'GR_crosscheck','settings':dict(rtol=2e-11,atol=1e-19,max_step=.04)})
        (out/'baseline.json').write_text(json.dumps({'primary':baseline,'crosscheck':other},indent=2)+'\n')
        for k,field in [('domain_peak','M'),('R14','R')]:
            if abs(other[k][field]/baseline[k][field]-1)>1e-8:raise RuntimeError('GR path mismatch')
        (out/'baseline.json').write_text(json.dumps({'primary':baseline,'crosscheck':other},indent=2)+'\n')
        for p in points:
            for ref in REFERENCES:
                for emb in EMBEDDINGS:
                    for fraction in FRACTIONS:
                        jobs.append({**p,'A':p['boundary_A']*fraction,'fraction':fraction,'rstar_cm':ref,'embedding':emb})
        def record(value,point):
            r=dict(value);r['point']=point
            r['delta_mass']=r['domain_peak']['M']/baseline['domain_peak']['M']-1
            r['delta_R14']=r['R14']['R']/baseline['R14']['R']-1
            r['reused_GR']=point['A']==0
            journal.write(json.dumps(r)+'\n');journal.flush();rows.append(r)
        for p in jobs:
            if p['A']==0:record(baseline,p)
        tasks=[p for p in jobs if p['A']>0]
        for p in tasks:
            value=metrics(p)
            record(value,p)
            if len(rows)%20==0:print(len(rows),'of',len(jobs),'complete',flush=True)
        journal.close()
        if len(rows)!=1980:raise RuntimeError('incomplete grid')
        # Store extrema by row index, not another copy of thousands of evaluations.
        extrema={key:{q:rows.index(row) for q,row in val.items()} for key,val in summarize(rows).items()}
        (out/'summary.json').write_text(json.dumps({'rows':len(rows),'nonzero_rows':len(tasks),'stellar_evaluations':sum(len(r['evaluations']) for r in rows if not r['reused_GR'])+len(baseline['evaluations'])+len(other['evaluations']),'extrema_indices':extrema,'global_envelope_certified':False},indent=2)+'\n')
        manifest.update(status='complete',rows=len(rows),sha256={name:sha(out/name) for name in ['evaluations.jsonl','baseline.json','summary.json']});save()
    except Exception as exc:
        journal.close();manifest.update(status='failed',error=repr(exc),traceback=traceback.format_exc(),rows=len(rows));save();raise
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--workers',type=int,default=1,choices=[1]);a=p.parse_args()
    if a.workers!=1:raise ValueError('Exactly one worker is permitted')
    run(a.output,a.workers)
