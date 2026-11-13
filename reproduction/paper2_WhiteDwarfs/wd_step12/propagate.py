"""Step 12: propagate accepted solar-a slices through the step 07 audited flux solver.
Samples are not a certified continuous-domain envelope. New output directory only.
"""
from pathlib import Path
import argparse,csv,json,sys
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT));sys.path.insert(0,str(ROOT/'paper2_WhiteDwarfs/audit_2026_11_11'))
from audit import flux,shifts
from numerical_guardrails import fingerprint,sha
SOLAR=ROOT/'paper1_SolarCore/solar_step11/results/solar_a'
PEAK=ROOT/'paper1_SolarCore/solar_step11/results/alpha_peak_checks.json'
HERE=Path(__file__).resolve().parent

def load_inputs():
    m=json.loads((SOLAR/'manifest.json').read_text())
    if m['status']!='complete':raise ValueError('Solar input incomplete')
    for f,h in m['sha256'].items():
        if sha(SOLAR/f)!=h:raise ValueError('Solar checksum mismatch')
    for f,h in m['contract']['source_sha256'].items():
        if sha(SOLAR.parent.parent/f)!=h:raise ValueError('Solar source mismatch')
    pinned=json.loads((ROOT/'2026_june_14_make_it_consistent/STEP11_ARTIFACTS.json').read_text())['sha256']
    if sha(PEAK)!=pinned[str(PEAK.relative_to(ROOT))]:raise ValueError('Peak checksum mismatch')
    region=json.loads((SOLAR/'region.json').read_text())
    points=[]
    for s in region['slices']:
        v=s['variants']['source_precision_diagonal_proxy_v1']
        if len(v['intervals'])!=1 or v['intervals'][0]['lower']!=0 or v['upper_endpoint'] is None:raise ValueError('Expected closed connected solar slices')
        points.append((s['alpha'],v['upper_endpoint'],'sampled_slice'))
    p=json.loads(PEAK.read_text())['checks']['source_precision_diagonal_proxy_v1']
    points.append((p['alpha'],p['endpoint'],'local_solar_peak'))
    return points

def run(output):
    out=Path(output);out.mkdir(parents=True,exist_ok=False)
    points=load_inputs();densities=np.geomspace(1e5,3e6,12);fractions=[0,.25,.5,.75,1.]
    meta={'status':'running','contract':fingerprint(ROOT,['paper2_WhiteDwarfs/wd_step12/propagate.py','paper2_WhiteDwarfs/audit_2026_11_11/audit.py','paper2_WhiteDwarfs/wd_v5/model.py'],{'solar_sha256':sha(SOLAR/'region.json'),'peak_sha256':sha(PEAK),'rho_c':densities.tolist(),'fractions':fractions,'points':points,'solver':'step07 flux RK45, no derivative cutoff','rtol':2e-11,'atol':2e-13,'step':.01,'eps':1e-6,'comparison_mass':'modified mass M_A','fractional_radius':'expm1(dlnR + dlnM/3)'})}
    def save(): (out/'manifest.json').write_text(json.dumps(meta,indent=2)+'\n')
    save();rows=[];point=None
    try:
        baseline={float(r):flux(0,.1,float(r)) for r in densities}
        b=baseline[float(densities[0])]
        if abs(b['xi']/3.65375374-1)>1e-8 or abs(b['I']/2.71405512-1)>1e-8:raise RuntimeError('Baseline gate failed')
        for alpha,limit,label in points:
            for fraction in fractions:
                A=limit*fraction
                for rho in densities:
                    point=dict(alpha=alpha,A=A,rho_c=float(rho))
                    s=baseline[float(rho)];m=s if A==0 else flux(A,alpha,float(rho));dr,dm,fixed=map(float,shifts(s,m))
                    row={**point,'boundary_A':limit,'fraction':fraction,'origin':label,'M_std':s['M']/1.98847e33,'M_mod':m['M']/1.98847e33,'R_std_km':s['R']/1e5,'R_mod_km':m['R']/1e5,'dlnR':dr,'dlnM':dm,'dlnR_fixedM':fixed,'fractional_fixedM':float(np.expm1(fixed))}
                    if not all(np.isfinite(v) for k,v in row.items() if k!='origin'):raise RuntimeError('Nonfinite row')
                    rows.append(row)
                    with (out/'evaluations.jsonl').open('a') as f:f.write(json.dumps(row)+'\n')
            print(alpha,'complete',len(rows),flush=True)
        selected=[]
        for alpha,limit,label in [points[0],points[32],points[-2],points[-1]]:
            for rho in [1e5,1e6,3e6]:
                s,m=flux(0,alpha,rho),flux(limit,alpha,rho)
                base=float(np.expm1(shifts(s,m)[2]))
                item=dict(alpha=alpha,A=limit,rho_c=rho,baseline=base,comparisons={})
                for name,kw in [('tight',dict(rtol=2e-12,atol=2e-14,step=.003)),('start',dict(eps=1e-4))]:
                    st,mt=flux(0,alpha,rho,**kw),flux(limit,alpha,rho,**kw)
                    value=float(np.expm1(shifts(st,mt)[2]));item['comparisons'][name]={'value':value,'difference':value-base}
                selected.append(item)
        worst=min(rows,key=lambda r:r['fractional_fixedM']);nonzero=[r for r in rows if r['A']>0]
        monotonic=True
        for alpha,_,_ in points:
            for rho in densities:
                vals=[r['fractional_fixedM'] for r in rows if r['alpha']==alpha and r['rho_c']==rho]
                monotonic &= bool(np.all(np.diff(vals)<=1e-12))
        summary={'rows':len(rows),'alpha_samples':len(points),'rho_samples':len(densities),'amplitude_fractions':fractions,'worst_sample':worst,'nonzero_fractional_range':[min(r['fractional_fixedM'] for r in nonzero),max(r['fractional_fixedM'] for r in nonzero)],'sampled_amplitude_monotonic':monotonic,'M_mod_range':[min(r['M_mod'] for r in rows),max(r['M_mod'] for r in rows)],'boundary_M_mod_range':[min(r['M_mod'] for r in rows if r['fraction']==1),max(r['M_mod'] for r in rows if r['fraction']==1)],'selected_convergence':selected,'max_tight_difference':max(abs(v['comparisons']['tight']['difference']) for v in selected),'max_start_difference':max(abs(v['comparisons']['start']['difference']) for v in selected),'global_envelope_certified':False}
        (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
        with (out/'predictions.csv').open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
        if len(rows)!=66*12*5:raise RuntimeError('Incomplete grid')
        meta.update(status='complete',rows=len(rows),sha256={f:sha(out/f) for f in ['evaluations.jsonl','predictions.csv','summary.json']});save()
    except Exception as e:
        meta.update(status='failed',point=point,error=str(e));save();raise
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);a=p.parse_args();run(a.output)
