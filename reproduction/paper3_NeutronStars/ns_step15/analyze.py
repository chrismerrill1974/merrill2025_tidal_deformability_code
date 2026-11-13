"""Bind sampled/refined NS rows to solar inputs; recompute extrema and sensitivity.
Failed/incomplete/mismatched runs are refused. Report finite sampling, not bounds.
"""
from pathlib import Path
import argparse,csv,json,sys
import numpy as np
from scan import ROOT,HERE,REFERENCES,EMBEDDINGS,FRACTIONS,inputs,sha
from refine import load

def refined(folder,parent):
    m=json.loads((folder/'manifest.json').read_text())
    if m['status']!='complete' or m['parent_manifest_sha256']!=sha(parent/'manifest.json'):raise ValueError('Incomplete/incompatible refinement')
    for f,h in m['source_sha256'].items():
        if sha(ROOT/f)!=h:raise ValueError('Refinement source mismatch')
    for f,h in m['sha256'].items():
        if sha(folder/f)!=h:raise ValueError('Refinement checksum mismatch')
    return [json.loads(v) for v in (folder/'evaluations.jsonl').read_text().splitlines()]

def run(parent,refinement,out):
    out.mkdir(parents=True,exist_ok=False);points,_,_=inputs();manifest,rows,baseline=load(parent)
    assert len(rows)==1980 and manifest['rows']==1980
    expected={(p['alpha'],p['boundary_A']*f,f,ref,emb) for p in points for f in FRACTIONS for ref in REFERENCES for emb in EMBEDDINGS}
    actual={(r['point']['alpha'],r['point']['A'],r['point']['fraction'],r['point']['rstar_cm'],r['point']['embedding']) for r in rows}
    assert actual==expected and len(actual)==len(rows)
    endpoints={p['alpha']:p['boundary_A'] for p in points}
    original=list(rows)
    if refinement:
        extra=refined(refinement,parent);eps=json.loads((refinement/'solar_endpoints.json').read_text())
        sr=json.loads((refinement/'solar_evaluations.json').read_text())
        sys.path.insert(0,str(ROOT/'paper1_SolarCore/solar_step11'))
        from inference import score
        scores={}
        for r in sr:
            assert r['A']==r['response']['A'] and r['alpha']==r['response']['alpha']
            assert abs(score(r['response'],manifest['contract']['variant'])-r['score'])<1e-12
            key=(r['alpha'],r['A']);assert key not in scores;scores[key]=r['score']
        for ep in eps:
            a=ep['alpha'];A=ep['boundary_A'];v=ep['solar_slice']
            assert abs(scores[(a,A)]-v['minimum_Q']-2.71)<1e-6
            assert abs(scores[(a,v['minimum_A'])]-v['minimum_Q'])<1e-10
            assert all((a,float(x)) in scores for x in np.linspace(0,.04,17))
            assert len(v['intervals'])==1 and v['intervals'][0]['lower']==0 and v['upper_endpoint']==A
            endpoints[a]=A
        rounds=json.loads((refinement/'rounds.json').read_text());expected_extra=set()
        for round in rounds:
            expected_extra |= {(a,endpoints[a]*f,f,ref,emb) for a in round['new_alphas'] for f in round['fractions'] for ref in REFERENCES for emb in EMBEDDINGS}
        found_extra={(r['point']['alpha'],r['point']['A'],r['point']['fraction'],r['point']['rstar_cm'],r['point']['embedding']) for r in extra}
        assert expected_extra==found_extra and len(extra)==len(found_extra)
        rows+=extra
    else:extra=[];rounds=[]
    flat=[]
    for r in rows:
        p=r['point'];assert p['boundary_A']==endpoints[p['alpha']] and p['A']==p['boundary_A']*p['fraction']
        assert np.isfinite([r['delta_mass'],r['delta_R14']]).all()
        assert abs(r['R14']['M']-1.4)<5e-7
        assert r['R14'] in r['evaluations'] and r['domain_peak'] in r['evaluations']
        assert len(r['pressure_grid'])==13
        assert r['domain_peak']['M']>=max(v['M'] for v in r['pressure_grid'])-1e-9
        for star in r['evaluations']:
            assert all(np.isfinite(v) and v>0 for v in star.values())
            assert star['P_surface']==baseline['R14']['P_surface']
        assert r['delta_mass']==r['domain_peak']['M']/baseline['domain_peak']['M']-1
        assert r['delta_R14']==r['R14']['R']/baseline['R14']['R']-1
        if p['A']==0:
            assert r['reused_GR'] and r['delta_mass']==r['delta_R14']==0
            assert r['domain_peak']==baseline['domain_peak'] and r['R14']==baseline['R14']
        flat.append({**p,'domain_peak_mass':r['domain_peak']['M'],'domain_peak_Pc':r['domain_peak']['Pc'],'peak_status':r['peak_status'],'R14_km':r['R14']['R'],'Pc14':r['R14']['Pc'],'delta_mass':r['delta_mass'],'delta_R14':r['delta_R14']})
    extrema={};monotonic=[];coarse={};sensitivity=[]
    coarse_alphas={p['alpha'] for p in points[:-1:2]}|{points[-1]['alpha']}
    for ref in REFERENCES:
        for emb in EMBEDDINGS:
            vals=[v for v in flat if v['rstar_cm']==ref and v['embedding']==emb];key=f'{ref:g}/{emb}'
            extrema[key]={q:{'min':min(vals,key=lambda v:v[q]),'max':max(vals,key=lambda v:v[q])} for q in ['delta_mass','delta_R14']}
            subset=[v for v in vals if v['alpha'] in coarse_alphas];main=[v for v in vals if v['alpha'] in {p['alpha'] for p in points}]
            coarse[key]={q:{'coarse_min':min(v[q] for v in subset),'original_fine_min':min(v[q] for v in main),'refined_min':min(v[q] for v in vals)} for q in ['delta_mass','delta_R14']}
            for p in points:
                group=sorted([r for r in original if r['point']['alpha']==p['alpha'] and r['point']['rstar_cm']==ref and r['point']['embedding']==emb],key=lambda r:r['point']['fraction'])
                assert [r['point']['fraction'] for r in group]==FRACTIONS
                monotonic.append({'alpha':p['alpha'],'rstar_cm':ref,'embedding':emb,**{q:{'largest_increment':float(max(np.diff([r[q] for r in group]))),'strictly_nonincreasing':bool(np.all(np.diff([r[q] for r in group])<=0))} for q in ['delta_mass','delta_R14']}})
    for p in points:
        for emb in EMBEDDINGS:
            for fraction in FRACTIONS[1:]:
                vals=[r for r in flat if r['alpha']==p['alpha'] and r['embedding']==emb and r['fraction']==fraction]
                assert len(vals)==3
                sensitivity.append({'alpha':p['alpha'],'fraction':fraction,'embedding':emb,**{q:max(v[q] for v in vals)-min(v[q] for v in vals) for q in ['delta_mass','delta_R14']}})
    summary={'status':'verified_sampled_results','original_rows':len(original),'refined_rows':len(extra),'alpha_count':len(endpoints),'extrema':extrema,'grid_comparison':coarse,'amplitude_traces':len(monotonic),'amplitude_monotonic':{q:all(v[q]['strictly_nonincreasing'] for v in monotonic) for q in ['delta_mass','delta_R14']},'amplitude_reversals':[v for v in monotonic if any(not v[q]['strictly_nonincreasing'] for q in ['delta_mass','delta_R14'])],'largest_reference_span':{q:max(sensitivity,key=lambda v:v[q]) for q in ['delta_mass','delta_R14']},'peak_status_counts':{v:sum(r['peak_status']==v for r in flat) for v in set(r['peak_status'] for r in flat)},'global_envelope_certified':False,'rounds':rounds}
    for name,value in [('summary.json',summary),('monotonicity.json',monotonic),('reference_sensitivity.json',sensitivity)]: (out/name).write_text(json.dumps(value,indent=2)+'\n')
    with (out/'predictions.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(flat[0]),lineterminator='\n');writer.writeheader();writer.writerows(flat)
    m={'status':'complete','scan_manifest_sha256':sha(parent/'manifest.json'),'refinement_manifest_sha256':sha(refinement/'manifest.json') if refinement else None,'source_sha256':{str(HERE.joinpath(n).relative_to(ROOT)):sha(HERE/n) for n in ['analyze.py','scan.py','refine.py']},'sha256':{n:sha(out/n) for n in ['summary.json','predictions.csv','monotonicity.json','reference_sensitivity.json']}}
    (out/'manifest.json').write_text(json.dumps(m,indent=2)+'\n');print(json.dumps({k:summary[k] for k in ['original_rows','refined_rows','alpha_count','amplitude_monotonic','peak_status_counts']},indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--scan',type=Path,required=True);p.add_argument('--refinement',type=Path);p.add_argument('--output',type=Path,required=True);a=p.parse_args();run(a.scan,a.refinement,a.output)
