"""Admit complete runs and recompute paired numerical sensitivities from bound rows."""
from pathlib import Path
import argparse,hashlib,json
import numpy as np
from model import ROOT,EOS

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def verify(folder,output):
    m=json.loads((folder/'manifest.json').read_text())
    if m['status']!='complete':raise ValueError('incomplete run')
    for f,h in m['contract']['source_sha256'].items():
        if sha(ROOT/f)!=h:raise ValueError('source hash '+f)
    for f,h in m['sha256'].items():
        if sha(folder/f)!=h:raise ValueError('artifact hash '+f)
    s=json.loads((folder/'summary.json').read_text());rows=[json.loads(v) for v in (folder/'evaluations.jsonl').read_text().splitlines()]
    assert len(rows)==m['evaluations']==s['evaluation_count']
    assert [r['id'] for r in rows]==list(range(len(rows)))
    for r in rows:
        assert r['status']=='complete'
        assert all(np.isfinite(v) and v>0 for v in r['result'].values())
        assert r['Pc']==r['result']['Pc']
        e=EOS(r['eos'],r['interpolation'])
        assert e.P_min<r['Pc']<=e.P_max*(1+1e-13)
        assert r['result']['P_surface']==e.P_min*r['settings'].get('surface_factor',1.001)
    def bound(value,tag,kind,interp,A,emb,settings):
        hits=[r for r in rows if r['tag']==tag and r['eos']==kind and r['interpolation']==interp and r['A']==A and r['embedding']==emb and r['settings']==settings and r['alpha']==.12 and r['Pc']==value['Pc']]
        assert len(hits)==1 and hits[0]['result']==value,(tag,value)
    required={'default','tight','tighter','start10','start1000','surface1','surface10','surface100','pchip','legacy_seed','RK45'}
    assert set(s['configurations'])==required
    for name,config in s['configurations'].items():
        rs=config['results'];assert set(rs)=={'GR','emb1','emb2'}
        for label,result in rs.items():
            assert len(result['grid'])==22
            for v in result['grid']+[result['domain_peak'],result['fixed_1p4']]:
                bound(v,name+'/'+label,'compose','pchip' if name=='pchip' else 'linear',0. if label=='GR' else .0087,'emb1' if label=='GR' else label,config['settings'])
            assert abs(result['fixed_1p4']['M']-1.4)<1e-7
            if result['peak_status']=='upper_EOS_boundary':
                assert result['domain_peak']==result['grid'][-1] and result['grid'][-1]['M']>result['grid'][-2]['M']
            if label!='GR':
                assert result['delta_domain_peak_fraction']==result['domain_peak']['M']/rs['GR']['domain_peak']['M']-1
                assert result['delta_R14_fraction']==result['fixed_1p4']['R']/rs['GR']['fixed_1p4']['R']-1
    for kind,result in s['dh_diagnosis'].items():
        for v in result['grid']+[result['domain_peak'],result['fixed_1p4']]:bound(v,kind,kind,'linear',0.,'emb1',s['configurations']['tight']['settings'])
    assert 2.04<s['dh_diagnosis']['dh']['domain_peak']['M']<2.06
    assert 1.78<s['dh_diagnosis']['dh_bug']['domain_peak']['M']<1.81
    for r in s['direct_comparisons']:
        for solver,tag in [('flux','cross/flux'),('direct','cross/direct')]:
            hits=[x for x in rows if x['tag']==tag and x['A']==r['A'] and x['alpha']==r['alpha'] and x['Pc']==r['Pc'] and x['embedding']==r['embedding'] and (solver=='flux' or x['settings']['derivative']==r['derivative'])]
            assert len(hits)==1 and hits[0]['result']==r[solver]
        for q in ['M','R']:assert r['relative_'+q+'_difference']==r['direct'][q]/r['flux'][q]-1
    for p in s['independent_paired_differences']:
        def get(A,al):return next(r for r in s['direct_comparisons'] if r['A']==A and r['alpha']==al and all(r[k]==p[k] for k in ['Pc','embedding','derivative']))
        r=get(p['A'],p['alpha']);g=get(0.,.12)
        for q in ['M','R']:
            for method in ['flux','direct']:assert p[method+'_delta_'+q]==r[method][q]/g[method][q]-1
            assert p['delta_'+q+'_difference']==p['direct_delta_'+q]-p['flux_delta_'+q]
    categories={'numerical_controls':['default','tighter','start10','start1000','RK45'],
                'seed_repair':['legacy_seed'],'interpolation':['pchip'],
                'nearest_table_floor':['surface1'],'raised_surface_floor':['surface10','surface100']}
    uncertainty={}
    for emb in ['emb1','emb2']:
        baseline=s['configurations']['tight']['results'][emb];uncertainty[emb]={}
        for group,names in categories.items():
            uncertainty[emb][group]={key:max(abs(s['configurations'][n]['results'][emb][key]-baseline[key]) for n in names) for key in ['delta_domain_peak_fraction','delta_R14_fraction']}
    independent={q:max(abs(v['delta_'+q+'_difference']) for v in s['independent_paired_differences']) for q in ['M','R']}
    assert max(independent.values())<1e-6
    ref=s['reference_comparison'];assert abs(ref['fractional_domain_mass_difference'])<.002 and abs(ref['fractional_R14_difference'])<.002
    result={'status':'pass_with_EOS_boundary_and_surface_scope_limits','evaluations_verified':len(rows),'configurations_verified':len(required),'independent_paired_cases':len(s['independent_paired_differences']),'observed_max_absolute_changes_in_fractional_shifts':uncertainty,'independent_paired_max_absolute_difference':independent,'meaning':'observed changes in this bounded suite, not certified errors or a solar-parameter envelope','reference_comparison':ref,'source_sha256':{str(Path(__file__).resolve().relative_to(ROOT)):sha(Path(__file__))}}
    output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--results',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();verify(a.results,a.output)
