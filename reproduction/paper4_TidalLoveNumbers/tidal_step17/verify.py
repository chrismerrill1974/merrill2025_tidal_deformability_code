"""Recompute row bindings, paired shifts and binary output from saved stars."""
from pathlib import Path
import json,hashlib
import numpy as np
from tides import ROOT,pair,love_k2
from resource_control import require_controls
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 require_controls();d=HERE/'results';m=json.loads((d/'manifest.json').read_text());assert m['status']=='complete'
 for p,h in m['source_sha256'].items():assert sha(ROOT/p)==h,p
 for p,h in m['sha256'].items():assert sha(d/p)==h,p
 assert sha(HERE/'gates/gates.json')==m['gates_sha256']
 parent=ROOT/'paper3_NeutronStars/ns_step15/resumed95_results/manifest.json'
 assert sha(parent)==m['parent_scan_manifest_sha256']
 points=json.loads(parent.read_text())['contract']['points']
 rows=[json.loads(line) for line in (d/'evaluations.jsonl').read_text().splitlines()];s=json.loads((d/'summary.json').read_text())
 assert len(rows)==m['rows'] and [v['id'] for v in rows]==list(range(len(rows)))
 assert len(s['pairs'])==21 and len(s['controls'])==24
 peak=next(p for p in points if p['origin']=='solar_amplitude_peak')
 permitted={'solar_peak':peak,'low_alpha':min(points,key=lambda p:p['alpha']),'high_alpha':max(points,key=lambda p:p['alpha']),'NS_extremum_sample':min(points,key=lambda p:abs(p['alpha']-.12859203321989024)),'historical_parameters':{'alpha':.122,'boundary_A':.0087},'reference_500000.0':peak,'reference_2000000.0':peak}
 for v in s['pairs']:
  p=permitted[v['label']];assert (v['alpha'],v['A'])==(p['alpha'],p['boundary_A'])
 assert len({(v['label'],v['mass'],v['embedding']) for v in s['pairs']})==21
 for v in rows:
  a=v['star'];bg=v['background'];assert bg in v['root_evaluations'];assert abs(a['M']-v['target_mass'])<5e-7
  assert all(np.isfinite(a[k]) for k in ['M','R','Pc','P_surface','C','k2','Lambda','y_in'])
  assert a['Pc']==bg['Pc'] and a['P_surface']==bg['P_surface']
  assert abs(a['M']/bg['M']-1)<5e-7 and abs(a['R']/bg['R']-1)<5e-7
  assert abs(a['segments'][-1]['t']-np.log(a['Pc']/a['P_surface']))<1e-12
  assert abs(a['segments'][-1]['state'][1]-a['M'])<1e-14
  for name,ex in a['extractions'].items():
   assert ex['y_out']==a['y_in']-ex['jump'];assert abs(ex['k2']-love_k2(a['C'],ex['y_out']))<1e-14
   assert abs(ex['Lambda']/(2/3*ex['k2']/a['C']**5)-1)<1e-14
  for k in ['k2','Lambda','y_out','jump']:assert a[k]==a['extractions']['GR_finite_pressure'][k]
 for v in s['pairs']+s['controls']:
  gr=rows[v['GR_id']];mod=rows[v['modified_id']]
  assert gr['A']==0 and mod['A']==v['A'] and mod['alpha']==v['alpha'] and mod['embedding']==v['embedding']
  assert gr['target_mass']==mod['target_mass']==v['mass']
  assert gr['settings']==mod['settings']=={k:x for k,x in v['settings'].items() if k!='interp'}
  assert gr['interpolation']==mod['interpolation']==v['settings'].get('interp','linear')
  for k,x in pair(gr,mod).items():assert abs(v[k]-x)<1e-14
  if 'mode' in v:
   original=next(q for q in s['pairs'] if q['label']==v['original_label'] and q['embedding']==v['embedding'] and q['mass']==v['mass'])
   for q in ['delta_Lambda','delta_k2','delta_R']:assert abs(v['change_'+q]-(v[q]-original[q]))<1e-14
 for mode,maxima in s['control_maxima'].items():
  group=[v for v in s['controls'] if v['mode']==mode];assert len(group)==3
  for q,val in maxima.items():assert val==max(abs(v['change_'+q]) for v in group)
  if mode in ['tight','linear_H']:assert maxima['delta_Lambda']<1e-6
 assert len(s['matching_sensitivity'])==21
 for v,sens in zip(s['pairs'],s['matching_sensitivity']):
  assert (v['label'],v['mass'],v['embedding'])==(sens['label'],sens['mass'],sens['embedding'])
  gr=rows[v['GR_id']]['star'];mod=rows[v['modified_id']]['star']
  expected={n:mod['extractions'][n]['Lambda']/gr['extractions'][n]['Lambda']-1 for n in gr['extractions']}
  assert expected==sens['shifts'] and sens['span']==max(expected.values())-min(expected.values())
 for emb,b in s['binary'].items():
  pp=[next(v for v in s['pairs'] if v['label']=='solar_peak' and v['embedding']==emb and v['mass']==mass) for mass in b['masses']]
  x,y=b['masses']
  for kind,key in [('GR','GR_id'),('modified','modified_id')]:
   l1,l2=[rows[v[key]]['star']['Lambda'] for v in pp]
   expected=16/13*((x+12*y)*x**4*l1+(y+12*x)*y**4*l2)/(x+y)**5
   assert abs(b[kind]-expected)<1e-12
  assert b['fractional_shift']==b['modified']/b['GR']-1
 result={'status':'pass','stars_verified':len(rows),'paired_results':len(s['pairs']),'controls':len(s['controls']),'binary_examples':len(s['binary']),'source_hashes':len(m['source_sha256'])}
 (HERE/'verification.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
if __name__=='__main__':main()
