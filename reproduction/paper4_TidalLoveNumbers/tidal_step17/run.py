"""Sequential retained conditional tides and paired controls; no silent skips."""
from pathlib import Path
import argparse,json,hashlib,platform,traceback
import numpy as np
import scipy
from tides import ROOT,EOS,fixed,pair
from resource_control import require_controls,measurement
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 require_controls();p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
 gates=json.loads((HERE/'gates/gates.json').read_text());assert gates['status']=='complete'
 prior=ROOT/'paper3_NeutronStars/ns_step15/resumed95_results';m=json.loads((prior/'manifest.json').read_text());assert m['status']=='complete'
 for f,h in m['contract']['source_sha256'].items():assert sha(ROOT/f)==h,f
 for f,h in m['sha256'].items():assert sha(prior/f)==h,f
 points=m['contract']['points'];peak=next(p for p in points if p['origin']=='solar_amplitude_peak')
 low=min(points,key=lambda p:p['alpha']);high=max(points,key=lambda p:p['alpha']);ns=min(points,key=lambda p:abs(p['alpha']-.12859203321989024))
 sources=[HERE/n for n in ['tides.py','run.py','gates.py','launch.py']]+[ROOT/'paper3_NeutronStars/ns_step14/model.py',ROOT/'paper3_NeutronStars/ns_step15/resource_control.py',ROOT/'paper3_NeutronStars/ns_v7/eos.nb',ROOT/'paper3_NeutronStars/ns_v7/eos.thermo',ROOT/'paper4_TidalLoveNumbers/tidal_v6/operator_tidal.py',ROOT/'paper4_TidalLoveNumbers/tidal_v5/tidal.py']
 manifest={'status':'running','scope':'Finite sampled conditional GR-form tides; prescribed Schwarzschild/cutoff extraction, not a covariant completion or EOS/global bound','source_sha256':{str(f.relative_to(ROOT)):sha(f) for f in sources},'parent_scan_manifest_sha256':sha(prior/'manifest.json'),'gates_sha256':sha(HERE/'gates/gates.json'),'runtime':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,'workers':1},'primary_matching':'GR_finite_pressure','settings':{'rtol':2e-9,'atol':1e-17,'max_step':.08,'surface_factor':1.001},'resource_start':measurement()}
 def save():(a.output/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
 save();rows=[];cache={};pairs=[];controls=[]
 journal=(a.output/'evaluations.jsonl').open('w')
 def calc(mass,A,alpha,emb,settings,label):
  key=(mass,A,alpha,emb,tuple(sorted(settings.items())))
  if key not in cache:
   kw=dict(settings);interp=kw.pop('interp','linear');v=fixed(EOS(interp=interp),mass,A,alpha,emb,**kw)
   v['interpolation']=interp;v['id']=len(rows);v['label']=label;rows.append(v);cache[key]=v
   journal.write(json.dumps(v)+'\n');journal.flush();print('completed',v['id'],label,mass,emb,'Lambda',v['star']['Lambda'],flush=True)
  return cache[key]
 def make_pair(mass,p,emb,settings,label):
  gr=calc(mass,0.,.12,'emb1',settings,'GR/'+label);mod=calc(mass,p['boundary_A'],p['alpha'],emb,settings,label)
  return {'label':label,'mass':mass,'alpha':p['alpha'],'A':p['boundary_A'],'embedding':emb,'settings':settings,'GR_id':gr['id'],'modified_id':mod['id'],**pair(gr,mod)}
 try:
  for mass in [1.2,1.27,1.4,1.46,1.6]:
   for emb in ['emb1','emb2']:pairs.append(make_pair(mass,peak,emb,{},'solar_peak'))
  for label,pnt in [('low_alpha',low),('high_alpha',high),('NS_extremum_sample',ns)]:
   for emb in ['emb1','emb2']:pairs.append(make_pair(1.4,pnt,emb,{},label))
  for ref in [5e5,2e6]:
   for emb in ['emb1','emb2']:pairs.append(make_pair(1.4,peak,emb,{'rstar':ref},'reference_'+str(ref)))
  pairs.append(make_pair(1.4,{'alpha':.122,'boundary_A':.0087},'emb1',{},'historical_parameters'))
  modes={'tight':{'rtol':2e-11,'atol':1e-19,'max_step':.03},'linear_H':{'form':'H','rtol':2e-11,'atol':1e-19,'max_step':.03},'start10':{'r0':10.},'start1000':{'r0':1000.},'pchip':{'interp':'pchip'},'surface10':{'surface_factor':10.},'surface100':{'surface_factor':100.},'Einstein_nu_diagnostic':{'metric':'einstein'}}
  for label,pnt,emb in [('solar_peak',peak,'emb1'),('solar_peak',peak,'emb2'),('low_alpha',low,'emb1')]:
   orig=next(v for v in pairs if v['label']==label and v['embedding']==emb and v['mass']==1.4)
   for mode,settings in modes.items():
    v=make_pair(1.4,pnt,emb,settings,label+'/'+mode);v['mode']=mode;v['original_label']=label
    for q in ['delta_Lambda','delta_k2','delta_R']:v['change_'+q]=v[q]-orig[q]
    controls.append(v)
  summary={'pairs':pairs,'controls':controls,'control_maxima':{mode:{q:max(abs(v['change_'+q]) for v in controls if v['mode']==mode) for q in ['delta_Lambda','delta_k2','delta_R']} for mode in modes}}
  for mode in ['tight','linear_H']:
   if summary['control_maxima'][mode]['delta_Lambda']>1e-6:raise RuntimeError('Paired tidal convergence/H check exceeds 1e-6')
  # Exact extraction sensitivity is recomputed from retained interior output.
  summary['matching_sensitivity']=[]
  for v in pairs:
   gr=rows[v['GR_id']]['star'];mod=rows[v['modified_id']]['star']
   shifts={name:mod['extractions'][name]['Lambda']/gr['extractions'][name]['Lambda']-1 for name in gr['extractions']}
   summary['matching_sensitivity'].append({'label':v['label'],'mass':v['mass'],'embedding':v['embedding'],'shifts':shifts,'span':max(shifts.values())-min(shifts.values())})
  def binary(m1,m2,l1,l2):return 16/13*((m1+12*m2)*m1**4*l1+(m2+12*m1)*m2**4*l2)/(m1+m2)**5
  summary['binary']={}
  for emb in ['emb1','emb2']:
   x=[next(v for v in pairs if v['label']=='solar_peak' and v['embedding']==emb and v['mass']==mass) for mass in [1.46,1.27]]
   lg=binary(1.46,1.27,*[rows[v['GR_id']]['star']['Lambda'] for v in x]);lm=binary(1.46,1.27,*[rows[v['modified_id']]['star']['Lambda'] for v in x]);summary['binary'][emb]={'masses':[1.46,1.27],'GR':lg,'modified':lm,'fractional_shift':lm/lg-1}
  journal.close();(a.output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
  manifest.update(status='complete',rows=len(rows),pairs=len(pairs),controls=len(controls),resource_end=measurement(),sha256={n:sha(a.output/n) for n in ['evaluations.jsonl','summary.json']});save()
 except BaseException as e:
  journal.close();manifest.update(status='failed',rows=len(rows),error=repr(e),traceback=traceback.format_exc(),resource_end=measurement());save();raise
if __name__=='__main__':main()
