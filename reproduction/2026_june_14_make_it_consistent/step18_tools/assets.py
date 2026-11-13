"""Generate III/IV tables and shared macros solely from verified stored results."""
from pathlib import Path
import json,hashlib,shutil
from resource_control import require_controls
ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'2026_june_14_make_it_consistent/step18_evidence/assets'
P3=ROOT/'paper3_NeutronStars/V8_Paper';P4=ROOT/'paper4_TidalLoveNumbers/V6_Paper'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 require_controls();OUT.mkdir(parents=True,exist_ok=True)
 directories=['paper3_NeutronStars/ns_step15/resumed95_results','paper3_NeutronStars/ns_step15/analysis','paper3_NeutronStars/ns_step15/checks','paper4_TidalLoveNumbers/tidal_step17/results']
 pins={}
 for name in directories:
  d=ROOT/name;m=json.loads((d/'manifest.json').read_text());assert m['status']=='complete'
  for f,h in m.get('contract',{}).get('source_sha256',m.get('source_sha256',{})).items():assert sha(ROOT/f)==h,f
  for f,h in m['sha256'].items():assert sha(d/f)==h,f
  pins[name+'/manifest.json']=sha(d/'manifest.json')
 ns=[json.loads(x) for x in (ROOT/directories[0]/'evaluations.jsonl').read_text().splitlines()]
 ss=json.loads((ROOT/directories[1]/'summary.json').read_text())
 ts=json.loads((ROOT/directories[3]/'summary.json').read_text());tr=[json.loads(x) for x in (ROOT/directories[3]/'evaluations.jsonl').read_text().splitlines()]
 pairs={emb:next(p for p in ts['pairs'] if p['label']=='solar_peak' and p['mass']==1.4 and p['embedding']==emb) for emb in ['emb1','emb2']}
 peak=pairs['emb1'];gr=tr[peak['GR_id']]['star'];baseline=json.loads((ROOT/directories[0]/'baseline.json').read_text())['primary']
 common={emb:next(v for v in ns if v['point']['origin']=='solar_amplitude_peak' and v['point']['fraction']==1 and v['point']['rstar_cm']==1e6 and v['point']['embedding']==emb) for emb in pairs}
 for emb,v in common.items():
  tide=tr[pairs[emb]['modified_id']]['star'];assert abs(tide['R']/v['R14']['R']-1)<5e-7
  assert v['point']['A']==pairs[emb]['A'] and v['point']['alpha']==pairs[emb]['alpha']
 macros={'SolarA':f"{peak['A']:.10f}",'SolarAlpha':f"{peak['alpha']:.10f}",'DimensionFloor':f"{3-3*peak['A']:.9f}",'GRPeakMass':f"{baseline['domain_peak']['M']:.7f}",'GRRadius':f"{gr['R']:.6f}",'GRLambda':f"{gr['Lambda']:.4f}",'GRLove':f"{gr['k2']:.8f}"}
 for emb,suffix in [('emb1','One'),('emb2','Two')]:
  p=pairs[emb];r=tr[p['modified_id']]['star'];b=ts['binary'][emb]
  for key,val in [('Lambda',r['Lambda']),('TidalShift',100*p['delta_Lambda']),('RadiusShift',100*p['delta_R']),('LoveShift',100*p['delta_k2']),('RadiusOnly',100*p['radius_only_exact']),('CrossTerm',100*p['multiplicative_cross_term']),('Mismatch',100*r['metric_derivative_relative_mismatch']),('Binary',b['modified'])]:macros[key+suffix]=f'{val:.5f}'
 macros['BinaryGR']=f"{ts['binary']['emb1']['GR']:.5f}"
 text=''.join('\\newcommand{\\'+k+'}{'+v+'}\n' for k,v in macros.items());(OUT/'shared_results.tex').write_text(text)
 def table(name,header,rows,cols):
  tex='\\begin{tabular}{'+cols+'}\n\\hline\\hline\n'+header+' \\\\\n\\hline\n'+'\n'.join(' & '.join(r)+' \\\\' for r in rows)+'\n\\hline\\hline\n\\end{tabular}\n'
  (OUT/name).write_text(tex)
 rows=[]
 for ref in [5e5,1e6,2e6]:
  for emb,n in [('emb1','1'),('emb2','2')]:
   v=ss['extrema'][f'{ref:g}/{emb}'];rows.append([f'{ref/1e5:g}',n,f"{100*v['delta_mass']['min']['delta_mass']:.6f}",f"{100*v['delta_R14']['min']['delta_R14']:.6f}"])
 table('ns_extrema.tex',r'$r_\star$ (km) & Emb. & $\Delta M_{\rm dom}/M$ (\%) & $\Delta R_{1.4}/R$ (\%)',rows,'cccc')
 rows=[]
 for emb,n in [('emb1','1'),('emb2','2')]:
  v=common[emb];p=pairs[emb];rows.append([n,f"{v['domain_peak']['M']:.7f}",f"{v['R14']['R']:.6f}",f"{tr[p['modified_id']]['star']['Lambda']:.4f}"])
 table('common_point.tex',r'Emb. & $M_{\rm dom}/M_\odot$ & $R_{1.4}$ (km) & $\Lambda_{1.4}$',rows,'cccc')
 rows=[]
 for mass in [1.2,1.27,1.4,1.46,1.6]:
  e1=next(p for p in ts['pairs'] if p['label']=='solar_peak' and p['mass']==mass and p['embedding']=='emb1');e2=next(p for p in ts['pairs'] if p['label']=='solar_peak' and p['mass']==mass and p['embedding']=='emb2')
  rows.append([f'{mass:.2f}',f"{tr[e1['GR_id']]['star']['Lambda']:.4f}",f"{tr[e1['modified_id']]['star']['Lambda']:.4f}",f"{100*e1['delta_Lambda']:.5f}",f"{tr[e2['modified_id']]['star']['Lambda']:.4f}",f"{100*e2['delta_Lambda']:.5f}"])
 table('tidal_samples.tex',r'$M/M_\odot$ & $\Lambda_{\rm GR}$ & $\Lambda_1$ & $\Delta\Lambda_1$ (\%) & $\Lambda_2$ & $\Delta\Lambda_2$ (\%)',rows,'cccccc')
 # Reuse verified scientific figures; no regenerated physics or interpolated envelope.
 figures={'ns_parameter_dependence.pdf':ROOT/'paper3_NeutronStars/ns_step15/figures/parameter_dependence.pdf','tidal_response.pdf':ROOT/'paper4_TidalLoveNumbers/tidal_step17/figures/tidal_response.pdf'}
 for name,p in figures.items():shutil.copyfile(p,OUT/name);pins[str(p.relative_to(ROOT))]=sha(p)
 for dest,files in [(P3,['shared_results.tex','common_point.tex','ns_extrema.tex','ns_parameter_dependence.pdf']),(P4,['shared_results.tex','common_point.tex','tidal_samples.tex','tidal_response.pdf'])]:
  for name in files:shutil.copyfile(OUT/name,dest/name)
 result={'status':'pass','input_manifest_sha256':pins,'macros':macros,'common_point':common,'tidal_point':pairs,'asset_sha256':{p.name:sha(p) for p in OUT.iterdir() if p.is_file()}}
 (OUT.parent/'asset_bindings.json').write_text(json.dumps(result,indent=2)+'\n');print('Shared parameters, backgrounds, tables and assets verified.')
if __name__=='__main__':main()
