"""Public entry point: verify the vendored capsule, regenerate assets, build papers.
No stellar parameter scan is launched. Results go to a new output directory.
"""
from pathlib import Path
import argparse,hashlib,json,os,shutil,signal,sys
BASE=Path(__file__).resolve().parent
p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--capsule',type=Path,default=BASE/'reproduction');p.add_argument('--work',action='append');a=p.parse_args()
ROOT=a.capsule.resolve();sys.path.insert(0,str(ROOT/'paper3_NeutronStars/ns_step15'))
import resource_control as rc
rc.start()
try:
 sys.path.insert(0,str(ROOT/'2026_june_14_make_it_consistent/step18_tools'))
 from build import run_child
 sys.path.insert(0,str(ROOT/'2026_june_14_make_it_consistent/step21_tools'))
 from config import CASES,pdf
 manifest=json.loads((ROOT/'CAPSULE.json').read_text())
 for n,h in manifest['sha256'].items():
  if hashlib.sha256((ROOT/n).read_bytes()).hexdigest()!=h:raise ValueError('Capsule mismatch: '+n)
 selected=a.work or json.loads((BASE/'release.json').read_text())['works']
 assert set(selected)<=set(c[0] for c in CASES)
 out=a.output.resolve();out.mkdir(parents=True,exist_ok=False);records=[];generated={}
 envfile=ROOT/'2026_june_14_make_it_consistent/step21_tools/execute.py'
 def child(args,label):
  rec=run_child([sys.executable,str(envfile),*map(str,args)],ROOT,out/(label+'.log'));records.append(rec)
  (out/'progress.json').write_text(json.dumps(records,indent=2)+'\n')
 def equal(new,old):
  assert new.read_bytes()==old.read_bytes(),f'Regeneration mismatch: {old}'
 if 'paper1' in selected:
  dest=out/'solar';child(['script','paper1_SolarCore/solar_step12/figures.py','--output',dest],'solar_assets')
  for f in ['score_map.png','slice_endpoints.png','slice_table.tex','summary.json']:equal(dest/f,ROOT/'paper1_SolarCore/solar_step12/figures'/f)
  generated['paper1']=dest
 if 'paper2' in selected:
  dest=out/'wd';child(['script','paper2_WhiteDwarfs/wd_step12/figures.py','--output',dest],'wd_assets')
  for f in ['radius_shifts.png','boundary_table.tex']:equal(dest/f,ROOT/'paper2_WhiteDwarfs/wd_step12/figures'/f)
  generated['paper2']=dest
 if any(n in selected for n in ['paper3','paper4','paper5']):
  dest=out/'shared';child(['shared',dest],'shared_tables')
  for n in ['shared_results.tex','common_point.tex','ns_extrema.tex','tidal_samples.tex']:equal(dest/n,ROOT/'2026_june_14_make_it_consistent/step18_evidence/assets'/n)
  generated.update({n:dest for n in ['paper3','paper4','paper5']})
 if 'paper3' in selected:
  dest=out/'ns';child(['script','paper3_NeutronStars/ns_step15/plot.py','--analysis',ROOT/'paper3_NeutronStars/ns_step15/analysis','--output',dest],'ns_figure')
  equal(dest/'parameter_dependence.png',ROOT/'paper3_NeutronStars/ns_step15/figures/parameter_dependence.png')
  shutil.copyfile(dest/'parameter_dependence.pdf',out/'shared/ns_parameter_dependence.pdf')
 if 'paper4' in selected:
  dest=out/'tides';child(['tidal',dest],'tidal_figure')
  equal(dest/'figures/tidal_response.png',ROOT/'paper4_TidalLoveNumbers/tidal_step17/figures/tidal_response.png')
  shutil.copyfile(dest/'figures/tidal_response.pdf',out/'shared/tidal_response.pdf')
 if any(n in selected for n in ['paper6','terrestrial','companion']):
  dest=out/'algebra';child(['algebra',dest],'scalar_lab_audit')
  for n in ['contrast_table.tex','material_table.tex']:equal(dest/'paper6/terrestrial_G_note'/n,ROOT/'paper6/terrestrial_G_note'/n)
  regenerated=json.loads((dest/'2026_june_14_make_it_consistent/step20_evidence/audit.json').read_text());saved=json.loads((ROOT/'2026_june_14_make_it_consistent/step20_evidence/audit.json').read_text())
  for k,v in saved.items():
   if k!='resource':assert regenerated[k]==v,k
  generated['terrestrial']=dest/'paper6/terrestrial_G_note'
 tex=shutil.which('pdflatex') or '/Library/TeX/texbin/pdflatex';extract=shutil.which('pdftotext') or '/opt/homebrew/bin/pdftotext'
 for name,source,assets,repo in CASES:
  if name not in selected:continue
  w=out/name;w.mkdir();shutil.copyfile(ROOT/source,w/'main.tex')
  for folder,files in assets.items():
   for f in files:
    path=generated.get(name,ROOT/folder)/f
    if not path.exists():path=ROOT/folder/f
    shutil.copyfile(path,w/f)
  for n in [1,2]:records.append(run_child([tex,'-interaction=nonstopmode','-halt-on-error','main.tex'],w,w/f'pass{n}.log'))
  log=(w/'main.log').read_text();assert not any(t in log for t in ['Overfull \\hbox','Overfull \\vbox','undefined references','undefined on input line'])
  records.append(run_child([extract,str(w/'main.pdf'),str(w/'text.txt')],w,w/'extract.log'))
  records.append(run_child([extract,str(ROOT/pdf(name)),str(w/'reference.txt')],w,w/'reference.log'))
  equal(w/'text.txt',w/'reference.txt')
 report={'status':'pass','work':selected,'files_verified':len(manifest['sha256']),'regenerated_assets_equal':True,'built_text_equal':True,'scope':'Saved-result asset regeneration, algebra checks, and PDF build; not a new full scientific scan','tools':records,'supervisor_resource':rc.measurement()}
 (out/'REPRODUCTION.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items() if k not in ['tools','supervisor_resource']},indent=2))
finally:signal.setitimer(signal.ITIMER_PROF,0.)
