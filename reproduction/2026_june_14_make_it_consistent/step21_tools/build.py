"""Build all eight release-candidate papers with sequential throttled tools."""
from pathlib import Path
import argparse,json,shutil,sys
from resource_control import require_controls
from config import ROOT,CASES,pdf
require_controls()
sys.path.insert(0,str(ROOT/'2026_june_14_make_it_consistent/step18_tools'))
from build import run_child
p=argparse.ArgumentParser();p.add_argument('--scratch',type=Path,required=True);p.add_argument('--work',action='append');a=p.parse_args();a.scratch.mkdir(parents=True,exist_ok=False)
e=ROOT/'2026_june_14_make_it_consistent/step21_evidence';records=[]
if a.work and (e/'build.json').exists():records=[r for r in json.loads((e/'build.json').read_text()) if r['work'] not in a.work]
for name,source,assets,repo in CASES:
 if a.work and name not in a.work:continue
 w=a.scratch/name;w.mkdir();shutil.copyfile(ROOT/source,w/'main.tex')
 for folder,names in assets.items():
  for f in names:shutil.copyfile(ROOT/folder/f,w/f)
 records.append({'work':name,'source':source,'repo':repo,'tools':[]});r=records[-1]
 for n in (1,2):r['tools'].append(run_child(['/Library/TeX/texbin/pdflatex','-interaction=nonstopmode','-halt-on-error','main.tex'],w,w/f'pass{n}.txt'))
 log=(w/'main.log').read_text()
 if any(t in log for t in ('Overfull \\hbox','Overfull \\vbox','undefined references','undefined on input line')):raise RuntimeError('Layout/reference failure: '+name)
 out=ROOT/pdf(name);shutil.copyfile(w/'main.pdf',out);shutil.copyfile(w/'main.log',e/f'{name}_build.log')
 render=w/'render';render.mkdir()
 r['tools'].append(run_child(['/opt/homebrew/bin/pdftoppm','-r','85','-png',str(out),str(render/'page')],w,w/'render.log'))
 r['tools'].append(run_child(['/opt/homebrew/bin/pdftotext',str(out),str(e/f'{name}_text.txt')],w,w/'text.log'))
 r.update(pdf=pdf(name),pages=len(list(render.glob('*.png'))),render_dir=str(render))
 (e/'build.json').write_text(json.dumps(records,indent=2)+'\n');print(name,r['pages'],'pages',flush=True)
