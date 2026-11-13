"""Sequential three-manuscript build with the retained controlled tool runner."""
from pathlib import Path
import argparse,json,shutil,sys
from resource_control import require_controls
require_controls()
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'2026_june_14_make_it_consistent/step18_tools'))
from build import run_child
p=argparse.ArgumentParser();p.add_argument('--scratch',type=Path,required=True);a=p.parse_args();a.scratch.mkdir(parents=True,exist_ok=False)
e=ROOT/'2026_june_14_make_it_consistent/step20_evidence'
cases=[('paper6','paper6/V1_Paper/paper6_v1.tex',[]),('terrestrial','paper6/terrestrial_G_note/note_terrestrial_G_v2.tex',['material_table.tex','contrast_table.tex']),('companion','companionpaper_GeometricCompletion/V3_Paper/paper6_v3.tex',[])]
records=[]
for name,source,assets in cases:
 w=a.scratch/name;w.mkdir();src=ROOT/source;shutil.copyfile(src,w/'main.tex')
 for asset in assets:shutil.copyfile(src.parent/asset,w/asset)
 tools=[]
 for n in (1,2):tools.append(run_child(['/Library/TeX/texbin/pdflatex','-interaction=nonstopmode','-halt-on-error','main.tex'],w,w/f'pass{n}.txt'))
 log=(w/'main.log').read_text()
 assert not any(t in log for t in ('Overfull \\hbox','undefined references','undefined on input line')),name
 pdf=ROOT/f'output/pdf/rdt_{name}_step20_review_2026-11-12.pdf';shutil.copyfile(w/'main.pdf',pdf);shutil.copyfile(w/'main.log',e/f'{name}_build.log')
 render=w/'render';render.mkdir()
 tools.append(run_child(['/opt/homebrew/bin/pdftoppm','-r','90','-png',str(pdf),str(render/'page')],w,w/'render.log'))
 tools.append(run_child(['/opt/homebrew/bin/pdftotext',str(pdf),str(e/f'{name}_text.txt')],w,w/'extract.log'))
 records.append({'paper':name,'source':source,'tools':tools,'pdf':str(pdf.relative_to(ROOT)),'render_dir':str(render),'pages':len(list(render.glob('*.png'))),'layout_reference_checks':'pass'})
 (e/'build.json').write_text(json.dumps(records,indent=2)+'\n');print(name,records[-1]['pages'],'pages',flush=True)
