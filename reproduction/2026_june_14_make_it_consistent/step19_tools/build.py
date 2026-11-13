"""Rebuild Paper V from retained sources with sequential controlled tools."""
from pathlib import Path
import argparse, hashlib, json, shutil, sys
from resource_control import require_controls
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'2026_june_14_make_it_consistent/step18_tools'))
from build import run_child
require_controls()
p=argparse.ArgumentParser();p.add_argument('--scratch',type=Path,required=True);a=p.parse_args()
a.scratch.mkdir(parents=True,exist_ok=False)
evid=ROOT/'2026_june_14_make_it_consistent/step19_evidence';evid.mkdir(exist_ok=True)
folder=ROOT/'paper5_BlackHoles/V3_Paper'
shared=(folder/'shared_results.tex').read_bytes()
assert shared==(ROOT/'paper3_NeutronStars/V8_Paper/shared_results.tex').read_bytes()==(ROOT/'paper4_TidalLoveNumbers/V6_Paper/shared_results.tex').read_bytes()
shutil.copyfile(folder/'paper5_v3.tex',a.scratch/'main.tex');shutil.copyfile(folder/'shared_results.tex',a.scratch/'shared_results.tex')
records=[]
for n in (1,2):records.append(run_child(['/Library/TeX/texbin/pdflatex','-interaction=nonstopmode','-halt-on-error','main.tex'],a.scratch,a.scratch/f'pass{n}.txt'))
log=(a.scratch/'main.log').read_text()
assert 'Overfull \\hbox' not in log and 'undefined references' not in log and 'undefined on input line' not in log
pdf=ROOT/'output/pdf/rdt_paper5_step19_review_2026-11-12.pdf';shutil.copyfile(a.scratch/'main.pdf',pdf)
shutil.copyfile(a.scratch/'main.log',evid/'paper5_build.log')
render=a.scratch/'render';render.mkdir()
records.append(run_child(['/opt/homebrew/bin/pdftoppm','-r','110','-png',str(pdf),str(render/'page')],a.scratch,a.scratch/'render.log'))
records.append(run_child(['/opt/homebrew/bin/pdftotext',str(pdf),str(evid/'paper5_text.txt')],a.scratch,a.scratch/'extract.log'))
(evid/'build.json').write_text(json.dumps({'tools':records,'pages':len(list(render.glob('*.png'))),'render_dir':str(render),'shared_macros_equal_III_IV':True,'shared_sha256':hashlib.sha256(shared).hexdigest(),'layout_reference_checks':'pass'},indent=2)+'\n')
print(render,flush=True)
