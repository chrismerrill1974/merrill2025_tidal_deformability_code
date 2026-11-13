"""Build isolated working PDFs from the retained I/II sources and step 12 figures.
Numerical outputs must already exist. Does not touch historical PDFs or aux files.
"""
from pathlib import Path
import argparse,json,shutil,subprocess,tempfile
ROOT=Path(__file__).resolve().parents[1]

def build(scratch,output):
    scratch=Path(scratch);scratch.mkdir(parents=True,exist_ok=False)
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    cases=[('paper1','paper1_SolarCore/V3_Paper/SolarCore_v3.tex','paper1_SolarCore/solar_step12/figures',['score_map.pdf','slice_endpoints.pdf','slice_table.tex']),('paper2','paper2_WhiteDwarfs/V5_Paper/main.tex','paper2_WhiteDwarfs/wd_step12/figures',['radius_shifts.pdf','boundary_table.tex'])]
    records=[]
    for name,source,figdir,files in cases:
        work=scratch/name;work.mkdir();shutil.copyfile(ROOT/source,work/'main.tex')
        for f in files:shutil.copyfile(ROOT/figdir/f,work/f)
        for passno in [1,2]:
            result=subprocess.run(['pdflatex','-interaction=nonstopmode','-halt-on-error','main.tex'],cwd=work,capture_output=True,text=True)
            (work/f'pass{passno}.txt').write_text(result.stdout+result.stderr)
            if result.returncode:raise RuntimeError(f'{name} LaTeX failed; see {work}')
        pdf=output/f'rdt_{name}_step12_review_2026-11-11.pdf';shutil.copyfile(work/'main.pdf',pdf)
        evidence=ROOT/'2026_june_14_make_it_consistent/step12_evidence';evidence.mkdir(exist_ok=True)
        shutil.copyfile(work/'main.log',evidence/f'{name}_build.log')
        render=work/'render';render.mkdir()
        subprocess.run(['pdftoppm','-r','110','-png',str(pdf.resolve()),str(render/'page')],check=True)
        records.append({'paper':name,'source':source,'figures':[figdir+'/'+f for f in files],'pdf':str(pdf),'render_dir':str(render)})
        print(pdf,render,flush=True)
    (scratch/'build.json').write_text(json.dumps(records,indent=2)+'\n')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--scratch',required=True);p.add_argument('--output',default=str(ROOT/'output/pdf'));a=p.parse_args();build(a.scratch,a.output)
