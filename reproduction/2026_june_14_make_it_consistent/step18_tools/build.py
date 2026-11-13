"""Isolated sequential TeX builds and renders with throttled child tools."""
from pathlib import Path
import argparse,json,os,resource,shutil,signal,subprocess,time
from resource_control import require_controls
ROOT=Path(__file__).resolve().parents[2]
EVID=ROOT/'2026_june_14_make_it_consistent/step18_evidence'
def run_child(cmd,cwd,log):
 require_controls();start=time.monotonic();before=resource.getrusage(resource.RUSAGE_CHILDREN)
 with log.open('w') as f:
  p=subprocess.Popen(cmd,cwd=cwd,stdout=f,stderr=subprocess.STDOUT)
  try:
   ni=os.getpriority(os.PRIO_PROCESS,p.pid)
   if ni!=19:raise RuntimeError('Child priority not 19')
   # 94% wall-duty cap leaves scheduling margin below the authorized 95%.
   while p.poll() is None:
    if time.monotonic()-start>180:raise RuntimeError('Tool deadline exceeded')
    time.sleep(.094)
    if p.poll() is not None:break
    os.kill(p.pid,signal.SIGSTOP);time.sleep(.006)
    if p.poll() is None:os.kill(p.pid,signal.SIGCONT)
   if p.returncode:raise RuntimeError(f'Command failed: {cmd}; {log}')
  finally:
   if p.poll() is None:
    os.kill(p.pid,signal.SIGCONT);p.terminate();p.wait(timeout=5)
 after=resource.getrusage(resource.RUSAGE_CHILDREN);wall=time.monotonic()-start;cpu=after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime
 return {'command':cmd,'pid':p.pid,'nice':ni,'returncode':p.returncode,'wall_seconds':wall,'cpu_seconds':cpu,'average_CPU_fraction':cpu/wall,'wall_duty_target':.94,'one_child_at_a_time':True}
def main():
 require_controls();a=argparse.ArgumentParser();a.add_argument('--scratch',type=Path,required=True);a.add_argument('--paper',choices=['paper3','paper4']);opt=a.parse_args();opt.scratch.mkdir(parents=True,exist_ok=False)
 out=ROOT/'output/pdf';out.mkdir(exist_ok=True,parents=True);records=[]
 cases=[('paper3','paper3_NeutronStars/V8_Paper','paper3_v8.tex',['shared_results.tex','common_point.tex','ns_extrema.tex','ns_parameter_dependence.pdf']),('paper4','paper4_TidalLoveNumbers/V6_Paper','paper4_v6.tex',['shared_results.tex','common_point.tex','tidal_samples.tex','tidal_response.pdf'])]
 if opt.paper:
  cases=[c for c in cases if c[0]==opt.paper]
  if (EVID/'build.json').exists():records=[r for r in json.loads((EVID/'build.json').read_text()) if r['paper']!=opt.paper]
 for name,folder,source,assets in cases:
  work=opt.scratch/name;work.mkdir();shutil.copyfile(ROOT/folder/source,work/'main.tex')
  for f in assets:shutil.copyfile(ROOT/folder/f,work/f)
  record={'paper':name,'source':folder+'/'+source,'tools':[]}
  for n in [1,2]:record['tools'].append(run_child(['/Library/TeX/texbin/pdflatex','-interaction=nonstopmode','-halt-on-error','main.tex'],work,work/f'pass{n}.txt'))
  text=(work/'main.log').read_text()
  if 'Overfull \\hbox' in text or 'undefined references' in text:raise RuntimeError('Unresolved layout/references in '+name)
  pdf=out/f'rdt_{name}_step18_review_2026-11-11.pdf';shutil.copyfile(work/'main.pdf',pdf);shutil.copyfile(work/'main.log',EVID/f'{name}_build.log')
  render=work/'render';render.mkdir()
  record['tools'].append(run_child(['/opt/homebrew/bin/pdftoppm','-r','110','-png',str(pdf),str(render/'page')],work,work/'render.log'))
  record.update(pdf=str(pdf.relative_to(ROOT)),render_dir=str(render),pages=len(list(render.glob('page-*.png'))))
  records.append(record);(EVID/'build.json').write_text(json.dumps(records,indent=2)+'\n');print(name,record['pages'],'pages',pdf,flush=True)
if __name__=='__main__':main()
