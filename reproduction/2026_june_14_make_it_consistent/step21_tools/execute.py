"""Controlled single scientific process used by the publication reproducer."""
from pathlib import Path
import json,os,runpy,shutil,signal,sys
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'paper3_NeutronStars/ns_step15'))
import resource_control as rc
rc.start();os.environ['MPLCONFIGDIR']=str(Path('/tmp')/f'rdt-mpl-{os.getpid()}');os.environ['PYTHONDONTWRITEBYTECODE']='1'
try:
 mode=sys.argv[1]
 if mode=='script':
  path=ROOT/sys.argv[2];assert path.resolve().is_relative_to(ROOT)
  sys.path.insert(0,str(path.parent));sys.path.insert(0,str(ROOT));sys.argv=[str(path),*sys.argv[3:]];runpy.run_path(str(path),run_name='__main__')
 elif mode=='shared':
  out=Path(sys.argv[2]);out.mkdir();path=ROOT/'2026_june_14_make_it_consistent/step18_tools/assets.py'
  ns=runpy.run_path(str(path),run_name='asset_module');g=ns['main'].__globals__;g['OUT']=out
  for key in ['P3','P4']:g[key]=out/key;g[key].mkdir()
  ns['main']()
 elif mode=='tidal':
  out=Path(sys.argv[2]);(out/'results').mkdir(parents=True);here=ROOT/'paper4_TidalLoveNumbers/tidal_step17'
  for n in ['results/summary.json','verification.json']:shutil.copyfile(here/n,out/n)
  ns=runpy.run_path(str(here/'plot.py'),run_name='plot_module');ns['main'].__globals__['HERE']=out;ns['main']()
 elif mode=='algebra':
  out=Path(sys.argv[2]);out.mkdir()
  for n in ['2026_june_14_make_it_consistent/step20_tools/audit.py','2026_june_14_make_it_consistent/step20_evidence/nist_tables.json','paper6/geff_check.py']:
   target=out/n;target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(ROOT/n,target)
  (out/'paper6/terrestrial_G_note').mkdir(parents=True)
  runpy.run_path(str(out/'2026_june_14_make_it_consistent/step20_tools/audit.py'),run_name='__main__')
 else:raise ValueError(mode)
finally:signal.setitimer(signal.ITIMER_PROF,0.)
