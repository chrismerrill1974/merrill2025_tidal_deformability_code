"""Retained one-process, nice-19, single-core resource launcher."""
from pathlib import Path
import argparse,json,runpy,signal,sys
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'paper3_NeutronStars/ns_step15'))
import resource_control as rc
p=argparse.ArgumentParser();p.add_argument('--log',type=Path,required=True);p.add_argument('script');p.add_argument('args',nargs=argparse.REMAINDER);a=p.parse_args();status='failed'
try:
 rc.start();a.log.parent.mkdir(parents=True,exist_ok=True);a.log.write_text(json.dumps({'status':'running',**rc.measurement()},indent=2)+'\n')
 script=Path(a.script).resolve()
 if script.parent!=Path(__file__).resolve().parent:raise ValueError('step20 scripts only')
 sys.path.insert(0,str(script.parent));sys.argv=[str(script),*a.args];runpy.run_path(str(script),run_name='__main__');status='complete'
finally:
 if rc.ACTIVE:a.log.write_text(json.dumps({'status':status,**rc.measurement()},indent=2)+'\n')
 signal.setitimer(signal.ITIMER_PROF,0.)
