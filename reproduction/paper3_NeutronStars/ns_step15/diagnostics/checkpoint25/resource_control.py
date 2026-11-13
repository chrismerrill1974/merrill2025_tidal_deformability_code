"""One process, verified nice 19, and a 25%-of-one-core CPU-time budget.
ITIMER_PROF accounts process CPU (including native threads); frequent ODE Python
callbacks deliver signals. No accumulated idle credit; each 10 ms CPU quantum
is followed by enough sleep for a <=25% duty cycle. This is a measured soft
process limit, not a claim to control other applications or guarantee temperature.
"""
import argparse,json,os,runpy,signal,sys,time
from pathlib import Path
ACTIVE=False
CPU_FRACTION=.25
QUANTUM=.01
THREAD_VARS=('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS')
STATE={}

def require_controls():
    if not ACTIVE or os.getpriority(os.PRIO_PROCESS,0)!=19 or signal.getitimer(signal.ITIMER_PROF)[1]!=QUANTUM:
        raise RuntimeError('Resource controls required: use resource_control.py; nice 19, one process, CPU throttle')
    if any(os.environ.get(k)!='1' for k in THREAD_VARS):raise RuntimeError('Native thread limit missing')

def start():
    global ACTIVE
    for k in THREAD_VARS:os.environ[k]='1'
    os.setpriority(os.PRIO_PROCESS,0,19)
    if os.getpriority(os.PRIO_PROCESS,0)!=19:raise RuntimeError('Could not establish nice 19')
    wall=time.monotonic();cpu=time.process_time()
    STATE.update(start_wall=wall,start_cpu=cpu,last_wall=wall,last_cpu=cpu,sleeps=0,sleep_seconds=0.,pid=os.getpid())
    def throttle(signum,frame):
        wall=time.monotonic();cpu=time.process_time()
        delay=(cpu-STATE['last_cpu'])/CPU_FRACTION-(wall-STATE['last_wall'])
        if delay>0:
            time.sleep(delay);STATE['sleep_seconds']+=delay;STATE['sleeps']+=1
        STATE['last_wall']=time.monotonic();STATE['last_cpu']=time.process_time()
    signal.signal(signal.SIGPROF,throttle);signal.setitimer(signal.ITIMER_PROF,QUANTUM,QUANTUM)
    ACTIVE=True;require_controls()

def measurement():
    wall=time.monotonic()-STATE['start_wall'];cpu=time.process_time()-STATE['start_cpu']
    return {'pid':os.getpid(),'nice':os.getpriority(os.PRIO_PROCESS,0),'worker_processes':1,'cpu_fraction_target':CPU_FRACTION,'cpu_quantum_seconds':QUANTUM,'wall_seconds':wall,'cpu_seconds':cpu,'measured_cpu_fraction':cpu/wall,'throttle_sleeps':STATE['sleeps'],'throttle_sleep_seconds':STATE['sleep_seconds'],'native_thread_environment':{k:os.environ[k] for k in THREAD_VARS}}

def main():
    p=argparse.ArgumentParser();p.add_argument('--log',type=Path,required=True);p.add_argument('--probe',action='store_true');p.add_argument('script',nargs='?');p.add_argument('args',nargs=argparse.REMAINDER);a=p.parse_args()
    # Register this module before runpy so scientific imports share ACTIVE state.
    sys.modules['resource_control']=sys.modules[__name__]
    start();a.log.parent.mkdir(parents=True,exist_ok=True)
    a.log.write_text(json.dumps({'status':'running',**measurement()},indent=2)+'\n')
    status='failed'
    try:
        if a.probe:
            goal=time.process_time()+2.;x=1
            while time.process_time()<goal:x=(x*1664525+1013904223)%4294967296
            r=measurement()
            if not .15<r['measured_cpu_fraction']<.29 or r['nice']!=19 or r['throttle_sleeps']<50:raise RuntimeError('Resource probe failed')
        else:
            if not a.script:raise ValueError('Script required')
            script=Path(a.script).resolve()
            if script.parent!=Path(__file__).resolve().parent:raise ValueError('Run only step 15 scripts')
            sys.argv=[str(script),*a.args];runpy.run_path(str(script),run_name='__main__')
        status='complete'
    finally:
        a.log.write_text(json.dumps({'status':status,**measurement()},indent=2)+'\n')
        signal.setitimer(signal.ITIMER_PROF,0.)
if __name__=='__main__':main()
