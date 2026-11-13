"""Stage exact sources in a NEW directory and optionally reproduce step 11.
No historical result is copied or used as a numerical checkpoint.
"""
from pathlib import Path
import argparse,hashlib,json,os,shutil,subprocess,sys
ROOT=Path(__file__).resolve().parents[2]
REL=Path('paper1_SolarCore/solar_step11')

def stage(destination):
    destination=Path(destination).resolve();destination.mkdir(parents=True,exist_ok=False)
    files=list((ROOT/REL).glob('*.py'))+[ROOT/REL/'README.md',ROOT/REL/'requirements.txt',ROOT/'numerical_guardrails.py',ROOT/'tests/test_step11_inference.py',ROOT/'paper1_SolarCore/solarcore_v3/model.py']
    files += [ROOT/'2026_june_14_make_it_consistent'/n for n in ['STEP08_RESPONSE_DECISION.json','STEP09_INFERENCE_CONTRACT.json','STEP09_INPUT_AUDIT.json']]
    hashes={}
    for source in files:
        relative=source.relative_to(ROOT);target=destination/relative;target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(source,target)
        hashes[str(relative)]=hashlib.sha256(target.read_bytes()).hexdigest()
    (destination/'source_manifest.json').write_text(json.dumps(hashes,indent=2)+'\n')
    return destination

def reproduce(destination):
    env={**os.environ,'PYTHONDONTWRITEBYTECODE':'1'};logs=destination/'logs';logs.mkdir()
    commands=[['tests/test_step11_inference.py'],[str(REL/'produce.py'),'--output',str(REL/'results/solar_a')],[str(REL/'reference_run.py'),str(REL/'results/solar_R')],[str(REL/'audit.py'),str(REL/'results/numerical_audit')]]
    commands += [[str(REL/f)] for f in ['local_checks.py','peak.py','peak_numerics.py','summarize.py','verify.py']]
    for index,args in enumerate(commands,1):
        print('Running',args,flush=True)
        with (logs/f'{index:02d}-{Path(args[0]).stem}.txt').open('w') as f:
            subprocess.run([sys.executable,*args],cwd=destination,env=env,stdout=f,stderr=subprocess.STDOUT,check=True)
    (destination/'REPRODUCTION_COMPLETE.json').write_text(json.dumps({'status':'complete','commands':commands},indent=2)+'\n')

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output-root',required=True);p.add_argument('--stage-only',action='store_true',help='Copy and hash sources without running the study')
    a=p.parse_args();destination=stage(a.output_root)
    if not a.stage_only:reproduce(destination)
    print(destination)
