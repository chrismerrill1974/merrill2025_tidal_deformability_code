"""Reproduce step 12 in a NEW directory using pinned step 11 solar inputs.
Retains the source/data closure needed for GitHub publication packaging.
"""
from pathlib import Path
import argparse,hashlib,json,os,shutil,subprocess,sys
ROOT=Path(__file__).resolve().parents[1]
BASE=Path('2026_june_14_make_it_consistent')

def stage(destination):
    destination=Path(destination).resolve();destination.mkdir(parents=True,exist_ok=False)
    prior=json.loads((ROOT/BASE/'STEP11_ARTIFACTS.json').read_text())
    files=set(prior['sha256'])
    for f,h in prior['sha256'].items():
        if hashlib.sha256((ROOT/f).read_bytes()).hexdigest()!=h:raise ValueError('Changed step 11 input: '+f)
    files.update([str(BASE/'STEP11_ARTIFACTS.json'),'paper2_WhiteDwarfs/audit_2026_11_11/audit.py','paper2_WhiteDwarfs/wd_v5/model.py','paper2_WhiteDwarfs/wd_v5/validate.py','paper1_SolarCore/solar_step12/figures.py','paper2_WhiteDwarfs/wd_step12/propagate.py','paper2_WhiteDwarfs/wd_step12/figures.py','paper1_SolarCore/V3_Paper/SolarCore_v3.tex','paper2_WhiteDwarfs/V5_Paper/main.tex',str(BASE/'build_step12_papers.py'),str(BASE/'verify_step12.py'),str(BASE/'reproduce_step12.py'),str(BASE/'step12_requirements.txt')])
    manifest={}
    for f in sorted(files):
        target=destination/f;target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(ROOT/f,target);manifest[f]=hashlib.sha256(target.read_bytes()).hexdigest()
    (destination/'step12_source_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    return destination

def run(destination):
    env={**os.environ,'PYTHONDONTWRITEBYTECODE':'1','MPLCONFIGDIR':str(destination/'tmp/mpl')}
    logs=destination/'step12_logs';logs.mkdir()
    commands=[['paper1_SolarCore/solar_step12/figures.py','--output','paper1_SolarCore/solar_step12/figures'],['paper2_WhiteDwarfs/wd_step12/propagate.py','--output','paper2_WhiteDwarfs/wd_step12/results'],['paper2_WhiteDwarfs/wd_step12/figures.py','--output','paper2_WhiteDwarfs/wd_step12/figures'],[str(BASE/'build_step12_papers.py'),'--scratch',str(destination/'tmp/build')],[str(BASE/'verify_step12.py')]]
    for i,args in enumerate(commands,1):
        print('Running',args,flush=True)
        with (logs/f'{i:02d}.txt').open('w') as f:subprocess.run([sys.executable,*args],cwd=destination,env=env,stdout=f,stderr=subprocess.STDOUT,check=True)
    (destination/'STEP12_REPRODUCTION_COMPLETE.json').write_text(json.dumps({'status':'complete','commands':commands,'solar_inputs':'pinned step 11 artifacts, not rerun here'},indent=2)+'\n')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output-root',required=True);p.add_argument('--stage-only',action='store_true');a=p.parse_args();d=stage(a.output_root)
    if not a.stage_only:run(d)
    print(d)
