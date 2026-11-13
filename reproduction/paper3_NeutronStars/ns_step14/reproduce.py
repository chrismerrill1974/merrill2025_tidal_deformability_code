"""Stage pinned step 14 sources/inputs in a fresh directory and rerun the audit."""
from pathlib import Path
import argparse,hashlib,json,os,shutil,subprocess,sys
ROOT=Path(__file__).resolve().parents[2];HERE=Path(__file__).resolve().parent

def main(dest):
    dest=dest.resolve();dest.mkdir(parents=True,exist_ok=False)
    m=json.loads((HERE/'results/manifest.json').read_text())
    if m['status']!='complete':raise ValueError('primary incomplete')
    paths=set(m['contract']['source_sha256'])
    paths.update(str(p.relative_to(ROOT)) for p in HERE.glob('*.py'))
    manifest={}
    for name in sorted(paths):
        source=ROOT/name;h=hashlib.sha256(source.read_bytes()).hexdigest()
        if name in m['contract']['source_sha256'] and h!=m['contract']['source_sha256'][name]:raise ValueError('changed source '+name)
        target=dest/name;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(source,target);manifest[name]=h
    (dest/'source_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    prefix='paper3_NeutronStars/ns_step14';env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1')
    cmds=[['test_model.py'],['run.py','--output',prefix+'/results'],['verify.py','--results',prefix+'/results','--output',prefix+'/verification.json'],['fixed_mass_crosscheck.py','--results',prefix+'/results','--output',prefix+'/fixed_mass_crosscheck.json']]
    cmds.append(['fixed_mass_crosscheck.py','--results',prefix+'/results','--output',prefix+'/fixed_mass_tight.json','--rtol','2e-13','--max-step','500'])
    cmds.append(['endpoint_probe.py','--output',prefix+'/endpoint_probe.json'])
    for i,command in enumerate(cmds):
        cmd=[sys.executable,prefix+'/'+command[0],*command[1:]];print('Running',command[0],flush=True)
        with (dest/f'command_{i+1}.log').open('w') as f:subprocess.run(cmd,cwd=dest,env=env,stdout=f,stderr=subprocess.STDOUT,check=True)
    compared={}
    for name in ['results/evaluations.jsonl','results/summary.json','verification.json','fixed_mass_crosscheck.json','fixed_mass_tight.json','endpoint_probe.json']:
        compared[name]=(HERE/name).read_bytes()==(dest/prefix/name).read_bytes()
    if not all(compared.values()):raise RuntimeError('reproduction differs: '+str(compared))
    (dest/'REPRODUCTION_COMPLETE.json').write_text(json.dumps({'status':'pass','byte_identical':compared},indent=2)+'\n');print('Complete',dest)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output-root',type=Path,required=True);a=p.parse_args();main(a.output_root)
