"""Verify shared III/IV evidence, generated assets and compiled PDF content."""
from pathlib import Path
import json,hashlib,re,subprocess
from build import run_child,ROOT,EVID
from resource_control import require_controls

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 require_controls();binding=json.loads((EVID/'asset_bindings.json').read_text());assert binding['status']=='pass'
 for p,h in binding['input_manifest_sha256'].items():assert sha(ROOT/p)==h,p
 for p,h in binding['asset_sha256'].items():assert sha(EVID/'assets'/p)==h,p
 builds=json.loads((EVID/'build.json').read_text());assert {b['paper'] for b in builds}=={'paper3','paper4'}
 sources=[ROOT/b['source'] for b in builds]
 assert (sources[0].parent/'shared_results.tex').read_bytes()==(sources[1].parent/'shared_results.tex').read_bytes()
 assert (sources[0].parent/'common_point.tex').read_bytes()==(sources[1].parent/'common_point.tex').read_bytes()
 records=[]
 for b in builds:
  source=ROOT/b['source'];s=source.read_text();assert 'Astra (OpenAI)' in s and 'prescribed Schwarzschild' in s
  assert '0.122, 0.0087' not in s and 'fig_mr_v8.pdf' not in s and 'fig_lambda_v6.pdf' not in s
  for name in re.findall(r'\\input\{([^}]+)\}',s):assert sha(source.parent/name)==binding['asset_sha256'][name]
  for name in re.findall(r'\\includegraphics(?:\[[^]]*\])?\{([^}]+)\}',s):assert sha(source.parent/name)==binding['asset_sha256'][name]
  keys=set(re.findall(r'\\bibitem\{([^}]+)\}',s))
  for refs in re.findall(r'\\cite\{([^}]+)\}',s):assert set(refs.split(','))<=keys
  log=(EVID/(b['paper']+'_build.log')).read_text();assert 'Overfull \\hbox' not in log and 'undefined references' not in log
  txt=EVID/(b['paper']+'_text.txt');tool=run_child(['/opt/homebrew/bin/pdftotext','-layout',str(ROOT/b['pdf']),str(txt)],ROOT,EVID/(b['paper']+'_extract.log'))
  raw=txt.read_text();normalized=raw.replace('−','-');flat=' '.join(normalized.split())
  for key in ['SolarA','SolarAlpha','DimensionFloor','GRLambda','TidalShiftOne','TidalShiftTwo','MismatchOne']:
   assert binding['macros'][key] in flat,(b['paper'],key)
  assert 'Astra (OpenAI)' in flat and 'Schwarzschild' in flat and 'GitHub' in flat
  expected=b['pages'];pages=raw.split('\f');assert len([p for p in pages if p.strip()])==expected
  renders=sorted(Path(b['render_dir']).glob('page-*.png'));assert len(renders)==expected
  records.append({'paper':b['paper'],'pages':expected,'source_sha256':sha(source),'pdf_sha256':sha(ROOT/b['pdf']),'text_sha256':sha(txt),'render_sha256':{str(p):sha(p) for p in renders},'extraction_resource':tool})
 result={'status':'pass','papers':records,'shared_inputs_identical':True,'scientific_inputs_unchanged':True,'new_stellar_runs':0}
 (EVID/'verification.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'status':'pass','pages':{b['paper']:b['pages'] for b in builds},'shared_inputs_identical':True,'new_stellar_runs':0}))
if __name__=='__main__':main()
