"""Retain EOS-scope review evidence without rerunning stellar simulations."""
from pathlib import Path
import hashlib,json,re,sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'paper3_NeutronStars/ns_step15'))
import resource_control
MANUSCRIPTS=[
'paper1_SolarCore/V3_Paper/SolarCore_v3.tex',
'paper2_WhiteDwarfs/V5_Paper/main.tex',
'paper3_NeutronStars/V8_Paper/paper3_v8.tex',
'paper4_TidalLoveNumbers/V6_Paper/paper4_v6.tex',
'paper5_BlackHoles/V3_Paper/paper5_v3.tex',
'companionpaper_GeometricCompletion/V3_Paper/paper6_v3.tex',
'paper6/V1_Paper/paper6_v1.tex',
'paper6/terrestrial_G_note/note_terrestrial_G_v2.tex']
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 resource_control.start()
 out=Path(__file__).parent/'step16_evidence';out.mkdir(exist_ok=True)
 evidence={'scope':'Static evidence and manuscript-scope audit; no new EOS or stellar calculation','manuscripts':{},'scientific_manifests':{}}
 for name in MANUSCRIPTS:
  s=(ROOT/name).read_text();assert 'Astra (OpenAI)' in s,name
  # Syntax smoke check for inserted prose; not a substitute for PDF compilation.
  stripped=re.sub(r'(?<!\\)%.*','',s)
  braces=re.findall(r'(?<!\\)[{}]',stripped);depth=0
  for b in braces:
   depth+=1 if b=='{' else -1
   assert depth>=0,name
  assert depth==0,name
  hits=[{'line':i,'text':line} for i,line in enumerate(s.splitlines(),1) if re.search(r'EOS|equation.of.state|APR|SLy|saturat|univers',line,re.I)]
  evidence['manuscripts'][name]={'sha256':sha(ROOT/name),'scope_review_lines':hits,'balanced_braces':True,'collaborator_entry_present':True}
 for folder in ['paper3_NeutronStars/ns_step14/results','paper3_NeutronStars/ns_step15/resumed95_results','paper3_NeutronStars/ns_step15/refinement','paper3_NeutronStars/ns_step15/analysis','paper3_NeutronStars/ns_step15/checks']:
  d=ROOT/folder;m=json.loads((d/'manifest.json').read_text());assert m['status']=='complete',folder
  sources=m.get('contract',{}).get('source_sha256',m.get('source_sha256',{}))
  for name,h in sources.items():assert sha(ROOT/name)==h,name
  for name,h in m.get('sha256',{}).items():assert sha(d/name)==h,name
  evidence['scientific_manifests'][folder]={'manifest_sha256':sha(d/'manifest.json'),'source_hashes_verified':len(sources),'artifact_hashes_verified':len(m.get('sha256',{}))}
 evidence['resource_controls']=resource_control.measurement()
 evidence['status']='pass'
 (out/'audit.json').write_text(json.dumps(evidence,indent=2)+'\n')
 print(json.dumps({'status':'pass','manuscripts':len(MANUSCRIPTS),'scientific_manifests':len(evidence['scientific_manifests']),'new_stellar_runs':0}))
if __name__=='__main__':
 try:main()
 finally:
  import signal
  signal.setitimer(signal.ITIMER_PROF,0.)
