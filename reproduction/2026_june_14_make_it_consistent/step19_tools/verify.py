"""Bind inherited values to the preceding working release; no science rerun."""
from pathlib import Path
import hashlib,json
from resource_control import require_controls
require_controls()
root=Path(__file__).resolve().parents[2]
old=json.loads((root/'2026_june_14_make_it_consistent/STEP18_ARTIFACTS.json').read_text())
checked=[]
for name,want in old['sha256'].items():
 if name.startswith(('paper3_NeutronStars/','paper4_TidalLoveNumbers/','2026_june_14_make_it_consistent/step18_')):
  assert hashlib.sha256((root/name).read_bytes()).hexdigest()==want,name
  checked.append(name)
s=(root/'paper5_BlackHoles/V3_Paper/paper5_v3.tex').read_text()
assert 'A_{95}' not in s and 'tab:ns_profiles' not in s
for item in ('Astra (OpenAI)','\\MismatchOne','\\MismatchTwo','\\LambdaOne','\\LambdaTwo','\\SolarAlpha','\\SolarA','g=C/r^2','alpha=0','arXiv:0806.0016v2','arXiv:2105.07101v2'):
 assert item in s,item
shared=(root/'paper5_BlackHoles/V3_Paper/shared_results.tex').read_bytes()
assert shared==(root/'paper3_NeutronStars/V8_Paper/shared_results.tex').read_bytes()==(root/'paper4_TidalLoveNumbers/V6_Paper/shared_results.tex').read_bytes()
(root/'2026_june_14_make_it_consistent/step19_evidence/verification.json').write_text(json.dumps({'status':'pass','inherited_files_checked':checked,'shared_macros_equal':True,'scope_checks':'pass','new_simulation':False},indent=2)+'\n')
print('Verified',len(checked),'inherited source/asset/evidence hashes and shared macros')
