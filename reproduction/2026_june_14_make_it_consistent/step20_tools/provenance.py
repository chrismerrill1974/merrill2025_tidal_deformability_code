"""Record immutable historical input provenance without executing legacy scripts."""
from pathlib import Path
import hashlib,json,subprocess
from resource_control import require_controls
require_controls()
root=Path(__file__).resolve().parents[2];pin='387cccc5da3e2fc582afb9b2ef7087a2ecdd9ad6'
files=['paper6/V1_Paper/paper6_v1.tex','paper6/terrestrial_G_note/note_terrestrial_G_v2.tex','companionpaper_GeometricCompletion/V3_Paper/paper6_v3.tex','paper6/u_sigma.py','paper6/geff_check.py','paper6/g_regression_v2.py','2026_June_14_Revisions/paper6_units_note.md']
hashes={f:hashlib.sha256(subprocess.check_output(['git','show',pin+':'+f],cwd=root)).hexdigest() for f in files}
related=root.parent/'rdt_related_repos/merrill2025_neutronstar_code'
external={'repo':subprocess.check_output(['git','remote','get-url','origin'],cwd=related,text=True).strip(),'head':subprocess.check_output(['git','rev-parse','HEAD'],cwd=related,text=True).strip(),'files':{}}
for f in ['companion_alpha1_bound.py','solar_model.py']:
 external['files'][f]=hashlib.sha256((related/f).read_bytes()).hexdigest()
external['path_status']=subprocess.check_output(['git','status','--short','--','companion_alpha1_bound.py','solar_model.py'],cwd=related,text=True)
out={'input_commit':pin,'historical_sha256':hashes,'alpha1_legacy':external,'alpha1_recomputed':False,'historical_scientific_scripts_modified':False}
(root/'2026_june_14_make_it_consistent/step20_evidence/input_provenance.json').write_text(json.dumps(out,indent=2)+'\n')
print('Input provenance retained')
