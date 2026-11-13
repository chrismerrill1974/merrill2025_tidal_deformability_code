"""Replay the archived V6 1.4-solar-mass example without modifying its files."""
from pathlib import Path
import sys,json,hashlib
from resource_control import require_controls
ROOT=Path(__file__).resolve().parents[2]
def main():
 require_controls();sys.path.insert(0,str(ROOT/'paper4_TidalLoveNumbers/tidal_v6'))
 from operator_tidal import at_mass,CompOSE_SLy4
 eos=CompOSE_SLy4();gr=at_mass(eos,1.4,0.,.122);modified=at_mass(eos,1.4,.0087,.122)
 import numpy as np
 assert np.all(np.isfinite([gr,modified])) and max(abs(v[0]-1.4) for v in [gr,modified])<1e-6
 result={'status':'complete','scope':'Historical V6 pipeline replay, not admission of its numerical/matching assumptions','columns':['M','R','C','k2','Lambda'],'A':.0087,'alpha':.122,'GR':gr,'modified':modified,'fractional_Lambda_shift':modified[-1]/gr[-1]-1}
 out=Path(__file__).resolve().parent/'legacy_replay.json';assert not out.exists();out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
if __name__=='__main__':main()
