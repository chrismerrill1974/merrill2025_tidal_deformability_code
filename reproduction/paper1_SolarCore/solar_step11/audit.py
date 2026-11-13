"""Numerical sensitivity evidence; failures retained, never replaced by zero."""
from pathlib import Path
import json,sys
from dataclasses import replace,asdict
import numpy as np
from scipy.integrate import simpson
from scipy.optimize import brentq,minimize_scalar
from solver import Settings,response,calibrate
from inference import score
from numerical_guardrails import fingerprint,sha
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'solarcore_v3'))
import model as legacy

def run(output):
    out=Path(output);out.mkdir(parents=True,exist_ok=False)
    settings={'baseline':Settings(),'tight':replace(Settings(),rtol=2e-13,atol=2e-15,max_step=.02),
      'start':replace(Settings(),xi0=1e-7),'profile':replace(Settings(),profile_n=3201),
      'si':replace(Settings(),si=True),'floor_1e-6':replace(Settings(),floor=1e-6),
      'floor_1e-10':replace(Settings(),floor=1e-10),'loose':replace(Settings(),rtol=1e-9,atol=1e-11,max_step=.05)}
    record={'contract':fingerprint(HERE,['solver.py','inference.py','audit.py'],{k:asdict(v) for k,v in settings.items()}),'points':[],'legacy_comparison':[],'endpoint_checks':[]}
    standard={k:calibrate(0,.1,v) for k,v in settings.items()}
    for alpha in [.02,.1,.15,.5]:
      for A in [.0042,.009,.04]:
        entry={'A':A,'alpha':alpha,'settings':{}}
        for name,s in settings.items():
          try:entry['settings'][name]={'status':'complete','response':response(A,alpha,s,standard[name])}
          except Exception as e:entry['settings'][name]={'status':'failed','error':str(e)}
        record['points'].append(entry)
    for A,alpha in [(.0042,.02),(.0334,.15),(.0087,.15)]:
      item={'A':A,'alpha':alpha}
      for label,kw in [('scan',dict(rtol=1e-9,atol=1e-11,max_step=.05)),('default',{})]:
        try:
          st=legacy.calibrate(legacy.Params(0,.1),**kw);m=legacy.calibrate(legacy.Params(A,alpha),**kw)
          _,h=legacy.sound_speed_deviation(m,st)
          item[label]={'status':'complete','rho_c':m.rho_c,'xi1':m.xi1,'ratios':legacy.flux_ratios(m,st),'h':float(max(abs(h)))}
        except Exception as e:item[label]={'status':'failed','error':str(e)}
      for label in ['loose','baseline']:
        m=calibrate(A,alpha,settings[label]);sol=m['sol'];quadrature=sol.y[2,-1]
        item['augmented_'+label]={'iterations':m['iterations'],'residual':m['residual'],'adaptive_simpson_mass_relative_difference':float(simpson(sol.t**2*np.clip(sol.y[0],0,None)**3,x=sol.t)/quadrature-1)}
      record['legacy_comparison'].append(item)
    # Tight/cutoff/profile/start effects on actual nominal endpoints at selected alphas.
    for alpha in np.geomspace(.02,.5,9):
      item={'alpha':float(alpha),'settings':{}}
      for name in ['baseline','tight','start','profile','floor_1e-6','floor_1e-10','si']:
        s=settings[name]
        try:
          f=lambda A:score(response(A,alpha,s,standard[name]),'source_precision_diagonal_proxy_v1')
          opt=minimize_scalar(f,bounds=(0,.001),method='bounded',options={'xatol':2e-10})
          if not opt.success:raise RuntimeError('Minimum failure')
          qm=min(f(0),float(opt.fun))
          end=brentq(lambda A:f(A)-qm-2.71,.001,.025,xtol=2e-10)
          item['settings'][name]={'status':'complete','endpoint':float(end)}
        except Exception as e:item['settings'][name]={'status':'failed','error':str(e)}
      record['endpoint_checks'].append(item)
      print('audited endpoint',alpha,flush=True)
    (out/'audit.json').write_text(json.dumps(record,indent=2)+'\n')

if __name__=='__main__':run(sys.argv[1])
