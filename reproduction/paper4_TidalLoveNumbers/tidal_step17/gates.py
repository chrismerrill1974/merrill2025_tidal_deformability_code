"""Numerical admission gates; preserves results and fails the process on failure."""
from pathlib import Path
import json,argparse
import numpy as np
from scipy.optimize import brentq
from unittest.mock import patch
import tides
from tides import EOS,C,love_k2,solve,fixed,pair
class Polytrope:
 analytic=True;P_min=1e15;P_max=1e38
 def values(self,P):
  if not self.P_min*(1-1e-12)<=P<=self.P_max*(1+1e-12):raise ValueError('polytrope range')
  rho=np.sqrt(P/4e4);return rho*C*C+P,rho
 def eps_of_P(self,P):return self.values(P)[0]
 def deps_dP(self,P):return C*C/(2*np.sqrt(4e4*P))+1

def main():
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
 result={'status':'running'}
 def save():(a.output/'gates.json').write_text(json.dumps(result,indent=2)+'\n')
 save()
 try:
  # Known Newtonian incompressible Love number; sensitive to cancellation.
  k=love_k2(1e-7,-1.);assert abs(k-.75)<1e-6
  result['incompressible_extraction']={'C':1e-7,'y_out':-1.,'k2':k,'limit':.75}
  rejects=[]
  for kw in [dict(Pc=1e10),dict(Pc=2e35,A=-.1),dict(Pc=2e35,form='bad')]:
   try:solve(EOS(),**kw)
   except ValueError:rejects.append(kw)
   else:raise AssertionError('Invalid input accepted')
  from types import SimpleNamespace
  with patch.object(tides,'solve_ivp',return_value=SimpleNamespace(success=False,message='injected failure')):
   try:solve(EOS(),2e35)
   except RuntimeError:rejects.append('injected integration failure')
   else:raise AssertionError('Failed integration accepted')
  result['rejected_inputs']=rejects;save()
  pol=Polytrope();stars=[]
  for target in [.002,.004,.008]:
   lp=brentq(lambda lp:solve(pol,np.exp(lp),rtol=2e-10,atol=1e-19,max_step=.05)['C']-target,np.log(1e27),np.log(1e33),xtol=2e-10)
   star=solve(pol,np.exp(lp),rtol=2e-10,atol=1e-19,max_step=.05);stars.append(star)
  extrap=float(np.polyfit([v['C'] for v in stars],[v['k2'] for v in stars],2)[2]);analytic=(15-np.pi**2)/(2*np.pi**2)
  assert abs(extrap-analytic)<5e-5
  result['n1_newtonian']={'stars':stars,'extrapolated_k2':extrap,'analytic':analytic,'tolerance':5e-5};save()
  gr=fixed(EOS(),1.4,0.,.12);g2=fixed(EOS(),1.4,0.,.5,'emb2');hy=fixed(EOS(),1.4,0.,.12,form='H')
  assert 270<gr['star']['Lambda']<330
  assert abs(g2['star']['Lambda']/gr['star']['Lambda']-1)<1e-6
  assert abs(hy['star']['Lambda']/gr['star']['Lambda']-1)<1e-6
  result['GR']={'primary':gr,'second_embedding':g2,'linear_H':hy};save()
  # EOS energy derivative agrees with centered finite differences inside intervals.
  eos=EOS();err=[]
  for P in np.sqrt(eos.P[:-1]*eos.P[1:])[::10]:
   fd=(eos.eps_of_P(P*1.000001)-eos.eps_of_P(P*.999999))/(P*.000002)
   err.append(abs(fd/tides.dedp(eos,P)-1))
  assert max(err)<1e-7
  result['derivative_max_relative_difference']=max(err)
  result['status']='complete';save();print(json.dumps({'status':'complete','n1_k2_limit':extrap,'GR_Lambda14':gr['star']['Lambda'],'GR_k2':gr['star']['k2']}),flush=True)
 except BaseException as e:result.update(status='failed',error=repr(e));save();raise
if __name__=='__main__':main()
