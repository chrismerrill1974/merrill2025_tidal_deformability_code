"""Independent radial-field fixed-mass solutions, with GR paired on same path."""
from pathlib import Path
import argparse,hashlib,json
from scipy.optimize import brentq
import numpy as np
from model import ROOT,EOS
from direct import solve_direct

def run(results,output,rtol=2e-11,max_step=1000.):
    m=json.loads((results/'manifest.json').read_text())
    if m['status']!='complete':raise ValueError('incomplete primary')
    for f,h in m['contract']['source_sha256'].items():
        if hashlib.sha256((ROOT/f).read_bytes()).hexdigest()!=h:raise ValueError('source mismatch')
    if hashlib.sha256((results/'summary.json').read_bytes()).hexdigest()!=m['sha256']['summary.json']:raise ValueError('summary hash')
    primary=json.loads((results/'summary.json').read_text())['configurations']['tight']['results']
    e=EOS();rows=[];sols={};settings={'rtol':rtol,'atol':rtol*5e-9,'max_step':max_step,'derivative':'exact'}
    for label,A,emb in [('GR',0.,'emb1'),('emb1',.0087,'emb1'),('emb2',.0087,'emb2')]:
        def ev(lp):
            v=solve_direct(e,np.exp(lp),A,.12,emb,**settings);rows.append({'label':label,'A':A,'alpha':.12,'embedding':emb,'result':v});return v
        lp=brentq(lambda p:ev(p)['M']-1.4,np.log(1e35),np.log(2e35),xtol=2e-10)
        sols[label]=ev(lp)
        if abs(sols[label]['M']-1.4)>1e-7:raise RuntimeError('fixed-mass residual')
    diffs={}
    for emb in ['emb1','emb2']:
        delta=sols[emb]['R']/sols['GR']['R']-1
        diffs[emb]={'direct_delta_R14_fraction':delta,'flux_delta_R14_fraction':primary[emb]['delta_R14_fraction'],'difference':delta-primary[emb]['delta_R14_fraction']}
        if abs(diffs[emb]['difference'])>1e-6:raise RuntimeError('independent fixed-mass mismatch')
    paths=[Path(__file__),Path(__file__).with_name('direct.py'),Path(__file__).with_name('model.py')]
    result={'status':'pass','settings':settings,'solutions':sols,'paired_comparison':diffs,'evaluations':rows,'source_sha256':{str(p.resolve().relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},'primary_summary_sha256':m['sha256']['summary.json']}
    output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'paired_comparison':diffs,'evaluations':len(rows)},indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--results',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--rtol',type=float,default=2e-11);p.add_argument('--max-step',type=float,default=1000.);a=p.parse_args();run(a.results,a.output,a.rtol,a.max_step)
