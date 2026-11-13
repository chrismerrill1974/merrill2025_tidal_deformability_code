"""Resolve the near-edge turnover that the initial coarse GR grid missed."""
from pathlib import Path
import argparse,hashlib,json
from model import ROOT,EOS,solve

def run(output):
    e=EOS();rows=[]
    for emb,A in [('emb1',0.),('emb1',.0087),('emb2',.0087)]:
        for f in [.98,.99,.995,.999,1.]:
            v=solve(e,e.P_max*f,A,.12,emb,rtol=2e-12,atol=1e-20,max_step=.02)
            rows.append({'embedding':emb,'A':A,'alpha':.12,'pressure_fraction':f,'result':v})
    gr=[r['result']['M'] for r in rows if r['A']==0]
    if not max(gr[1:-1])>max(gr[0],gr[-1]):raise RuntimeError('GR near-edge turnover not reproduced')
    for emb in ['emb1','emb2']:
        ms=[r['result']['M'] for r in rows if r['A']>0 and r['embedding']==emb]
        if not all(b>a for a,b in zip(ms,ms[1:])):raise RuntimeError('modified endpoint classification changed')
    files=[Path(__file__),Path(__file__).with_name('model.py'),ROOT/'paper3_NeutronStars/ns_v7/eos.nb',ROOT/'paper3_NeutronStars/ns_v7/eos.thermo']
    result={'status':'pass','settings':{'rtol':2e-12,'atol':1e-20,'max_step':.02},'rows':rows,'finding':'GR has a shallow interior turnover; both nonzero examples increase through the checked near-edge points','source_sha256':{str(p.resolve().relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}}
    output.write_text(json.dumps(result,indent=2)+'\n')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();run(a.output)
