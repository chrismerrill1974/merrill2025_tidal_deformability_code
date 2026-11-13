"""Retained bounded step 14 numerical audit; fresh output directory required."""
from pathlib import Path
import argparse, hashlib, json, platform, traceback
import numpy as np
import scipy
from scipy.optimize import brentq, minimize_scalar
from model import EOS,solve,ROOT
from direct import solve_direct

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(out):
    out.mkdir(parents=True,exist_ok=False)
    files=[Path(__file__),Path(__file__).with_name('model.py'),Path(__file__).with_name('direct.py')]
    files += [ROOT/'paper3_NeutronStars/ns_v7'/f for f in ['eos.nb','eos.thermo','eos.mr']]
    files += [ROOT/'paper3_NeutronStars/037/eos_sly4_dh2001_v2.dat']
    contract={'source_sha256':{str(p.relative_to(ROOT)):sha(p) for p in files},'runtime':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__},'reference_sample':{'A':.0087,'alpha':.12,'meaning':'historical canonical sample, not a new solar-domain envelope'},'rstar_cm':1e6}
    manifest={'status':'running','contract':contract}
    def save(): (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    save();journal=(out/'evaluations.jsonl').open('w');counter=0
    def evaluate(eos,pc,A,alpha,emb,cfg,tag,solver=solve):
        nonlocal counter
        row={'id':counter,'tag':tag,'eos':eos.kind,'interpolation':eos.interp,'Pc':float(pc),'A':A,'alpha':alpha,'embedding':emb,'settings':cfg,'solver':solver.__name__};counter+=1
        try:row.update(result=solver(eos,float(pc),A,alpha,emb,**cfg),status='complete')
        except Exception as exc:
            row.update(status='failed',error=repr(exc));journal.write(json.dumps(row)+'\n');journal.flush();raise
        journal.write(json.dumps(row)+'\n');journal.flush();return row['result']
    def stats(e,A,alpha,emb,cfg,tag):
        cache={}
        def ev(lp):
            key=float(lp)
            if key not in cache:cache[key]=evaluate(e,np.exp(key),A,alpha,emb,cfg,tag)
            return cache[key]
        grid=np.unique(np.r_[np.linspace(np.log(5e34),np.log(e.P_max),18),np.log(e.P_max*np.array([.98,.99,.995,.999]))]);curve=[ev(t) for t in grid];i=int(np.argmax([v['M'] for v in curve]))
        status='upper_EOS_boundary' if i==len(grid)-1 else 'interior_refined'
        if i==0:raise RuntimeError('maximum at lower survey boundary')
        peak=curve[i]
        if i<len(grid)-1:
            opt=minimize_scalar(lambda t:-ev(t)['M'],bounds=(grid[i-1],grid[i+1]),method='bounded',options={'xatol':2e-7})
            if not opt.success:raise RuntimeError('maximum optimizer')
            peak=ev(opt.x)
        lo=np.log(2e34);hi=np.log(peak['Pc'])
        if not ev(lo)['M']<1.4<ev(hi)['M']:raise RuntimeError('fixed mass not bracketed')
        lp=brentq(lambda t:ev(t)['M']-1.4,lo,hi,xtol=2e-10)
        r14=ev(lp)
        if abs(r14['M']-1.4)>1e-7:raise RuntimeError('fixed mass residual')
        return {'domain_peak':peak,'peak_status':status,'fixed_1p4':r14,'grid':curve}
    try:
        tight={'rtol':2e-11,'atol':1e-19,'max_step':.04}
        configs={'default':{},'tight':tight,'tighter':dict(rtol=2e-12,atol=1e-20,max_step=.02),
                 'start10':dict(tight,r0=10.),'start1000':dict(tight,r0=1000.),
                 'surface1':dict(tight,surface_factor=1.),'surface10':dict(tight,surface_factor=10.),
                 'surface100':dict(tight,surface_factor=100.),'pchip':tight,
                 'legacy_seed':dict(tight,legacy_seed=True),'RK45':dict(tight,method='RK45')}
        summary={'configurations':{}}
        for name,cfg in configs.items():
            e=EOS(interp='pchip' if name=='pchip' else 'linear');rows={}
            for label,A,emb in [('GR',0.,'emb1'),('emb1',.0087,'emb1'),('emb2',.0087,'emb2')]:
                rows[label]=stats(e,A,.12,emb,cfg,name+'/'+label)
            for emb in ['emb1','emb2']:
                rows[emb]['delta_domain_peak_fraction']=rows[emb]['domain_peak']['M']/rows['GR']['domain_peak']['M']-1
                rows[emb]['delta_R14_fraction']=rows[emb]['fixed_1p4']['R']/rows['GR']['fixed_1p4']['R']-1
            summary['configurations'][name]={'settings':cfg,'results':rows}
            (out/'summary.partial.json').write_text(json.dumps(summary,indent=2)+'\n')
            print(name,[(k,rows[k].get('delta_domain_peak_fraction'),rows[k].get('delta_R14_fraction')) for k in rows],flush=True)
        # Source-convention diagnosis: same pressure/rho table, different epsilon column.
        summary['dh_diagnosis']={}
        for kind in ['dh','dh_bug']:
            summary['dh_diagnosis'][kind]=stats(EOS(kind),0.,.12,'emb1',tight,kind)
        # Direct fields, analytic/legacy-stencil derivatives, and two low/high-alpha probes.
        cross=[]
        for alpha,A in [(.12,0.),(.12,.0087),(.02,.0042),(.5,.0067)]:
            for pc in [2e35,1.3e36]:
                for emb in ['emb1','emb2']:
                    e=EOS();f=evaluate(e,pc,A,alpha,emb,tight,'cross/flux')
                    for deriv in (['exact','fd'] if emb=='emb1' else ['exact']):
                        d=evaluate(e,pc,A,alpha,emb,{'rtol':2e-11,'atol':1e-19,'max_step':1000.,'derivative':deriv},'cross/direct',solve_direct)
                        cross.append({'alpha':alpha,'A':A,'Pc':pc,'embedding':emb,'derivative':deriv,'flux':f,'direct':d,'relative_M_difference':d['M']/f['M']-1,'relative_R_difference':d['R']/f['R']-1})
        summary['direct_comparisons']=cross
        # Paired fixed-Pc corrections compared across independent implementations.
        paired=[]
        for row in cross:
            if row['A']==0:continue
            gr=next(x for x in cross if x['A']==0 and x['Pc']==row['Pc'] and x['embedding']==row['embedding'] and x['derivative']==row['derivative'])
            p={k:row[k] for k in ['alpha','A','Pc','embedding','derivative']}
            for qty in ['M','R']:
                p['flux_delta_'+qty]=row['flux'][qty]/gr['flux'][qty]-1
                p['direct_delta_'+qty]=row['direct'][qty]/gr['direct'][qty]-1
                p['delta_'+qty+'_difference']=p['direct_delta_'+qty]-p['flux_delta_'+qty]
            paired.append(p)
        summary['independent_paired_differences']=paired
        # Historical 60-point maximum was a finite-grid result, not a new bound.
        e=EOS();old=evaluate(e,1.3e36,0.,.12,'emb1',tight,'historical_GR_endpoint')
        ref=np.loadtxt(ROOT/'paper3_NeutronStars/ns_v7/eos.mr');imax=int(np.argmax(ref[:,1]));rr=float(np.interp(1.4,ref[:imax+1,1],ref[:imax+1,0]))
        gr=summary['configurations']['tight']['results']['GR']
        summary['reference_comparison']={'historical_GR_at_1p3e36':old,'reference_tabulated_max_M':float(ref[:,1].max()),'reference_max_at_last_row':imax==len(ref)-1,'reference_R14_linear_interp':rr,'fractional_domain_mass_difference':gr['domain_peak']['M']/ref[:,1].max()-1,'fractional_R14_difference':gr['fixed_1p4']['R']/rr-1}
        # Derivative diagnostic away from knots and spanning knots.
        e=EOS();mids=np.exp((e.lp[:-1]+e.lp[1:])/2);errs=[];knoterrs=[]
        for p in mids:
            fd=(e.rho_of_P(p*1.000001)-e.rho_of_P(p*.999999))/(p*2e-6);errs.append(abs(fd/e.drho_dP(p)-1))
        for p in e.P[1:-1]:
            fd=(e.rho_of_P(p*1.001)-e.rho_of_P(p*.999))/(p*.002);knoterrs.append(abs(fd/e.drho_dP(p)-1))
        summary['EOS_derivative']={'midinterval_1e6_stencil_max_relative_difference':max(errs),'knot_legacy_stencil_vs_right_derivative_max_relative_difference':max(knoterrs),'interpretation':'linear-log interpolant derivative is discontinuous at knots; symmetric stencil mixes adjacent slopes'}
        summary['evaluation_count']=counter
        (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');(out/'summary.partial.json').unlink();journal.close()
        manifest.update(status='complete',evaluations=counter,sha256={p.name:sha(p) for p in out.iterdir() if p.name!='manifest.json'});save()
    except Exception as exc:
        journal.close();manifest.update(status='failed',error=repr(exc),traceback=traceback.format_exc(),evaluations=counter);save();raise
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);args=p.parse_args();main(args.output)
