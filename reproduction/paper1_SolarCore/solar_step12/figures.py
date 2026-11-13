"""Paper I step 12 figures and tables from verified step 11 artifacts.
Writes to a new directory. Original figures and caches remain unchanged.
"""
from pathlib import Path
import argparse,json,sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar
ROOT=Path(__file__).resolve().parents[2];CODE=ROOT/'paper1_SolarCore/solar_step11'
sys.path.insert(0,str(CODE))
from summarize import load_region
from inference import score,VARIANTS,C
from solver import response
from numerical_guardrails import sha,fingerprint


def run(output):
    out=Path(output);out.mkdir(parents=True,exist_ok=False)
    region=load_region(CODE/'results/solar_a');reference=load_region(CODE/'results/solar_R')
    peak_path=CODE/'results/alpha_peak_checks.json'
    pinned=json.loads((ROOT/'2026_june_14_make_it_consistent/STEP11_ARTIFACTS.json').read_text())['sha256']
    if sha(peak_path)!=pinned[str(peak_path.relative_to(ROOT))]:raise ValueError('Peak checksum mismatch')
    peak=json.loads(peak_path.read_text());rows=[json.loads(l) for l in (CODE/'results/solar_a/evaluations.jsonl').read_text().splitlines()]
    cache={(r['A'],r['alpha']):r for r in rows};m=json.loads((CODE/'results/solar_a/manifest.json').read_text())
    aa=np.array(m['contract']['settings']['A']);al=np.array(m['contract']['settings']['alpha'])
    plt.rcParams.update({'font.size':10,'axes.labelsize':11,'legend.fontsize':8,'pdf.fonttype':42})
    fig,axes=plt.subplots(1,2,figsize=(10,3.6),layout='constrained')
    for ax,variant,title in zip(axes,['neutrino_only',VARIANTS[0]],['Neutrino-only score','Combined proxy score']):
        z=np.array([[score(cache[(float(A),float(a))],variant) for a in al] for A in aa])
        qmin=min(s['variants'][variant]['minimum_Q'] for s in region['slices']);delta=z-qmin
        im=ax.pcolormesh(al,aa,delta,shading='auto',cmap='viridis',vmin=0,vmax=12,rasterized=True)
        cs=ax.contour(al,aa,delta,levels=[2.30,6.18],colors=['white','orange'],linewidths=1)
        ax.clabel(cs,fmt={2.30:'2.30',6.18:'6.18'},fontsize=8)
        endpoints=[s['variants'][variant]['upper_endpoint'] for s in region['slices']]
        ax.plot(al,[np.nan if a is None else a for a in endpoints],'r--',lw=1.5,label=r'Slice endpoint ($\Delta Q=2.71$)')
        ax.set(xscale='log',xlabel=r'$\alpha$',ylabel='$A$',title=title,ylim=(0,.04));ax.legend(loc='upper left')
    fig.colorbar(im,ax=axes,extend='max',label=r'$Q-\min Q$ over sampled $\alpha$, refined $A$')
    fig.savefig(out/'score_map.pdf');fig.savefig(out/'score_map.png',dpi=160);plt.close(fig)
    fig,ax=plt.subplots(figsize=(5.2,3.6),layout='constrained')
    for variant,label,color,style in [(VARIANTS[0],'Primary combined','black','-'),('neutrino_only','Neutrinos only','C0','--'),('measurement_only','Measurement errors + proxy','C1',':')]:
        vals=[s['variants'][variant]['upper_endpoint'] for s in region['slices']]
        ax.plot(al,[np.nan if a is None else a for a in vals],color=color,ls=style,label=label)
        if any(v is None for v in vals):ax.scatter(al[[v is None for v in vals]],np.full(sum(v is None for v in vals),.0395),marker='^',s=10,c=color,label='Upper end grid-limited')
    ax.plot([s['alpha'] for s in reference['slices']],[s['variants'][VARIANTS[0]]['upper_endpoint'] for s in reference['slices']],color='C2',ls='-.',label=r'Combined, $r_\star=R_\odot$')
    ax.set(xscale='log',xlabel=r'$\alpha$',ylabel=r'Nominal slice endpoint $A_{\Delta Q=2.71}$',ylim=(0,.041));ax.grid(alpha=.25);ax.legend(loc='upper right');fig.savefig(out/'slice_endpoints.pdf');fig.savefig(out/'slice_endpoints.png',dpi=160);plt.close(fig)
    chosen=[0,16,32,40,48,56,64]
    table=[]
    for i in chosen:
        s=region['slices'][i];table.append({'alpha':s['alpha'],'combined':s['variants'][VARIANTS[0]]['upper_endpoint'],'neutrino':s['variants']['neutrino_only']['upper_endpoint']})
    lines=[r'\begin{tabular}{ccc}',r'\hline\hline',r'$\alpha$ & Neutrino only & Combined \\',r'\hline']
    for r in table:lines.append(f"{r['alpha']:.5f} & "+('grid-limited' if r['neutrino'] is None else f"{r['neutrino']:.6f}")+f" & {r['combined']:.6f}"+' '+chr(92)*2)
    lines += [r'\hline\hline',r'\end{tabular}'];(out/'slice_table.tex').write_text('\n'.join(lines)+'\n')
    density_points=[]
    for s in region['slices']:
        A=s['variants'][VARIANTS[0]]['upper_endpoint'];a=s['alpha'];r=cache[(A,a)];x=(r['rho_c']/150)**a
        density_points.append(dict(A=A,alpha=a,rho_c=r['rho_c'],d_core=3-3*A*x/(1+x),asymptote=3-3*A))
    p=peak['checks'][VARIANTS[0]];pr=next(r for r in peak['evaluations'] if r['A']==p['endpoint'] and r['alpha']==p['alpha']);x=(pr['rho_c']/150)**p['alpha']
    fiducials=[]
    for a in [.15,.30]:
        row=response(.0334,a);f=lambda A:score(response(A,a),VARIANTS[0]);opt=minimize_scalar(f,bounds=(0,.004),method='bounded',options={'xatol':2e-10})
        if not opt.success:raise RuntimeError('Fiducial minimum failed')
        fiducials.append({'response':row,'Q':score(row,VARIANTS[0]),'delta_Q':score(row,VARIANTS[0])-min(f(0),opt.fun)})
    summary={'slice_table':table,'peak':{**p,'rho_c':pr['rho_c'],'d_core':3-3*p['endpoint']*x/(1+x),'asymptote':3-3*p['endpoint']},'minimum_sampled_boundary_core_dimension':min(density_points,key=lambda r:r['d_core']),'fiducials':fiducials,'inputs_y':C['y'],'inputs_variance':C['variance']}
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    meta={'status':'complete','contract':fingerprint(ROOT,['paper1_SolarCore/solar_step12/figures.py','paper1_SolarCore/solar_step11/solver.py','paper1_SolarCore/solar_step11/inference.py'],{'solar_region_sha256':sha(CODE/'results/solar_a/region.json'),'reference_region_sha256':sha(CODE/'results/solar_R/region.json'),'peak_sha256':sha(peak_path)}),'sha256':{p.name:sha(p) for p in out.iterdir() if p.is_file()}}
    (out/'manifest.json').write_text(json.dumps(meta,indent=2)+'\n')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);a=p.parse_args();run(a.output)
