"""Paper II figure at the actual modified comparison mass, with exact fractional shifts."""
from pathlib import Path
import argparse,json,sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
from numerical_guardrails import sha,fingerprint
HERE=Path(__file__).resolve().parent

def run(output):
    out=Path(output);out.mkdir(parents=True,exist_ok=False)
    meta=json.loads((HERE/'results/manifest.json').read_text())
    if meta['status']!='complete':raise ValueError('Incomplete WD run')
    for f,h in meta['sha256'].items():
        if sha(HERE/'results'/f)!=h:raise ValueError('WD checksum mismatch')
    for f,h in meta['contract']['source_sha256'].items():
        if sha(ROOT/f)!=h:raise ValueError('WD source mismatch')
    rows=[json.loads(l) for l in (HERE/'results/evaluations.jsonl').read_text().splitlines()]
    boundary=[r for r in rows if r['fraction']==1]
    alphas=sorted(set(r['alpha'] for r in boundary if r['origin']=='sampled_slice'))
    selected=alphas[::8];peak=next(r['alpha'] for r in boundary if r['origin']=='local_solar_peak');selected=sorted([*selected,peak])
    plt.rcParams.update({'font.size':10,'legend.fontsize':8,'pdf.fonttype':42})
    fig,ax=plt.subplots(figsize=(7,4),layout='constrained')
    for alpha,color in zip(selected,plt.cm.viridis(np.linspace(0,.95,len(selected)))):
        r=sorted([r for r in boundary if r['alpha']==alpha],key=lambda r:r['M_mod'])
        ax.plot([v['M_mod'] for v in r],[100*v['fractional_fixedM'] for v in r],color=color,label=fr'$\alpha={alpha:.4f}$')
        extra=[v for v in r if v['rho_c']>2e6]
        ax.scatter([v['M_mod'] for v in extra],[100*v['fractional_fixedM'] for v in extra],color=color,marker='s',s=13)
    ax.set(xlabel=r'Comparison mass $M_A(\rho_c)/M_\odot$',ylabel=r'Radius shift at $M_A$, $\Delta R/R$ (%)',title='Solar-a nominal boundary: selected polytropic samples')
    ax.grid(alpha=.25);ax.legend(ncol=2,loc='center left',bbox_to_anchor=(1,0.5))
    fig.savefig(out/'radius_shifts.pdf');fig.savefig(out/'radius_shifts.png',dpi=180);plt.close(fig)
    lines=[r'\begin{tabular}{ccc}',r'\hline\hline',r'$\alpha$ & $A_{\Delta Q=2.71}$ & $\Delta R/R$ (\%) '+chr(92)*2,r'\hline']
    for a in alphas[::16]:
        r=max((r for r in boundary if r['alpha']==a),key=lambda r:r['rho_c'])
        lines.append(f"{a:.5f} & {r['A']:.6f} & {100*r['fractional_fixedM']:.5f}"+' '+chr(92)*2)
    lines += [r'\hline\hline',r'\end{tabular}'];(out/'boundary_table.tex').write_text('\n'.join(lines)+'\n')
    (out/'manifest.json').write_text(json.dumps({'status':'complete','contract':fingerprint(ROOT,['paper2_WhiteDwarfs/wd_step12/figures.py'],{'input_manifest_sha256':sha(HERE/'results/manifest.json')}),'sha256':{p.name:sha(p) for p in out.iterdir()}},indent=2)+'\n')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);a=p.parse_args();run(a.output)
