"""Plot verified finite samples; no interpolated confidence envelope."""
from pathlib import Path
import argparse,csv,json
from resource_control import require_controls

def main():
    require_controls()
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from scan import sha
    p=argparse.ArgumentParser();p.add_argument('--analysis',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    m=json.loads((a.analysis/'manifest.json').read_text());assert m['status']=='complete'
    for f,h in m['sha256'].items():assert sha(a.analysis/f)==h
    rows=list(csv.DictReader((a.analysis/'predictions.csv').open()))
    a.output.mkdir(parents=True,exist_ok=False)
    fig,axes=plt.subplots(2,2,figsize=(10,7),sharex=True)
    for col,emb in enumerate(['emb1','emb2']):
        for ref,color in [(5e5,'#4477aa'),(1e6,'#228833'),(2e6,'#cc6677')]:
            group=sorted((r for r in rows if r['embedding']==emb and float(r['rstar_cm'])==ref and float(r['fraction'])==1),key=lambda r:float(r['alpha']))
            for row,q in enumerate(['delta_mass','delta_R14']):
                axes[row,col].plot([float(r['alpha']) for r in group],[100*float(r[q]) for r in group],'.-',ms=3,lw=1,color=color,label=f'{ref/1e5:g} km reference')
        axes[0,col].set_title(f'Embedding {col+1}')
        axes[1,col].set_xlabel('Density exponent alpha')
        for ax in axes[:,col]:ax.set_xscale('log');ax.grid(alpha=.25);ax.axhline(0,color='gray',lw=.6)
    axes[0,0].set_ylabel('EOS-domain peak mass shift (%)')
    axes[1,0].set_ylabel('Radius shift at 1.4 solar masses (%)')
    axes[0,1].legend(fontsize=8)
    fig.suptitle('Full density law at sampled solar-slice endpoints\nRG(SLy4); finite-pressure surface; lines connect samples',fontsize=12)
    fig.tight_layout();fig.savefig(a.output/'parameter_dependence.png',dpi=160);fig.savefig(a.output/'parameter_dependence.pdf');plt.close(fig)
    (a.output/'manifest.json').write_text(json.dumps({'status':'complete','analysis_manifest_sha256':sha(a.analysis/'manifest.json'),'script_sha256':sha(Path(__file__)),'sha256':{n:sha(a.output/n) for n in ['parameter_dependence.png','parameter_dependence.pdf']}},indent=2)+'\n')
if __name__=='__main__':main()
