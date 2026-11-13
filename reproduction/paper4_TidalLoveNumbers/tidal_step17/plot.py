"""Plot sampled conditional tides and exact radius/k2 decomposition."""
from pathlib import Path
import json
from resource_control import require_controls
HERE=Path(__file__).resolve().parent

def main():
 require_controls()
 import matplotlib
 matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 s=json.loads((HERE/'results/summary.json').read_text());assert json.loads((HERE/'verification.json').read_text())['status']=='pass'
 fig,ax=plt.subplots(1,2,figsize=(10,4.5))
 for emb,color in [('emb1','#4477aa'),('emb2','#cc6677')]:
  group=sorted((r for r in s['pairs'] if r['label']=='solar_peak' and r['embedding']==emb),key=lambda r:r['mass'])
  ax[0].plot([v['mass'] for v in group],[100*v['delta_Lambda'] for v in group],'o-',color=color,label=emb.replace('emb','Embedding '))
 group=next(v for v in s['pairs'] if v['label']=='solar_peak' and v['embedding']=='emb1' and v['mass']==1.4)
 labels=['Full tidal','Radius only','Love number','Cross term']
 values=[group[q]*100 for q in ['delta_Lambda','radius_only_exact','delta_k2','multiplicative_cross_term']]
 ax[1].bar(labels,values,color=['#4477aa','#228833','#cc6677','#aaaaaa'])
 ax[1].set_title('Embedding 1 at 1.4 solar masses');ax[1].tick_params(axis='x',rotation=15)
 ax[0].set_xlabel('Mass (solar masses)');ax[0].legend();ax[0].set_ylabel('Fractional tidal shift (%)')
 ax[1].set_ylabel('Contribution (%)')
 for a in ax:a.axhline(0,color='gray',lw=.7);a.grid(axis='y',alpha=.2)
 fig.suptitle('Sampled solar-amplitude-peak case; CompOSE RG(SLy4)\nConditional GR-form tides, prescribed finite-cutoff matching',fontsize=11)
 fig.tight_layout();out=HERE/'figures';out.mkdir(exist_ok=False);fig.savefig(out/'tidal_response.png',dpi=160);fig.savefig(out/'tidal_response.pdf');plt.close(fig)
if __name__=='__main__':main()
