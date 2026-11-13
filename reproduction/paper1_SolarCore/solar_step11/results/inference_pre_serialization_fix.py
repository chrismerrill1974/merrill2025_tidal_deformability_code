"""Versioned step 09 scores and all-component nominal slice sets."""
import json
from pathlib import Path
import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.optimize import minimize_scalar, brentq
ROOT=Path(__file__).resolve().parents[2]
CONTRACT=ROOT/'2026_june_14_make_it_consistent/STEP09_INFERENCE_CONTRACT.json'
INPUT=ROOT/'2026_june_14_make_it_consistent/STEP09_INPUT_AUDIT.json'
C=json.loads(CONTRACT.read_text());ROWS=json.loads(INPUT.read_text())['source_rows']
BRANCHES=['pp','pep','be7','b8']
VARIANTS=['source_precision_diagonal_proxy_v1','legacy_diagonal_proxy','measurement_only','helioseismic_scale_0.001','helioseismic_scale_0.004','omit_pep','all_lower_error_sides','all_upper_error_sides','neutrino_only']

def score(row,variant):
    if variant not in VARIANTS:raise ValueError('Unknown statistic')
    y=np.array([C['y'][b] for b in BRANCHES]);v=np.array([C['variance'][b] for b in BRANCHES])
    scale=.002;mask=np.ones(4,dtype=bool)
    if variant=='legacy_diagonal_proxy':
        y=np.array([r['legacy_ratio'] for r in ROWS]);v=np.array([r['legacy_sigma_total']**2 for r in ROWS])
    if variant=='measurement_only':v=np.array([r['sigma_ratio_rms']**2 for r in ROWS])
    if variant in ('all_lower_error_sides','all_upper_error_sides'):
        side='lower' if variant=='all_lower_error_sides' else 'upper'
        v=np.array([r['sigma_ratio_'+side]**2+r['theory_fraction']**2 for r in ROWS])
    if variant=='omit_pep':mask[1]=False
    if variant.startswith('helioseismic_scale_'):scale=float(variant.rsplit('_',1)[1])
    r=np.array([row['R_'+b] for b in BRANCHES])
    nu=float(np.sum(((r-y)**2/v)[mask]));h=0 if variant=='neutrino_only' else (row['h']/scale)**2
    return nu+float(h)


def intervals(x,values,threshold=2.71,evaluate=None):
    """Locate every component on the supplied mesh; refine each sign change.
    Does not certify absence of narrower components than the checked mesh.
    evaluate, when supplied, is a fresh ODE-based score, not the interpolant.
    """
    x=np.asarray(x);values=np.asarray(values)
    if len(x)<3 or np.any(np.diff(x)<=0) or not np.isfinite(values).all():raise ValueError('Invalid slice')
    interp=PchipInterpolator(x,values);f=interp if evaluate is None else evaluate
    candidates=[(float(values[0]),float(x[0])),(float(values[-1]),float(x[-1]))]
    brackets=[]
    if float(f(x[0]+.01*(x[1]-x[0])))<values[0]:brackets.append((x[0],x[1]))
    if float(f(x[-1]-.01*(x[-1]-x[-2])))<values[-1]:brackets.append((x[-2],x[-1]))
    for i in range(1,len(x)-1):
        if values[i]<=values[i-1] and values[i]<=values[i+1]:brackets.append((x[i-1],x[i+1]))
    for bounds in brackets:
            fit=minimize_scalar(f,bounds=bounds,method='bounded',options={'xatol':2e-10})
            if not fit.success:raise RuntimeError('Minimum refinement failed')
            candidates.append((float(fit.fun),float(fit.x)))
    qmin,amin=min(candidates);level=qmin+threshold
    # Include refined minima so a component around a minimum between knots survives.
    nodes=sorted(set(map(float,x))|{a for q,a in candidates})
    vals=[float(f(a))-level for a in nodes]
    roots=[]
    for a,b,fa,fb in zip(nodes[:-1],nodes[1:],vals[:-1],vals[1:]):
        if fa==0:roots.append(a)
        if fa*fb<0:roots.append(float(brentq(lambda z:float(f(z))-level,a,b,xtol=2e-10)))
    if vals[-1]==0:roots.append(nodes[-1])
    cuts=sorted(set([float(x[0]),float(x[-1]),*roots]));accepted=[]
    for lo,hi in zip(cuts[:-1],cuts[1:]):
        if float(f((lo+hi)/2))<=level:
            accepted.append({'lower':lo,'upper':hi,'upper_grid_limited':hi==x[-1]})
    # Isolated exactly tangent threshold points are kept explicitly.
    for r in roots:
        if not any(a['lower']<=r<=a['upper'] for a in accepted):accepted.append({'lower':r,'upper':r,'upper_grid_limited':r==x[-1]})
    accepted.sort(key=lambda a:a['lower'])
    return {'minimum_A':amin,'minimum_Q':qmin,'intervals':accepted,
            'upper_endpoint':None if not accepted or accepted[-1]['upper_grid_limited'] else accepted[-1]['upper'],
            'maximum_root_score_residual':max([abs(float(f(r))-level) for r in roots] or [0.])}
