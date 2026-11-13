"""Manifested WD samples on the historical nominal solar boundary.
Requires --output; optionally --figures with a new directory. Existing archives
are never overwritten. log-radius shifts compare at modified mass, while the
historical horizontal coordinate is baseline mass (see step 07).
"""

import argparse
import math
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from model import Params, mr_pair, MSUN_CGS, THETA_FLOOR
from numerical_guardrails import Artifact, fingerprint, positive
import validate

HERE = os.path.dirname(os.path.abspath(__file__))


# Paper 1 V3 legacy nominal boundary samples (no coverage claim),
# from solarcore_v3/analyze.py output, 2026-06-10.
A95 = [
    (0.020, 0.0042), (0.030, 0.0047), (0.045, 0.0055), (0.067, 0.0068),
    (0.100, 0.0083), (0.150, 0.0083), (0.224, 0.0066), (0.334, 0.0061),
    (0.500, 0.0067),
]

# Mathematical polytrope grid; its upper portion exceeds the
# nonrelativistic applicability estimate (step 07).
RHO_C = np.geomspace(1e5, 3e6, 12)
FAST = dict(rtol=1e-9, atol=1e-11, max_step=0.05)


FIELDS=['alpha','A','rho_c','M_std','R_std_km','dlnR','dlnM','dlnR_fixedM']


def check_row(r):
    positive(rho_c=r['rho_c'],mass=r['M_std'],radius=r['R_std_km'])
    if not math.isclose(r['dlnR_fixedM'],r['dlnR']+r['dlnM']/3,rel_tol=1e-12,abs_tol=1e-12):
        raise ValueError('Inconsistent fixed-mass shift')


def artifact(output, points, densities):
    points=[list(map(float,p)) for p in points];densities=list(map(float,densities))
    if (not points or not densities or len(points)!=len({tuple(p) for p in points})
            or any(len(p)!=2 or not .02<=p[0]<=.5 or not 0<=p[1]<1 for p in points)):
        raise ValueError('Invalid white-dwarf parameter sample')
    positive(**{str(i):r for i,r in enumerate(densities)})
    settings={'points':points,'rho_c':densities,'solver':FAST,'theta_floor':THETA_FLOOR,
              'source':'legacy nominal solar-boundary samples; not new inference'}
    contract=fingerprint(HERE,['model.py','predict.py','validate.py'],settings)
    return Artifact(output,contract,FIELDS,['alpha','A','rho_c'],
                    [(a,A,r) for a,A in points for r in densities],check_row)


def main(output, points=A95, densities=RHO_C):
    data=artifact(output,points,densities);data.load()
    point={'phase':'validation'}
    try:
        if validate.run_all() is not True: raise RuntimeError('Validation did not pass')
        done={(r['alpha'],r['A'],r['rho_c']) for r in data.rows}
        for alpha,A in points:
            for rc in densities:
                if (alpha,A,rc) in done: continue
                point={'alpha':float(alpha),'A':float(A),'rho_c':float(rc)}
                std,mod=mr_pair(Params(A,alpha),rc,**FAST)
                positive(M_std=std.M,R_std=std.R,M_mod=mod.M,R_mod=mod.R)
                dr=float(np.log(mod.R/std.R));dm=float(np.log(mod.M/std.M))
                data.append(dict(alpha=alpha,A=A,rho_c=rc,M_std=std.M/MSUN_CGS,
                                 R_std_km=std.R/1e5,dlnR=dr,dlnM=dm,dlnR_fixedM=dr+dm/3))
        data.save('complete')
    except Exception as exc:
        data.fail(point,exc)
        raise
    print(f'Complete: {len(data.rows)} white-dwarf pairs; {output}')
    return data.rows


def plot_complete(output, figure_dir):
    # Only the canonical full sample is plotted by this helper.
    data=artifact(output,A95,RHO_C);rows=data.load(complete=True)
    if os.path.exists(figure_dir): raise ValueError('Use a new figure directory')
    FIGDIR=figure_dir;os.makedirs(FIGDIR)
    # Selected samples over alpha at each rho_c
    print("\n=== Paper-numbers block (legacy nominal-boundary samples) ===")
    print(f"{'rho_c':>9} {'M_std':>7} {'max|dlnR|_M|':>13} {'at alpha':>9}")
    env = []
    for rc in sorted({r['rho_c'] for r in rows}):
        sub = [r for r in rows if r["rho_c"] == rc]
        worst = max(sub, key=lambda r: abs(r["dlnR_fixedM"]))
        env.append((rc, worst["M_std"], worst["dlnR_fixedM"], worst["alpha"]))
        print(f"{rc:9.1e} {worst['M_std']:7.3f} {worst['dlnR_fixedM']:13.4e} "
              f"{worst['alpha']:9.3f}")
    gmax = max(env, key=lambda e: abs(e[2]))
    print(f"\nSampled maximum: |dR/R|_M <= {abs(gmax[2]):.3e} "
          f"(at rho_c={gmax[0]:.1e}, M={gmax[1]:.2f} Msun, alpha={gmax[3]})")

    # Figure: dlnR at fixed mass vs mass, one curve per alpha on the boundary
    plt.figure(figsize=(7, 4.8))
    cmap = plt.cm.viridis(np.linspace(0.15, 0.9, len(A95)))
    for (alpha, A), c in zip(A95, cmap):
        sub = sorted((r for r in rows if r["alpha"] == alpha),
                     key=lambda r: r["M_std"])
        plt.plot([r["M_std"] for r in sub],
                 [1e3 * r["dlnR_fixedM"] for r in sub],
                 color=c, lw=1.5, label=fr"$\alpha={alpha:g},\ A=A_{{95}}$")
    plt.axhline(0, color="k", lw=0.6)
    plt.xscale("log")
    plt.xlabel(r"$M_0$ ($M_\odot$), baseline coordinate")
    plt.ylabel(r"$\Delta\ln R$ at modified mass ($\times 10^{-3}$)")
    plt.title("Legacy nominal-boundary polytrope samples")
    plt.legend(fontsize=7, ncol=2)
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(FIGDIR, "fig_dR_bands.pdf"))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description='Generate manifested legacy nominal WD samples')
    parser.add_argument('--output',required=True)
    parser.add_argument('--figures',help='Optional new figure directory')
    args=parser.parse_args()
    main(args.output)
    if args.figures: plot_complete(args.output,args.figures)
