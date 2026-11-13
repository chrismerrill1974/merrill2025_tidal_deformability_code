"""
NS V7 paper figures. Run after validate.py passes.

  fig_mr_v7.pdf     : two panels — (a) M-R: GR baseline, RDT at the solar
                      ceiling (A=0.0087), CompOSE reference curve, D&H
                      benchmarks; (b) fixed-mass radius shift dR/R(M) for
                      boundary cases.
"""

import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from eos import CompOSE_SLy4, compose_reference
from tov import mr_curve

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "figures")
os.makedirs(FIG, exist_ok=True)

PC = np.geomspace(1e34, 1.3e36, 60)
CASES = [(0.120, 0.0087), (0.150, 0.0083), (0.334, 0.0061)]


def rising(curve):
    ok = ~np.isnan(curve[:, 1])
    M, R = curve[ok, 1], curve[ok, 2]
    i = M.argmax()
    return M[: i + 1], R[: i + 1]


eos = CompOSE_SLy4()
Mg, Rg = rising(mr_curve(eos, PC))

fig, ax = plt.subplots(1, 2, figsize=(11, 4.6))

ax[0].plot(Rg, Mg, "b-", lw=2, label="GR (SLy)")
Mr_, Rr_ = None, None
Rref, Mref = compose_reference()
ax[0].plot(Rref, Mref, "k:", lw=1.2, label="CompOSE reference")
Mx, Rx = rising(mr_curve(eos, PC, A=0.0087, alpha=0.12))
ax[0].plot(Rx, Mx, "r--", lw=2, label=r"RDT ceiling ($A=0.0087$)")
ax[0].plot([11.7, 9.99], [1.4, 2.05], "ks", ms=6, mfc="none",
           label="D&H 2001 benchmarks")
ax[0].axhspan(2.01, 2.15, color="g", alpha=0.12, label=r"heavy pulsars")
ax[0].set_xlabel("R (km)"); ax[0].set_ylabel(r"M ($M_\odot$)")
ax[0].set_xlim(9, 13); ax[0].set_ylim(0.4, 2.2)
ax[0].legend(fontsize=8, loc="lower left")
ax[0].set_title("(a) Mass–radius: GR vs maximum solar-allowed RDT")

for (alpha, A), st in zip(CASES, ("r--", "C1-.", "C4:")):
    Mx, Rx = rising(mr_curve(eos, PC, A=A, alpha=alpha))
    Mgrid = np.linspace(0.6, 1.9, 60)  # stop below M_max turning point
    dR = np.interp(Mgrid, Mx, Rx) / np.interp(Mgrid, Mg, Rg) - 1
    ax[1].plot(Mgrid, 1e3 * dR, st, lw=1.8,
               label=fr"$\alpha={alpha:g},\ A={A:.4f}$")
ax[1].axhline(0, color="k", lw=0.6)
ax[1].set_xlabel(r"M ($M_\odot$)")
ax[1].set_ylabel(r"$\Delta R/R$ at fixed mass ($\times10^{-3}$)")
ax[1].legend(fontsize=8)
ax[1].set_title("(b) Radius shift along the solar exclusion boundary")

plt.tight_layout()
out = os.path.join(FIG, "fig_mr_v7.pdf")
plt.savefig(out)
print("wrote", out)
