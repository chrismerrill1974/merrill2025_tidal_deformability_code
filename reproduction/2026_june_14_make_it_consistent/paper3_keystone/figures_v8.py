"""
Paper III V8 figure: operator-form CONTRACTION (replaces V7's expansion fig).
(a) M-R: GR baseline, CompOSE reference, RDT at the solar ceiling (inside GR).
(b) fixed-mass radius contraction along the A95(alpha) boundary (Embedding 1).
Writes paper3_NeutronStars/V8_Paper/fig_mr_v8.pdf.
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(
    HERE, "..", "..", "paper3_NeutronStars", "ns_v7")))
from eos import CompOSE_SLy4, compose_reference
from tov_operator import mr_curve_op

OUT = os.path.normpath(os.path.join(
    HERE, "..", "..", "paper3_NeutronStars", "V8_Paper", "fig_mr_v8.pdf"))

eos = CompOSE_SLy4()
PC = np.geomspace(1e34, 1.3e36, 60)
BOUNDARY = [(0.122, 0.0087), (0.15, 0.0083), (0.334, 0.0061)]


def rising(curve):
    ok = ~np.isnan(curve[:, 1])
    M, R = curve[ok, 1], curve[ok, 2]
    i = int(M.argmax())
    return M[:i + 1], R[:i + 1]


gr = mr_curve_op(eos, PC, 0.0, 0.122)
rdt = mr_curve_op(eos, PC, 0.0087, 0.122, embedding="emb1")
Mg, Rg = rising(gr)
Mr, Rr = rising(rdt)
Rref, Mref = compose_reference()

fig, ax = plt.subplots(1, 2, figsize=(11, 4.6))

ax[0].plot(Rg, Mg, "b-", lw=2, label="GR (SLy)")
ax[0].plot(Rref, Mref, "k:", lw=1.2, label="CompOSE reference")
ax[0].plot(Rr, Mr, "r--", lw=2, label=r"RDT operator form ($A=0.0087$)")
ax[0].plot([11.7, 9.99], [1.4, 2.05], "ks", ms=6, mfc="none",
           label="D\\&H 2001 benchmarks")
ax[0].axhspan(2.01, 2.15, color="g", alpha=0.12, label="heavy pulsars")
ax[0].set_xlabel("R (km)")
ax[0].set_ylabel(r"M ($M_\odot$)")
ax[0].set_xlim(9, 13)
ax[0].set_ylim(0.4, 2.2)
ax[0].legend(fontsize=8, loc="lower left")
ax[0].set_title("(a) Mass--radius: GR vs solar-allowed RDT (operator form)")

Mgrid = np.linspace(1.0, 1.9, 60)
for (alpha, A), st in zip(BOUNDARY, ("r--", "C1-.", "C4:")):
    c = mr_curve_op(eos, PC, A, alpha, embedding="emb1")
    Mc, Rc = rising(c)
    dR = np.interp(Mgrid, Mc, Rc) / np.interp(Mgrid, Mg, Rg) - 1.0
    ax[1].plot(Mgrid, 1e2 * dR, st, lw=1.8,
               label=fr"$\alpha={alpha:g},\ A={A:.4f}$")
ax[1].axhline(0, color="k", lw=0.6)
ax[1].set_xlabel(r"M ($M_\odot$)")
ax[1].set_ylabel(r"$\Delta R/R$ at fixed mass (\%)")
ax[1].legend(fontsize=8)
ax[1].set_title("(b) Radius contraction along the solar exclusion boundary")

plt.tight_layout()
plt.savefig(OUT)
print("wrote", OUT)
