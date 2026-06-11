"""
Tidal V5 paper figure: (a) Lambda(M) for GR and the solar-allowed RDT
ceiling, with the GW170817 measured band; (b) fractional ceiling shift
dLambda/Lambda vs mass for the three boundary cases.
"""

import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from tidal import solve_tidal

HERE = os.path.dirname(os.path.abspath(__file__))
from eos import CompOSE_SLy4

FIG = os.path.join(HERE, "figures")
os.makedirs(FIG, exist_ok=True)

CASES = [("GR", 0.0, 0.2, "b-"),
         (r"RDT ceiling ($A_{95}$@$\alpha$=0.12)", 0.0087, 0.12, "r--"),
         (r"$A_{95}$@$\alpha$=0.15", 0.0083, 0.15, "C1-."),
         (r"$A_{95}$@$\alpha$=0.334", 0.0061, 0.334, "C4:")]
LOGPC = np.linspace(34.35, 35.95, 26)


def curve(A, alpha):
    eos = CompOSE_SLy4()
    out = []
    for lp in LOGPC:
        try:
            M, R, Cc, k2, L = solve_tidal(eos, 10.0**lp, A=A, alpha=alpha)
            out.append((M, L))
        except RuntimeError:
            pass
    arr = np.array(out)
    i = arr[:, 0].argmax()
    return arr[: i + 1, 0], arr[: i + 1, 1]


curves = {}
for name, A, alpha, st in CASES:
    curves[name] = curve(A, alpha)
    print(name, "done")

fig, ax = plt.subplots(1, 2, figsize=(11, 4.6))

for name, A, alpha, st in CASES[:2]:
    M, L = curves[name]
    ax[0].semilogy(M, L, st, lw=2, label=name)
ax[0].errorbar([1.186], [300], yerr=[[230], [420]], fmt="ks", ms=5,
               capsize=4, label=r"GW170817 $\tilde\Lambda$ (90%)")
ax[0].set_xlabel(r"M ($M_\odot$)")
ax[0].set_ylabel(r"$\Lambda$")
ax[0].set_xlim(1.0, 2.05)
ax[0].legend(fontsize=8)
ax[0].set_title("(a) Tidal deformability: GR vs maximum solar-allowed RDT")

Mg, Lg = curves["GR"]
Mgrid = np.linspace(1.0, 1.9, 50)
Lg_i = np.interp(Mgrid, Mg, Lg)
for name, A, alpha, st in CASES[1:]:
    M, L = curves[name]
    dL = np.interp(Mgrid, M, L) / Lg_i - 1
    ax[1].plot(Mgrid, 100 * dL, st, lw=1.8, label=name)
ax[1].axhline(0, color="k", lw=0.6)
ax[1].set_xlabel(r"M ($M_\odot$)")
ax[1].set_ylabel(r"$\Delta\Lambda/\Lambda$ (%)")
ax[1].legend(fontsize=8)
ax[1].set_title("(b) Ceiling shift along the solar exclusion boundary")

plt.tight_layout()
out = os.path.join(FIG, "fig_lambda_v5.pdf")
plt.savefig(out)
print("wrote", out)
