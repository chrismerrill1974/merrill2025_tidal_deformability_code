"""
Tidal V6 paper figure: tidal deformability on the CONTRACTED operator
backgrounds (Embedding 1), the V6 replacement for the V5 figures.py
(archive_v5/). Uses the canonical divergence-operator mechanism in place of
the V5 F-multiplier.

(a) Lambda(M) for GR and the solar-allowed RDT ceiling, with the GW170817
    measured band; the RDT curve lies BELOW GR (a reduction), the mirror of
    V5's figure where it lay above.
(b) Fractional ceiling shift dLambda/Lambda vs mass for the three boundary
    cases; now NEGATIVE (a reduction), bounded by the ceiling.

Reuses solve_op_tidal from operator_tidal.py (Emb1, R_op = a_op/a_N on both
dP/dr and nu'). A=0 reduces to GR (gate: Lambda(1.4) ~ 298).
"""

import os
import sys

import numpy as np
from scipy.optimize import brentq
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from operator_tidal import solve_op_tidal, CompOSE_SLy4  # noqa: E402

FIG = os.path.join(HERE, "figures")
os.makedirs(FIG, exist_ok=True)

# alpha=0.122 ceiling matches Paper III V8 and operator_tidal.py
CASES = [("GR", 0.0, 0.122, "b-"),
         (r"RDT ceiling ($A_{95}$@$\alpha$=0.122)", 0.0087, 0.122, "r--"),
         (r"$A_{95}$@$\alpha$=0.15", 0.0083, 0.15, "C1-."),
         (r"$A_{95}$@$\alpha$=0.334", 0.0061, 0.334, "C4:")]
LOGPC = np.linspace(34.4, 35.9, 24)


def curve(A, alpha):
    eos = CompOSE_SLy4()
    out = []
    for lp in LOGPC:
        try:
            M, R, Cc, k2, L = solve_op_tidal(eos, 10.0**lp, A, alpha)
            out.append((M, L))
        except RuntimeError:
            pass
    arr = np.array(out)
    i = arr[:, 0].argmax()
    return arr[: i + 1, 0], arr[: i + 1, 1]


def lambda_at_mass(eos, Mt, A, alpha):
    """Lambda at a fixed target mass, by root-finding on central pressure.
    Clean fractional shifts (mirrors operator_tidal.at_mass) — avoids the
    interpolation noise of differencing two independently-sampled curves."""
    f = lambda lp: solve_op_tidal(eos, 10.0**lp, A, alpha)[0] - Mt
    lp = brentq(f, 34.3, 35.85, xtol=1e-6)
    return solve_op_tidal(eos, 10.0**lp, A, alpha)[4]


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

eos = CompOSE_SLy4()
Mgrid = np.linspace(1.0, 1.8, 9)
Lg_i = np.array([lambda_at_mass(eos, Mt, 0.0, 0.122) for Mt in Mgrid])
for name, A, alpha, st in CASES[1:]:
    Lr = np.array([lambda_at_mass(eos, Mt, A, alpha) for Mt in Mgrid])
    dL = Lr / Lg_i - 1
    ax[1].plot(Mgrid, 100 * dL, st, lw=1.8, label=name)
ax[1].axhline(0, color="k", lw=0.6)
ax[1].set_xlabel(r"M ($M_\odot$)")
ax[1].set_ylabel(r"$\Delta\Lambda/\Lambda$ (%)")
ax[1].legend(fontsize=8)
ax[1].set_title("(b) Ceiling shift along the solar exclusion boundary")

plt.tight_layout()
out = os.path.join(FIG, "fig_lambda_v6.pdf")
plt.savefig(out)
print("wrote", out)
