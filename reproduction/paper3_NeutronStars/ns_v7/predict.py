"""
NS V7 predictions: RDT neutron-star shifts along the Paper 1 (V3) combined
95% CL boundary A_95(alpha) — the largest signatures allowed by solar data.

At NS densities the opening law is fully saturated for all allowed alpha
(x = (rho/rho0)^alpha >> 1), so F = (d_eff - 1)/2 ~ 1 - (3/2)A = const:
an interior G_eff rescaling. EOS-independence and alpha-saturation are
therefore BUILT IN, not findings. The envelope is set by max A_95 = 0.0087
(alpha ~ 0.10-0.15) => F >= 0.98695.

Writes results_predictions.csv and prints the paper-numbers block.
"""

import csv
import os

import numpy as np

from eos import CompOSE_SLy4
from tov import mr_curve

HERE = os.path.dirname(os.path.abspath(__file__))

# Paper 1 V3 combined A_95(alpha) — representative points incl. the global max
BOUNDARY = [(0.100, 0.0083), (0.150, 0.0083), (0.334, 0.0061)]
A_MAX = 0.0087  # global bound; quoted envelope uses this at alpha = 0.12

PC = np.geomspace(1e34, 1.3e36, 60)


def stats(curve):
    ok = ~np.isnan(curve[:, 1])
    M, R = curve[ok, 1], curve[ok, 2]
    i = M.argmax()
    return M[i], R[i], np.interp(1.4, M[:i + 1], R[:i + 1]), M[:i + 1], R[:i + 1]


def main():
    eos = CompOSE_SLy4()
    g = mr_curve(eos, PC)
    Mm0, Rm0, R140, Mg, Rg = stats(g)
    print(f"GR baseline: M_max = {Mm0:.4f} Msun @ {Rm0:.3f} km, "
          f"R(1.4) = {R140:.3f} km")

    rows = []
    print(f"\n{'alpha':>6} {'A':>7} {'F_floor':>8} {'M_max':>7} {'dM/M %':>7} "
          f"{'R(1.4)':>7} {'dR/R %':>7} {'dLam/Lam %':>10}")
    cases = BOUNDARY + [(0.120, A_MAX)]
    for alpha, A in cases:
        c = mr_curve(eos, PC, A=A, alpha=alpha)
        Mm, Rm, R14, _, _ = stats(c)
        dM = Mm / Mm0 - 1
        dR = R14 / R140 - 1
        F = 1 - 1.5 * A
        print(f"{alpha:6.3f} {A:7.4f} {F:8.5f} {Mm:7.4f} {100*dM:+7.2f} "
              f"{R14:7.3f} {100*dR:+7.2f} {100*5*dR:+10.1f}")
        rows.append(dict(alpha=alpha, A=A, F_floor=F, M_max=Mm,
                         dMmax_frac=dM, R14=R14, dR14_frac=dR))

    out = os.path.join(HERE, "results_predictions.csv")
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"\n(GR: M_max={Mm0:.4f}, R14={R140:.3f}; dLam/Lam ~ 5 dR/R)")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
