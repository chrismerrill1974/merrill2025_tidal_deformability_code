"""
Tidal V5 predictions: maximum solar-allowed RDT tidal signatures.

For GR and points on Paper 1's combined 95% CL exclusion boundary
A_95(alpha), compute Lambda(M) at fixed masses, decompose the shift into
Love-number and R^5 contributions, and evaluate the GW170817 binary
combined deformability Lambda-tilde (m1 = 1.46, m2 = 1.27 Msun, low-spin).

Writes results_tidal.csv and prints the paper-numbers block. ~2 min.
"""

import csv
import os
import sys

import numpy as np
from scipy.optimize import brentq

from tidal import solve_tidal, lambda_tilde

HERE = os.path.dirname(os.path.abspath(__file__))
from eos import CompOSE_SLy4

# Paper 1 V3 combined boundary (same representative points as Papers 2-3)
CASES = [("GR", 0.0, 0.2),
         ("A95@0.12", 0.0087, 0.12),
         ("A95@0.15", 0.0083, 0.15),
         ("A95@0.334", 0.0061, 0.334)]
MASSES = [1.0, 1.2, 1.4, 1.6]
GW_M1, GW_M2 = 1.46, 1.27


def star_at_mass(eos, M_target, A, alpha):
    f = lambda lp: solve_tidal(eos, 10.0**lp, A=A, alpha=alpha)[0] - M_target
    logPc = brentq(f, 34.2, 36.0, xtol=1e-7)
    return solve_tidal(eos, 10.0**logPc, A=A, alpha=alpha)


def main():
    eos = CompOSE_SLy4()
    results = {}
    for name, A, alpha in CASES:
        for M in MASSES + [GW_M1, GW_M2]:
            results[(name, M)] = star_at_mass(eos, M, A, alpha)
        print(f"{name} done")

    rows = []
    print("\n=== Paper-numbers block ===")
    print(f"{'case':>10} {'M':>4} {'R km':>7} {'k2':>7} {'Lambda':>8} "
          f"{'dR/R %':>7} {'dk2/k2 %':>8} {'dLam/Lam %':>10} {'5dR/R %':>8}")
    for M in MASSES:
        _, Rg, _, k2g, Lg = results[("GR", M)]
        for name, A, alpha in CASES:
            _, R, _, k2, L = results[(name, M)]
            dR, dk2, dL = R / Rg - 1, k2 / k2g - 1, L / Lg - 1
            print(f"{name:>10} {M:4.1f} {R:7.3f} {k2:7.4f} {L:8.1f} "
                  f"{100*dR:+7.2f} {100*dk2:+8.2f} {100*dL:+10.2f} "
                  f"{100*5*dR:+8.2f}")
            rows.append(dict(case=name, A=A, alpha=alpha, M=M, R=R,
                             k2=k2, Lam=L, dR_frac=dR, dk2_frac=dk2,
                             dLam_frac=dL))

    print("\nGW170817 binary (m1=1.46, m2=1.27 Msun, low-spin):")
    for name, A, alpha in CASES:
        _, _, _, _, L1 = results[(name, GW_M1)]
        _, _, _, _, L2 = results[(name, GW_M2)]
        Lt = lambda_tilde(GW_M1, GW_M2, L1, L2)
        flag = "compatible" if Lt < 720 else "EXCEEDS"
        print(f"  {name:>10}: Lambda1 = {L1:6.1f}, Lambda2 = {L2:6.1f}, "
              f"Lambda-tilde = {Lt:6.1f}  ({flag} GW170817 bound 720)")

    out = os.path.join(HERE, "results_tidal.csv")
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"\nWrote {out}")


if __name__ == "__main__":
    main()
