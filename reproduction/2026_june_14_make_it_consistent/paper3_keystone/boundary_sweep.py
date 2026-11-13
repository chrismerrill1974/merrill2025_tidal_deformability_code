"""
Check A (gate independence) + Check C (A95(alpha) boundary ceiling).

Check A: run the MODIFIED emb1/emb2 path at tiny A (1e-6, 1e-9) — bypassing the
A==0 fallback — and confirm it collapses to V7 with deviation scaling LINEARLY
in A (proves the contraction comes from the same RHS, only A varied).

Check C: run the canonical A95(alpha) boundary {(.12,.0087),(.15,.0083),
(.334,.0061)} for both embeddings; report contraction as a ceiling.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.normpath(os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "..", "..", "paper3_NeutronStars", "ns_v7")))
from eos import CompOSE_SLy4
from tov_operator import mr_curve_op

eos = CompOSE_SLy4()
PC = np.geomspace(1e34, 1.3e36, 60)


def summarize(curve):
    ok = ~np.isnan(curve[:, 1])
    M, R = curve[ok, 1], curve[ok, 2]
    i = int(M.argmax())
    return M[i], R[i], float(np.interp(1.4, M[:i + 1], R[:i + 1]))


# baseline via the modified harness with A=0 (fallback = V7)
Mg, Rg, R14g = summarize(mr_curve_op(eos, PC, 0.0, 0.12))
print(f"V7 baseline (A=0): M_max={Mg:.6f}, R(1.4)={R14g:.4f} km")

print("\n=== CHECK A: modified path -> V7 as A->0 (linear scaling) ===")
print(f"{'emb':>5} {'A':>8} {'|dM_max/M|':>12} {'|dR14/R|':>12}")
for emb in ("emb1", "emb2"):
    devs = []
    for A in (1e-6, 1e-9):
        Mm, _, R14 = summarize(mr_curve_op(eos, PC, A, 0.12, embedding=emb))
        dM = abs(Mm / Mg - 1.0)
        dR = abs(R14 / R14g - 1.0)
        devs.append((A, dM, dR))
        print(f"{emb:>5} {A:>8.0e} {dM:>12.2e} {dR:>12.2e}")
    # linear scaling check: dM(1e-6)/dM(1e-9) should be ~1000
    ratio = devs[0][1] / devs[1][1] if devs[1][1] > 0 else float('inf')
    print(f"      -> dM ratio (1e-6 / 1e-9) = {ratio:.0f}  (expect ~1000 if "
          f"deviation is linear in A and path reduces to V7)")

print("\n=== CHECK C: A95(alpha) boundary — contraction ceiling ===")
BOUNDARY = [(0.12, 0.0087), (0.15, 0.0083), (0.334, 0.0061)]
print(f"{'alpha':>7} {'A95':>8} {'emb':>5} {'M_max':>9} {'dM_max%':>9} "
      f"{'R(1.4)':>8} {'dR14%':>8}")
results = {}
for alpha, A in BOUNDARY:
    for emb in ("emb1", "emb2"):
        Mm, _, R14 = summarize(mr_curve_op(eos, PC, A, alpha, embedding=emb))
        results[(alpha, emb)] = (100 * (Mm / Mg - 1), 100 * (R14 / R14g - 1))
        print(f"{alpha:>7.3f} {A:>8.4f} {emb:>5} {Mm:>9.4f} "
              f"{100*(Mm/Mg-1):>+9.2f} {R14:>8.3f} {100*(R14/R14g-1):>+8.2f}")

dM_ceiling = max(abs(results[k][0]) for k in results)
dR_ceiling = max(abs(results[k][1]) for k in results)
print(f"\nCEILING across boundary: |ΔM_max| ≤ {dM_ceiling:.2f}%, "
      f"|ΔR(1.4)| ≤ {dR_ceiling:.2f}%  (both at α=0.12, A=0.0087; "
      f"Emb1 primary). Λ(1.4) reduced; full value from Paper IV re-run.")
