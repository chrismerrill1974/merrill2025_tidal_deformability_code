"""
Check A done right — gate independence at FIXED central pressure (no grid/argmax
artifacts), baselining the modified path against ITSELF at delta->0.

For each embedding, single star at P_c = 5e35:
  (1) path reproduction: modified path at A=1e-15 vs the V7 2-ODE solver
      -> how well the modified code path reproduces V7 at delta=0;
  (2) linear scaling: modified path at A=1e-6, 1e-9 vs the modified A=1e-15
      baseline -> deviation must scale ~linearly in A (ratio ~1000), proving
      the A-dependence is clean within the SAME path;
  (3) contraction: modified path at A=0.0087 vs the modified A=1e-15 baseline
      -> the contraction must appear relative to the path's OWN delta=0 limit,
      i.e. "identical code path, only A varied."
A=1e-15 != 0.0 so it bypasses the A==0 V7 fallback and exercises the 4-ODE RHS.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.normpath(os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "..", "..", "paper3_NeutronStars", "ns_v7")))
from eos import CompOSE_SLy4
from tov import solve_star as solve_star_v7
from tov_operator import solve_star_op

eos = CompOSE_SLy4()
PC = 5.0e35
ALPHA = 0.12

M_v7, R_v7 = solve_star_v7(eos, PC, A=0.0)
print(f"V7 (2-ODE) at P_c=5e35: M={M_v7:.8f}, R={R_v7:.6f}")

for emb in ("emb1", "emb2"):
    print(f"\n--- {emb} ---")
    Mb, Rb = solve_star_op(eos, PC, 1e-15, ALPHA, embedding=emb)  # forced modified path, delta~0
    repM, repR = abs(Mb / M_v7 - 1), abs(Rb / R_v7 - 1)
    print(f"(1) path reproduction (A=1e-15 vs V7): |dM|={repM:.2e}, |dR|={repR:.2e}"
          f"  -> {'OK (integrator level)' if repM < 1e-5 else 'INVESTIGATE'}")

    devs = {}
    for A in (1e-6, 1e-9):
        M, R = solve_star_op(eos, PC, A, ALPHA, embedding=emb)
        devs[A] = abs(M / Mb - 1)
        print(f"(2) A={A:.0e} vs modified baseline: |dM/M|={devs[A]:.2e}")
    ratio = devs[1e-6] / devs[1e-9] if devs[1e-9] > 0 else float('inf')
    print(f"    linear-scaling ratio (1e-6/1e-9) = {ratio:.0f} "
          f"(expect ~1000) -> {'PASS' if 500 < ratio < 2000 else 'CHECK'}")

    Mc, Rc = solve_star_op(eos, PC, 0.0087, ALPHA, embedding=emb)
    print(f"(3) contraction A=0.0087 vs modified baseline: "
          f"dM/M={100*(Mc/Mb-1):+.2f}%, dR/R={100*(Rc/Rb-1):+.2f}% "
          f"(same path, only A varied)")
