"""
Tidal V5 validation gates. Run BEFORE any RDT tidal output.

  Gate 1: Hinderer polytrope benchmark — n=1 polytrope, k2 ~ 0.26 at
          C = 0.15 (Hinderer 2008, corrected; Paper 4 V4's own check gave
          0.2598 with its independent implementation).
  Gate 2: background consistency — (M, R) from the tidal integration must
          match ns_v7 tov.solve_star to ~1e-4 (same physics, augmented
          state vector).
  Gate 3: published-SLy regression — Lambda(1.4) for SLy in GR must land in
          the literature band (~270-330; commonly quoted ~297). This is the
          gate Paper 4 V4 would have failed by a factor ~3 (it claimed
          'literature ~950-1000', an artifact of its broken EOS baseline).
  Gate 4: numerics — k2 stable under tightened tolerances.
"""

import numpy as np
from scipy.optimize import brentq

from tidal import solve_tidal, Polytrope1
from eos import CompOSE_SLy4
from tov import solve_star


def gate1():
    """n=1 polytrope, Newtonian limit: k2(C->0) = 0.2599 analytically.
    NOTE: Paper 4 V4 quoted 'k2 = 0.2598 at C = 0.15' as its validation —
    that is the NEWTONIAN value mislabeled; the relativistic k2 at C=0.15
    is ~0.071. We extrapolate C -> 0 from three small-C models."""
    eos = Polytrope1()
    Cs, k2s = [], []
    for Ct in (0.005, 0.01, 0.02):
        logPc = brentq(lambda lp: solve_tidal(eos, 10.0**lp)[2] - Ct,
                       29.0, 36.5, xtol=1e-7)
        _, _, Cc, k2, _ = solve_tidal(eos, 10.0**logPc)
        Cs.append(Cc)
        k2s.append(k2)
    k2_0 = np.polyfit(Cs, k2s, 1)[1]
    ok = abs(k2_0 - 0.2599) < 0.003
    print(f"Gate 1 (n=1 polytrope, Newtonian limit): extrapolated "
          f"k2(C->0) = {k2_0:.4f} (analytic 0.2599) "
          f"-> {'PASS' if ok else 'FAIL'}")
    return ok


def gate2():
    eos = CompOSE_SLy4()
    Pc = 2.0e35
    M1, R1 = solve_star(eos, Pc)
    M2, R2, _, _, _ = solve_tidal(eos, Pc)
    ok = abs(M1 - M2) < 1e-4 and abs(R1 - R2) < 1e-3
    print(f"Gate 2 (background consistency): dM = {abs(M1-M2):.1e} Msun, "
          f"dR = {abs(R1-R2):.1e} km -> {'PASS' if ok else 'FAIL'}")
    return ok


def gate3():
    eos = CompOSE_SLy4()

    def mass_diff(logPc):
        M, _, _, _, _ = solve_tidal(eos, 10.0**logPc)
        return M - 1.4

    logPc = brentq(mass_diff, 34.5, 35.8, xtol=1e-7)
    M, R, Cc, k2, Lam = solve_tidal(eos, 10.0**logPc)
    ok = 270.0 < Lam < 330.0
    print(f"Gate 3 (published-SLy regression): M = {M:.4f}, R = {R:.3f} km, "
          f"k2 = {k2:.4f}, Lambda(1.4) = {Lam:.1f} "
          f"(literature ~297) -> {'PASS' if ok else 'FAIL'}")
    return ok, 10.0**logPc


def gate4(Pc14):
    eos = CompOSE_SLy4()
    _, _, _, k2a, La = solve_tidal(eos, Pc14)
    _, _, _, k2b, Lb = solve_tidal(eos, Pc14, rtol=1e-10, max_step=5e3)
    ok = abs(k2a - k2b) < 1e-4 and abs(La / Lb - 1) < 1e-3
    print(f"Gate 4 (numerics): dk2 = {abs(k2a-k2b):.1e}, "
          f"dLambda/Lambda = {abs(La/Lb-1):.1e} -> {'PASS' if ok else 'FAIL'}")
    return ok


if __name__ == "__main__":
    r1 = gate1()
    r2 = gate2()
    r3, Pc14 = gate3()
    r4 = gate4(Pc14)
    print("\nALL GATES PASS" if all([r1, r2, r3, r4])
          else "\nGATE FAILURE — fix before physics")
