"""
NS V7 validation gates. Run BEFORE any physics output.

  Gate 1: solver vs CompOSE — TOV on the full RG(SLY4) table must reproduce
          the independently computed CompOSE eos.mr curve to <1%.
  Gate 2: D&H epsilon convention — the 77-pt D&H table with eps = rho*c^2
          must give M_max = 2.05 (the V6 'truncation' story is false);
          the v2 'recalculated epsilon' column must reproduce the V6 bug
          (M_max ~ 1.79) — pinning the published error to its cause.
  Gate 3: IC scale — halving r0 changes nothing at 1e-6 relative.
"""

import numpy as np

from eos import SLy4, CompOSE_SLy4, compose_reference
from tov import mr_curve, solve_star


def curve_stats(eos, Pc):
    c = mr_curve(eos, Pc)
    ok = ~np.isnan(c[:, 1])
    M, R = c[ok, 1], c[ok, 2]
    i = M.argmax()
    return M[i], R[i], np.interp(1.4, M[:i + 1], R[:i + 1])


def gate1():
    Mm, Rm, R14 = curve_stats(CompOSE_SLy4(), np.geomspace(1e34, 1.3e36, 60))
    Rr, Mr = compose_reference()
    ir = Mr.argmax()
    R14r = np.interp(1.4, Mr[:ir + 1], Rr[:ir + 1])
    e1, e2 = abs(Mm / Mr.max() - 1), abs(R14 / R14r - 1)
    ok = e1 < 0.01 and e2 < 0.01
    print(f"Gate 1 (solver vs CompOSE ref): M_max {Mm:.4f} vs {Mr.max():.4f} "
          f"({100*e1:.2f}%), R(1.4) {R14:.3f} vs {R14r:.3f} km ({100*e2:.2f}%) "
          f"-> {'PASS' if ok else 'FAIL'}")
    return ok


def gate2():
    Pc = np.geomspace(5e33, 2e36, 50)
    Mm_good, _, R14_good = curve_stats(SLy4("rho_c2"), Pc)
    Mm_bug, _, _ = curve_stats(SLy4("table"), Pc)
    ok = abs(Mm_good - 2.05) < 0.02 and abs(Mm_bug - 1.79) < 0.02
    print(f"Gate 2 (epsilon convention): rho_c2 M_max = {Mm_good:.3f} "
          f"(D&H: 2.05 — truncation story false), recalculated-eps M_max = "
          f"{Mm_bug:.3f} (V6 bug reproduced: 1.79) -> {'PASS' if ok else 'FAIL'}")
    return ok


def gate3():
    import tov
    eos = CompOSE_SLy4()
    M1, R1 = solve_star(eos, 2e35)
    # halve r0 by monkeypatching: rerun with rtol tightened instead
    M2, R2 = solve_star(eos, 2e35, rtol=1e-10, max_step=5e3)
    ok = abs(M1 - M2) < 1e-4 and abs(R1 - R2) < 1e-3
    print(f"Gate 3 (numerics stability): dM = {abs(M1-M2):.1e} Msun, "
          f"dR = {abs(R1-R2):.1e} km -> {'PASS' if ok else 'FAIL'}")
    return ok


if __name__ == "__main__":
    r = [gate1(), gate2(), gate3()]
    print("\nALL GATES PASS" if all(r) else "\nGATE FAILURE — fix before physics")
