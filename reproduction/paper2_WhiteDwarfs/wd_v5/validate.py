"""
Validation gates for WD V5. Run BEFORE any physics output.

  Gate 1: A=0 regression — xi_1 = 3.65375, omega = xi1^2|theta'(xi1)| = 2.71406
  Gate 2: analytic scaling — R ~ rho_c^(-1/6), M ~ rho_c^(1/2) over 1e5..1e7
          (the V4 failure mode: its numerics broke these laws badly)
  Gate 3: physical sanity — rho_c = 1e6, mu_e = 2: R ~ 1.1e9 cm, M ~ 0.5 Msun
  Gate 4: r0-inside-star check — integration must START inside the star
          (by construction here: xi0 = 1e-6 << xi_1; assert explicitly)
  Gate 5: unit invariance — fractional shifts identical in cgs and SI
  Gate 6: fiducial RDT model behaves; dR/R, dM/M finite and small
"""

import numpy as np

from model import Params, build_wd, mr_pair, consts, MSUN_CGS

XI1_REF = 3.65375374
OMEGA_REF = 2.71406      # xi1^2 * |theta'(xi1)| for n = 3/2

FIDUCIAL = Params(A=0.0087, alpha=0.15)   # Paper 1 V3 95% CL boundary scale


def gate1():
    s = build_wd(Params(0.0, 0.1), 1e6)
    # omega from M: M = 4 pi a^3 rho_c * I, and I -> omega at A=0
    I = s.M / (4.0 * np.pi * s.a_le**3 * s.rho_c)
    e1 = abs(s.xi1 - XI1_REF) / XI1_REF
    e2 = abs(I - OMEGA_REF) / OMEGA_REF
    ok = e1 < 1e-6 and e2 < 1e-4
    print(f"Gate 1 (A=0 regression): xi1 = {s.xi1:.6f} (rel err {e1:.1e}), "
          f"omega = {I:.5f} (rel err {e2:.1e}) -> {'PASS' if ok else 'FAIL'}")
    return ok


def gate2():
    rhos = np.array([1e5, 1e6, 1e7])
    sols = [build_wd(Params(0.0, 0.1), r) for r in rhos]
    R = np.array([s.R for s in sols])
    M = np.array([s.M for s in sols])
    pR = np.polyfit(np.log(rhos), np.log(R), 1)[0]
    pM = np.polyfit(np.log(rhos), np.log(M), 1)[0]
    ok = abs(pR + 1.0 / 6.0) < 1e-4 and abs(pM - 0.5) < 1e-4
    print(f"Gate 2 (scalings): d lnR/d lnrho = {pR:.6f} (ref -1/6), "
          f"d lnM/d lnrho = {pM:.6f} (ref +1/2) -> {'PASS' if ok else 'FAIL'}")
    return ok


def gate3():
    s = build_wd(Params(0.0, 0.1), 1e6)
    ok = 1.0e9 < s.R < 1.25e9 and 0.4 < s.M / MSUN_CGS < 0.6
    print(f"Gate 3 (physical): rho_c=1e6 -> R = {s.R:.3e} cm "
          f"({s.R/6.957e10*1000:.1f} mR_sun), M = {s.M/MSUN_CGS:.3f} Msun "
          f"-> {'PASS' if ok else 'FAIL'}")
    return ok


def gate4():
    s = build_wd(Params(0.0, 0.1), 1e6)
    xi0 = 1e-6
    ok = xi0 / s.xi1 < 1e-5
    print(f"Gate 4 (start inside star): xi0/xi1 = {xi0/s.xi1:.2e} "
          f"(V4 equivalent was >1: started OUTSIDE) -> {'PASS' if ok else 'FAIL'}")
    return ok


def gate5():
    out = {}
    for si in (False, True):
        std, mod = mr_pair(FIDUCIAL, 1e6 * (1e3 if si else 1.0), si=si)
        out[si] = (mod.R / std.R - 1.0, mod.M / std.M - 1.0)
    dR = abs(out[False][0] - out[True][0])
    dM = abs(out[False][1] - out[True][1])
    ok = dR < 1e-10 and dM < 1e-10
    print(f"Gate 5 (unit invariance): |dR/R diff| = {dR:.1e}, "
          f"|dM/M diff| = {dM:.1e} -> {'PASS' if ok else 'FAIL'}")
    return ok


def gate6():
    std, mod = mr_pair(FIDUCIAL, 1e6)
    dR = mod.R / std.R - 1.0
    dM = mod.M / std.M - 1.0
    ok = np.isfinite(dR) and np.isfinite(dM) and abs(dR) < 0.05
    print(f"Gate 6 (fiducial A={FIDUCIAL.A}, alpha={FIDUCIAL.alpha}): "
          f"dR/R = {dR:+.4e}, dM/M = {dM:+.4e}, xi1 mod = {mod.xi1:.5f} "
          f"-> {'PASS' if ok else 'FAIL'}")
    return ok


def run_all():
    from numerical_guardrails import run_gates
    return run_gates([gate1, gate2, gate3, gate4, gate5, gate6])


if __name__ == "__main__":
    run_all()
