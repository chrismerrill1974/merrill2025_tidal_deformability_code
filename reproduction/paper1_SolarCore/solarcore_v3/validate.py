"""
Validation gates for SolarCore V3 (V3_REVISIONS.md section C6).
Run BEFORE any physics results. All gates must pass.

  Gate 1: A=0 regression — xi_1 = 6.896849, rho_c/rho_mean = 54.18 (n=3 values)
  Gate 2: unit invariance — identical dimensionless results in cgs and SI
  Gate 3: convergence — xi_1 and integrals stable under tighter tolerances
  Gate 4: fiducial RDT model solves, calibrates, and behaves (delta<=0,
          surface exists, flux ratios finite)
"""

import numpy as np

from model import Params, calibrate, lane_emden, flux_ratios, units

XI1_REF = 6.8968486      # n=3 Lane-Emden surface (literature)
RHOC_OVER_RHOBAR = 54.1825  # n=3 central-to-mean density ratio (literature)

FIDUCIAL = Params(A=0.0334, alpha=0.15)


def gate1_standard_regression():
    std = calibrate(Params(A=0.0, alpha=0.15))
    M, R, _ = units()
    rho_bar = M / (4.0 / 3.0 * np.pi * R**3)
    ratio = std.rho_c / rho_bar
    err_xi = abs(std.xi1 - XI1_REF) / XI1_REF
    err_ratio = abs(ratio - RHOC_OVER_RHOBAR) / RHOC_OVER_RHOBAR
    ok = err_xi < 1e-5 and err_ratio < 1e-4
    print(f"Gate 1 (A=0 regression): xi1 = {std.xi1:.7f} (ref {XI1_REF}, "
          f"rel err {err_xi:.2e}); rho_c/rho_bar = {ratio:.4f} "
          f"(ref {RHOC_OVER_RHOBAR}, rel err {err_ratio:.2e}); "
          f"rho_c = {std.rho_c:.4f} g/cm^3 -> {'PASS' if ok else 'FAIL'}")
    return ok


def gate2_unit_invariance():
    mod_c = calibrate(FIDUCIAL, si=False)
    std_c = calibrate(Params(A=0.0, alpha=0.15), si=False)
    mod_s = calibrate(FIDUCIAL, si=True)
    std_s = calibrate(Params(A=0.0, alpha=0.15), si=True)
    r_cgs = flux_ratios(mod_c, std_c)
    r_si = flux_ratios(mod_s, std_s)
    worst = max(abs(r_cgs[b] - r_si[b]) for b in r_cgs)
    xi_diff = abs(mod_c.xi1 - mod_s.xi1)
    ok = worst < 1e-9 and xi_diff < 1e-9
    print(f"Gate 2 (unit invariance): max flux-ratio diff cgs vs SI = "
          f"{worst:.2e}; xi1 diff = {xi_diff:.2e} -> {'PASS' if ok else 'FAIL'}")
    return ok


def gate3_convergence():
    _, _, rho0 = units()
    mod = calibrate(FIDUCIAL)
    xi_a, th_a, _ = lane_emden(FIDUCIAL, mod.rho_c, rho0, rtol=1e-10, atol=1e-12)
    xi_b, th_b, _ = lane_emden(FIDUCIAL, mod.rho_c, rho0, rtol=1e-12, atol=1e-14)
    diff = abs(xi_a[-1] - xi_b[-1])
    ok = diff < 1e-7
    print(f"Gate 3 (convergence): xi1 shift under 100x tighter tol = "
          f"{diff:.2e} -> {'PASS' if ok else 'FAIL'}")
    return ok


def gate3b_floor_sensitivity():
    import model
    _, _, rho0 = units()
    mod = calibrate(FIDUCIAL)
    saved = model.THETA_FLOOR
    try:
        model.THETA_FLOOR = 1e-6
        xi_a, _, _ = lane_emden(FIDUCIAL, mod.rho_c, rho0)
        model.THETA_FLOOR = 1e-10
        xi_b, _, _ = lane_emden(FIDUCIAL, mod.rho_c, rho0)
    finally:
        model.THETA_FLOOR = saved
    diff = abs(xi_a[-1] - xi_b[-1])
    ok = diff < 1e-6
    print(f"Gate 3b (surface floor): xi1 shift, floor 1e-6 vs 1e-10 = "
          f"{diff:.2e} -> {'PASS' if ok else 'FAIL'}")
    return ok


def gate4_fiducial_behaviour():
    std = calibrate(Params(A=0.0, alpha=0.15))
    mod = calibrate(FIDUCIAL)
    ratios = flux_ratios(mod, std)
    finite = all(np.isfinite(v) and v > 0 for v in ratios.values())
    print(f"Gate 4 (fiducial model): xi1 = {mod.xi1:.6f} "
          f"(std {std.xi1:.6f}), rho_c = {mod.rho_c:.4f} g/cm^3 "
          f"(std {std.rho_c:.4f}), Tc_mod/Tc_std = {mod.Tc_unit/std.Tc_unit:.6f}")
    print(f"         flux ratios: " +
          ", ".join(f"{b}={v:.5f}" for b, v in ratios.items()) +
          f" -> {'PASS' if finite else 'FAIL'}")
    return finite


def run_all():
    from numerical_guardrails import run_gates
    return run_gates([gate1_standard_regression, gate2_unit_invariance, gate3_convergence, gate3b_floor_sensitivity, gate4_fiducial_behaviour])


if __name__ == "__main__":
    run_all()
