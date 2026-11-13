"""
Companion-paper support: the solar bound at the scalar-tensor companion's
alpha = 1 extension point.

The minimal chameleon companion model (a scalar-tensor companion to the RDT
series I-V, NOT part of the load-bearing chain; it reproduces the density law
rather than deriving the canonical divergence-operator mechanism) reproduces
the series saturating opening law with density exponent alpha = 1
(phi_eq ~ beta0*rho/(m^2*Mpl) => eps ~ rho/(rho + rho_*)). This alpha = 1 is an
EXTENSION POINT outside the canonical range alpha in [0.02, 0.5] scanned in
Paper I (which bounds A < 0.0087). This script reruns the Paper I V3 analysis
(exact modified
Lane-Emden solar models, Borexino flux ratios with B16-GS98 theory errors
in quadrature, helioseismic proxy penalty) at alpha = 1 and reports the
95% CL upper limit on the amplitude A.

Result (the alpha = 1 variant bound, reported separately from the canonical
A < 0.0087): A_95(alpha=1) = 0.0120  =>  eps_sat = 3A < 0.036,
d_eff > 2.964 at saturation (distinct from the canonical floor 2.974).

Uses solar_model.py (the Paper I V3 model module, copied verbatim from the
solarcore_code repository). Runs in ~30 s single-threaded.
"""

import numpy as np

from solar_model import Params, calibrate, flux_ratios, sound_speed_deviation

FAST = dict(rtol=1e-9, atol=1e-11, max_step=0.05)
SIGMA_HELIO = 2e-3
EXP_R = {"pp": 1.02, "pep": 0.88, "be7": 1.01, "b8": 1.04}
SIG = {"pp": 0.10, "pep": 0.15,
       "be7": np.hypot(0.03, 0.06), "b8": np.hypot(0.08, 0.12)}
ALPHA = 1.0


def chi2_combined(A, std):
    if A == 0.0:
        return sum(((1.0 - EXP_R[b]) / SIG[b]) ** 2 for b in EXP_R)
    mod = calibrate(Params(A, ALPHA), **FAST)
    r = flux_ratios(mod, std)
    _, dcs = sound_speed_deviation(mod, std)
    chi2_nu = sum(((r[b] - EXP_R[b]) / SIG[b]) ** 2 for b in EXP_R)
    return chi2_nu + (np.max(np.abs(dcs)) / SIGMA_HELIO) ** 2


def main():
    std = calibrate(Params(0.0, ALPHA), **FAST)
    A_grid = np.arange(0.0, 0.0201, 0.001)
    prof = np.array([chi2_combined(A, std) for A in A_grid])
    j = int(prof.argmin())
    target = prof[j] + 2.71  # one-sided 95% CL
    a95 = np.nan
    for i in range(j, len(A_grid) - 1):
        if prof[i] <= target <= prof[i + 1]:
            f = (target - prof[i]) / (prof[i + 1] - prof[i])
            a95 = A_grid[i] + f * (A_grid[i + 1] - A_grid[i])
            break
    print(f"alpha = {ALPHA}: chi2 min = {prof[j]:.3f} at A = {A_grid[j]:.3f}")
    print(f"A_95(alpha=1) = {a95:.4f}  (95% CL, one-sided)")
    print(f"=> eps_sat = 3*A_95 = {3*a95:.4f};  d_eff > {3*(1-a95):.4f} at saturation")


if __name__ == "__main__":
    main()
