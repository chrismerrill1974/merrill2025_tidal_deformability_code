"""
paper6 decisive check: is the divergence operator's stellar-regime content a
density-dependent effective coupling G_eff(rho)?

We integrate TOV with dP/dr multiplied by a LOCAL density factor mult(rho):
  GR        : mult = 1
  G_eff(op) : mult = 1/(1 - A x/(1+x))   = 3/(3+delta), the operator field
              ratio a_op/a_N evaluated locally; > 1 => stronger gravity =>
              CONTRACTION.
  F-mult    : mult = (3 Omega - 1)/2     = (d_eff-1)/2, the V1-V7 ansatz;
              < 1 => weaker gravity => EXPANSION (the published wrong sign).

If G_eff(op) reproduces Paper III V8's operator result (M_max -0.71%,
R(1.4) -0.42%) and is OPPOSITE to the F-multiplier (+0.89%), then the
operator's measurable content IS G_eff(rho) = G/(1 - A x/(1+x)). A consistent
LOCAL field-theory completion then only needs to deliver this factor with the
correct (>1) sign; the residual genuinely-dimensional r-scaling piece is the
negligible ln(r/r*) term (III V8: ~0.15 pp).

Gate: A=0 (mult=1) reproduces the GR baseline (M_max ~ 2.0495, R(1.4) ~ 11.69).
Reuses the validated ns_v7 EOS/TOV machinery.  Runs in seconds, single-thread.
"""
import os
import sys

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

NS = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                   "..", "paper3_NeutronStars", "ns_v7"))
sys.path.insert(0, NS)
from eos import CompOSE_SLy4                      # noqa: E402
from tov import G, C, MSUN, KM, RHO0_RDT          # noqa: E402

A, ALPHA = 0.0087, 0.122


def _x(rho):
    return (rho / RHO0_RDT) ** ALPHA


def gr(rho):
    return 1.0


def geff_op(rho):
    # operator field ratio a_op/a_N = 3/(3+delta), delta = -3A x/(1+x)
    return 1.0 / (1.0 - A * _x(rho) / (1.0 + _x(rho)))


def fmult(rho):
    omega = 1.0 - A * _x(rho) / (1.0 + _x(rho))
    return (3.0 * omega - 1.0) / 2.0


def solve_mult(eos, P_c, mult, r_max=3.0e6, rtol=1e-8, max_step=2.0e4):
    eps_c = eos.eps_of_P(P_c)
    r0 = 1.0e2
    m0 = 4.0 / 3.0 * np.pi * r0**3 * eps_c / C**2
    P0 = P_c - (2.0 * np.pi * G / 3.0) * (eps_c + P_c) \
        * (eps_c + 3.0 * P_c) / C**4 * r0**2
    P_surf = eos.P_min * 1.001

    def rhs(r, y):
        P, m = y
        if P <= P_surf:
            return [0.0, 0.0]
        eps = eos.eps_of_P(P)
        denom = 1.0 - 2.0 * G * m / (r * C**2)
        if denom <= 0:
            return [0.0, 0.0]
        dP = -G * (eps + P) * (m + 4.0 * np.pi * r**3 * P / C**2) \
            / (C**2 * r**2 * denom)
        dP *= mult(eos.rho_of_P(P))
        dm = 4.0 * np.pi * r**2 * eps / C**2
        return [dP, dm]

    def surf(r, y):
        return y[0] - P_surf
    surf.terminal = True
    surf.direction = -1.0

    sol = solve_ivp(rhs, (r0, r_max), [P0, m0], events=surf, rtol=rtol,
                    atol=[P_surf * 1e-3, 1e25], max_step=max_step,
                    method="DOP853")
    if sol.t_events[0].size == 0:
        raise RuntimeError("no surface")
    return sol.y_events[0][0][1] / MSUN, sol.t_events[0][0] / KM


def mmax(eos, mult):
    lps = np.linspace(34.6, 36.1, 30)
    Ms = [(_safe(eos, 10.0**lp, mult)) for lp in lps]
    return np.nanmax(Ms)


def _safe(eos, Pc, mult):
    try:
        return solve_mult(eos, Pc, mult)[0]
    except RuntimeError:
        return np.nan


def r_at(eos, mult, Mt=1.4):
    f = lambda lp: solve_mult(eos, 10.0**lp, mult)[0] - Mt
    lp = brentq(f, 34.5, 35.7, xtol=1e-6)
    return solve_mult(eos, 10.0**lp, mult)[1]


if __name__ == "__main__":
    eos = CompOSE_SLy4()
    print(f"saturated factors (A={A}): G_eff(op)={geff_op(1e15):.5f} (>1, "
          f"contraction); F-mult={fmult(1e15):.5f} (<1, expansion)")
    Mg, Rg = mmax(eos, gr), r_at(eos, gr)
    print(f"\n{'model':>10} {'M_max':>8} {'dMmax%':>8} {'R(1.4)km':>9} {'dR%':>7}")
    for name, mult in [("GR", gr), ("G_eff(op)", geff_op), ("F-mult", fmult)]:
        Mm, Rr = mmax(eos, mult), r_at(eos, mult)
        print(f"{name:>10} {Mm:8.4f} {100*(Mm/Mg-1):+8.2f} {Rr:9.3f} "
              f"{100*(Rr/Rg-1):+7.2f}")
    print("\nExpect: G_eff(op) ~ -0.7%/-0.4% (= III V8 operator, contraction);"
          " F-mult ~ +0.9% (expansion). If so, operator content = G_eff(rho).")
