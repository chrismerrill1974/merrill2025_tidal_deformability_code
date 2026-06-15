"""
Paper III keystone — relativistic divergence-operator modified TOV.

Canonical Papers I-II mechanism (NOT the III/IV F-multiplier): the gravitational
field obeys a (3+delta)-dimensional Gauss law, matter/mass kept 3D, delta<=0.
Two inequivalent relativistic embeddings (see derivation.md); both reduce to
exact V7 TOV at delta=0 and to the I-II Newtonian Gauss law in weak field.

Reuses eos.py (EOS) and constants/baseline from tov.py (this repo). CGS units
throughout: P,eps in dyne/cm^2 = erg/cm^3, rho in g/cm^3, m in g, r in cm.
A=0 routes to the exact V7 solver (guarantees the baseline).
"""
import os
import sys

import numpy as np
from scipy.integrate import solve_ivp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from tov import G, C, MSUN, KM, RHO0_RDT, solve_star as solve_star_v7  # noqa: E402


def delta(rho, A, alpha):
    """delta = d_eff - 3 = -3A x/(1+x) <= 0."""
    x = (rho / RHO0_RDT) ** alpha
    return -3.0 * A * x / (1.0 + x)


def ddelta_drho(rho, A, alpha):
    x = (rho / RHO0_RDT) ** alpha
    return -3.0 * A * alpha * x / (rho * (1.0 + x) ** 2)


def _drho_dP(eos, P):
    dP = 1e-3 * P
    return (eos.rho_of_P(P + dP) - eos.rho_of_P(P - dP)) / (2.0 * dP)


def solve_star_op(eos, P_c, A, alpha, embedding="emb1", rstar_cm=1.0e6,
                  r_max=3.0e6, rtol=1e-8, max_step=2.0e4):
    """Integrate one star with the operator modification. Returns (M[Msun], R[km])."""
    if A == 0.0:
        return solve_star_v7(eos, P_c, A=0.0, alpha=alpha,
                             r_max=r_max, rtol=rtol, max_step=max_step)

    eps_c = eos.eps_of_P(P_c)
    rho_c = eos.rho_of_P(P_c)
    r0 = 1.0e2
    m0 = 4.0 / 3.0 * np.pi * r0**3 * eps_c / C**2
    P0 = P_c - (2.0 * np.pi * G / 3.0) * (eps_c + P_c) \
        * (eps_c + 3.0 * P_c) / C**4 * r0**2
    P_surf = eos.P_min * 1.001
    d0 = delta(rho_c, A, alpha)

    if embedding == "emb1":
        # state = [P, m, a_op, a_N]   (a_* = Newtonian field, cm/s^2)
        aN0 = 4.0 * np.pi * G * rho_c / 3.0 * r0
        aop0 = 4.0 * np.pi * G * rho_c / (3.0 + d0) * r0
        y0 = [P0, m0, aop0, aN0]

        def rhs(r, y):
            P, m, a_op, a_N = y
            if P <= P_surf:
                return [0.0, 0.0, 0.0, 0.0]
            eps = eos.eps_of_P(P)
            rho = eos.rho_of_P(P)
            denom = 1.0 - 2.0 * G * m / (r * C**2)
            if denom <= 0 or a_N <= 0:
                return [0.0, 0.0, 0.0, 0.0]
            dP_tov = -G * (eps + P) * (m + 4.0 * np.pi * r**3 * P / C**2) \
                / (C**2 * r**2 * denom)
            R_op = a_op / a_N
            dP = dP_tov * R_op
            dm = 4.0 * np.pi * r**2 * eps / C**2
            src = 4.0 * np.pi * G * rho
            d = delta(rho, A, alpha)
            dprime = ddelta_drho(rho, A, alpha) * _drho_dP(eos, P) * dP
            da_N = src - (2.0 / r) * a_N
            da_op = src - ((2.0 + d) / r + dprime * np.log(r / rstar_cm)) * a_op
            return [dP, dm, da_op, da_N]

    elif embedding == "emb2":
        # state = [P, m, m_delta]   (effective-dimension weighted mass)
        w0 = (r0 / rstar_cm) ** d0
        md0 = 4.0 / 3.0 * np.pi * r0**3 * w0 * eps_c / C**2
        y0 = [P0, m0, md0]

        def rhs(r, y):
            P, m, m_delta = y
            if P <= P_surf:
                return [0.0, 0.0, 0.0]
            eps = eos.eps_of_P(P)
            rho = eos.rho_of_P(P)
            denom = 1.0 - 2.0 * G * m / (r * C**2)
            if denom <= 0:
                return [0.0, 0.0, 0.0]
            d = delta(rho, A, alpha)
            w = (r / rstar_cm) ** d
            dP = -G * (eps + P) * (m_delta + 4.0 * np.pi * r**3 * w * P / C**2) \
                / (C**2 * r**2 * w * denom)
            dm = 4.0 * np.pi * r**2 * eps / C**2
            dm_delta = 4.0 * np.pi * r**2 * w * eps / C**2
            return [dP, dm, dm_delta]
    else:
        raise ValueError(embedding)

    def surface(r, y):
        return y[0] - P_surf
    surface.terminal = True
    surface.direction = -1.0

    atol = [P_surf * 1e-3, 1e25] + [1e-3] * (len(y0) - 2)
    sol = solve_ivp(rhs, (r0, r_max), y0, events=surface,
                    rtol=rtol, atol=atol, max_step=max_step, method="DOP853")
    if sol.t_events[0].size == 0:
        raise RuntimeError(f"no surface for P_c={P_c:.3e} ({embedding})")
    R = sol.t_events[0][0]
    M = sol.y_events[0][0][1]
    return M / MSUN, R / KM


def mr_curve_op(eos, P_c_grid, A, alpha, embedding="emb1", **kw):
    out = []
    for P_c in P_c_grid:
        try:
            M, R = solve_star_op(eos, P_c, A, alpha, embedding=embedding, **kw)
            out.append((P_c, M, R))
        except RuntimeError:
            out.append((P_c, np.nan, np.nan))
    return np.array(out)
