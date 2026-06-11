"""
NS V7 — TOV solver (GR baseline and RDT-modified).

Standard TOV:
    dP/dr = -G (eps + P)(m + 4 pi r^3 P/c^2) / [c^2 r^2 (1 - 2Gm/(r c^2))]
    dm/dr = 4 pi r^2 eps / c^2

RDT (Paper III ansatz, documented as such): dP/dr multiplied by
F(rho) = (d_eff(rho) - 1)/2 with the series opening law
Omega(rho) = 1 - A x/(1+x), x = (rho/rho0)^alpha, rho0 = 150 g/cm^3.
At NS densities the law is saturated for all allowed alpha, so
F ~ 1 - (3/2) A = const: an interior G_eff rescaling.

Series IC near center (uniform-density expansion), not a seed blob:
    m(r0) = (4/3) pi r0^3 eps_c / c^2,
    P(r0) = P_c - (2 pi G / 3) (eps_c/c^2 + P_c/c^2)
            (eps_c/c^2 + 3 P_c/c^2) ... (standard expansion); r0 chosen
    so the IC terms are O(1e-8) relative.
"""

import numpy as np
from scipy.integrate import solve_ivp

G = 6.67430e-8
C = 2.99792458e10
MSUN = 1.98847e33
KM = 1.0e5

RHO0_RDT = 150.0


def F_rdt(rho, A, alpha):
    x = (rho / RHO0_RDT) ** alpha
    omega = 1.0 - A * x / (1.0 + x)
    return (3.0 * omega - 1.0) / 2.0


def solve_star(eos, P_c, A=0.0, alpha=0.2, r_max=3.0e6,
               rtol=1e-8, atol_P=None, max_step=2.0e4):
    """Integrate one star. Returns (M [Msun], R [km])."""
    eps_c = eos.eps_of_P(P_c)
    r0 = 1.0e2  # cm; star is ~1e6 cm, so r0/R ~ 1e-4 and IC error ~ (r0/R)^2
    m0 = 4.0 / 3.0 * np.pi * r0**3 * eps_c / C**2
    # P(r0) from the series expansion
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
        if A != 0.0:
            dP *= F_rdt(eos.rho_of_P(P), A, alpha)
        dm = 4.0 * np.pi * r**2 * eps / C**2
        return [dP, dm]

    def surface(r, y):
        return y[0] - P_surf
    surface.terminal = True
    surface.direction = -1.0

    sol = solve_ivp(rhs, (r0, r_max), [P0, m0], events=surface,
                    rtol=rtol, atol=[P_surf * 1e-3, 1e25],
                    max_step=max_step, method="DOP853")
    if sol.t_events[0].size == 0:
        raise RuntimeError(f"no surface for P_c={P_c:.3e}")
    R = sol.t_events[0][0]
    M = sol.y_events[0][0][1]
    return M / MSUN, R / KM


def mr_curve(eos, P_c_grid, A=0.0, alpha=0.2, **kw):
    out = []
    for P_c in P_c_grid:
        try:
            M, R = solve_star(eos, P_c, A=A, alpha=alpha, **kw)
            out.append((P_c, M, R))
        except RuntimeError:
            out.append((P_c, np.nan, np.nan))
    return np.array(out)
