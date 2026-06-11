"""
Tidal V5 — relativistic tidal deformability on ns_v7 backgrounds.

Hinderer (2008) / Postnikov-Prakash-Lattimer (2010) formalism: integrate the
quadrupole metric-perturbation variable y(r) = r H'/H alongside the TOV
equations, extract the Love number k2 at the surface, and
Lambda = (2/3) k2 C^-5.

Minimal-coupling assumption (carried over from Paper 4 V4, where it was one
of the few sound elements): the RDT F-multiplier modifies the equilibrium
background only; the linearized perturbation equations keep their GR form
on that background, and the metric potential nu'(r) is computed from the
standard Einstein relation (not from the modified dP/dr).

Background EOS/TOV reused from Paper 3's ns_v7 package (imported, not
copied). Geometrized internally: P_g = GP/c^4, eps_g = G eps/c^4,
m_g = Gm/c^2; y is dimensionless.
"""

import os
import sys

import numpy as np
from scipy.integrate import solve_ivp

HERE = os.path.dirname(os.path.abspath(__file__))
from tov import F_rdt, G, C, MSUN, KM  # noqa: E402

KG = G / C**4   # converts cgs P or eps to geometrized [cm^-2]


def cs2_of_P(eos, P, rel=1e-4):
    """Sound speed squared dP/deps (units of c^2), finite difference on the
    EOS interpolant."""
    P1, P2 = P * (1 - rel), P * (1 + rel)
    e1, e2 = eos.eps_of_P(P1), eos.eps_of_P(P2)
    return (P2 - P1) / max(e2 - e1, 1e-300)


def solve_tidal(eos, P_c, A=0.0, alpha=0.2, r_max=3.0e6,
                rtol=1e-8, max_step=2.0e4):
    """
    Integrate [P, m, y] from the series-expansion center to the surface.
    Returns (M [Msun], R [km], C, k2, Lambda).
    """
    eps_c = eos.eps_of_P(P_c)
    r0 = 1.0e2  # cm
    m0 = 4.0 / 3.0 * np.pi * r0**3 * eps_c / C**2
    P0 = P_c - (2.0 * np.pi * G / 3.0) * (eps_c + P_c) \
        * (eps_c + 3.0 * P_c) / C**4 * r0**2
    y0 = 2.0

    P_surf = eos.P_min * 1.001

    def rhs(r, s):
        P, m, y = s
        if P <= P_surf:
            return [0.0, 0.0, 0.0]
        eps = eos.eps_of_P(P)
        Pg, eg = KG * P, KG * eps
        mg = G * m / C**2
        denom = 1.0 - 2.0 * mg / r
        if denom <= 0:
            return [0.0, 0.0, 0.0]
        elam = 1.0 / denom
        dPg = -(eg + Pg) * (mg + 4.0 * np.pi * r**3 * Pg) * elam / r**2
        if A != 0.0:
            dPg *= F_rdt(eos.rho_of_P(P), A, alpha)
        # metric potential from the standard Einstein relation
        nu_p = 2.0 * elam * (mg + 4.0 * np.pi * r**3 * Pg) / r**2
        cs2 = cs2_of_P(eos, P)
        Q = (4.0 * np.pi * elam * (5.0 * eg + 9.0 * Pg + (eg + Pg) / cs2)
             - 6.0 * elam / r**2 - nu_p**2)
        dy = -(y**2 + y * elam * (1.0 + 4.0 * np.pi * r**2 * (Pg - eg))
               + r**2 * Q) / r
        return [dPg * C**4 / G, 4.0 * np.pi * r**2 * eps / C**2, dy]

    def surface(r, s):
        return s[0] - P_surf
    surface.terminal = True
    surface.direction = -1.0

    sol = solve_ivp(rhs, (r0, r_max), [P0, m0, y0], events=surface,
                    rtol=rtol, atol=[P_surf * 1e-3, 1e25, 1e-8],
                    max_step=max_step, method="DOP853")
    if sol.t_events[0].size == 0:
        raise RuntimeError(f"no surface for P_c={P_c:.3e}")
    R = sol.t_events[0][0]
    M = sol.y_events[0][0][1]
    yR = sol.y_events[0][0][2]

    # surface energy-density correction for tabulated EOS (eps does not
    # vanish at the surface cut): y_R -> y_R - 4 pi R^3 eps_s / (M c^2)
    eps_s = eos.eps_of_P(P_surf * 1.01)
    yR -= 4.0 * np.pi * R**3 * (KG * eps_s) / (G * M / C**2)

    Cc = G * M / (R * C**2)
    k2 = love_k2(Cc, yR)
    Lam = (2.0 / 3.0) * k2 / Cc**5
    return M / MSUN, R / KM, Cc, k2, Lam


def love_k2(C_, yR):
    """Hinderer (2008) k2(C, y_R)."""
    c1 = 1.0 - 2.0 * C_
    num = (8.0 / 5.0) * C_**5 * c1**2 * (2.0 + 2.0 * C_ * (yR - 1.0) - yR)
    den = (2.0 * C_ * (6.0 - 3.0 * yR + 3.0 * C_ * (5.0 * yR - 8.0))
           + 4.0 * C_**3 * (13.0 - 11.0 * yR + C_ * (3.0 * yR - 2.0)
                            + 2.0 * C_**2 * (1.0 + yR))
           + 3.0 * c1**2 * (2.0 - yR + 2.0 * C_ * (yR - 1.0))
           * np.log(c1))
    return num / den


def lambda_tilde(m1, m2, L1, L2):
    """Mass-weighted binary tidal deformability (Flanagan-Hinderer)."""
    return (16.0 / 13.0) * ((m1 + 12.0 * m2) * m1**4 * L1
                            + (m2 + 12.0 * m1) * m2**4 * L2) \
        / (m1 + m2)**5


class Polytrope1:
    """n=1 relativistic polytrope (P = K rho^2, eps = rho c^2 + P) for the
    Hinderer benchmark gate. Mimics the eos.py interface."""

    def __init__(self, K=4.0e4):
        self.K = K

    @property
    def P_min(self):
        return 1.0e20

    def rho_of_P(self, P):
        return np.sqrt(np.maximum(P, 1e-300) / self.K)

    def eps_of_P(self, P):
        return self.rho_of_P(P) * C**2 + P
