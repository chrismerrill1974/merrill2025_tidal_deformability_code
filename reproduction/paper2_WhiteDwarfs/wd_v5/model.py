"""
WD V5 — RDT-modified white dwarf structure (Paper 2 revision).

Same canonical RDT law and modified-operator equation as Paper 1 V3
(solarcore_v3/model.py), applied to n = 3/2 polytropic white dwarfs:

    Omega(rho) = 1 - A*x/(1+x),  x = (rho/rho0)^alpha,  rho0 = 150 g/cm^3
    delta(rho) = 3*(Omega - 1) <= 0
    theta'' + [ (2+delta)/xi + (d delta/d xi) ln(xi) ] theta' + theta^n = 0

with delta = delta(rho_c * theta^n), reference scale r_star = a (the
Lane-Emden length), physical polytropic constant for a non-relativistic
degenerate electron gas (mu_e = 2):

    P = K rho^(5/3),  K = 1.0036e13 / mu_e^(5/3)  [cgs]

This fixes the V4 implementation, which used K = 1 in cgs (making ~10 m
"stars") and then started integration OUTSIDE the stellar surface
(r0 = 1e3 cm), inflating masses 9-79x (verified; see V3_REVISIONS.md G4).

No solar-style calibration here: white dwarfs are a one-parameter family
in rho_c. For each rho_c we return (M, R) for standard and RDT models.
"""

from dataclasses import dataclass
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from numerical_guardrails import check_problem, check_surface, positive

import numpy as np
from scipy.integrate import solve_ivp, simpson

N_POLY = 1.5
MU_E = 2.0

G_CGS = 6.67430e-8
MSUN_CGS = 1.98847e33
K_CGS = 1.0036e13 / MU_E ** (5.0 / 3.0)
RHO0_CGS = 150.0

# SI equivalents for the unit-invariance gate
G_SI = 6.67430e-11
MSUN_SI = 1.98847e30
# P = K rho^(5/3): [K] = Pa / (kg/m^3)^(5/3); cgs->SI factor = 0.1 / 1000^(5/3)
K_SI = K_CGS * 0.1 / 1000.0 ** (5.0 / 3.0)
RHO0_SI = 150.0e3

# Surface regularization (same role as Paper 1 V3): d(delta)/d(theta)
# ~ theta^(n*alpha - 1) diverges integrably at the surface for alpha < 2/3.
THETA_FLOOR = 1e-8


@dataclass(frozen=True)
class Params:
    """RDT density-law parameters. A=0 recovers the standard polytrope."""
    A: float
    alpha: float


def consts(si=False):
    if si:
        return G_SI, MSUN_SI, K_SI, RHO0_SI
    return G_CGS, MSUN_CGS, K_CGS, RHO0_CGS


def delta_of_theta(theta, rho_c, p, rho0):
    th = np.clip(theta, 1e-30, None)
    x = (rho_c * th**N_POLY / rho0) ** p.alpha
    return -3.0 * p.A * x / (1.0 + x)


def ddelta_dtheta(theta, rho_c, p, rho0):
    th = np.clip(theta, 1e-30, None)
    x = (rho_c * th**N_POLY / rho0) ** p.alpha
    return -3.0 * N_POLY * p.A * p.alpha * x / (th * (1.0 + x) ** 2)


def lane_emden(p, rho_c, rho0, xi0=1e-6, xi_max=10.0, rtol=1e-10, atol=1e-12,
               max_step=0.01):
    """Exact modified Lane-Emden (n = 3/2). Returns (xi, theta, dtheta)."""

    check_problem(p, rho_c, rho0, xi0, xi_max, rtol, atol, max_step)

    def rhs(xi, y):
        theta, thetap = y
        th = max(theta, 0.0)
        d = delta_of_theta(th, rho_c, p, rho0)
        dd_dxi = (ddelta_dtheta(th, rho_c, p, rho0) * thetap
                  if th > THETA_FLOOR else 0.0)
        coeff = (2.0 + d) / xi + dd_dxi * np.log(xi)
        return [thetap, -(th**N_POLY) - coeff * thetap]

    d0 = delta_of_theta(1.0, rho_c, p, rho0)
    theta0 = 1.0 - xi0**2 / (2.0 * (3.0 + d0))
    thetap0 = -xi0 / (3.0 + d0)

    def surface(xi, y):
        return y[0]
    surface.terminal = True
    surface.direction = -1.0

    sol = solve_ivp(rhs, (xi0, xi_max), [theta0, thetap0], events=surface,
                    max_step=max_step, rtol=rtol, atol=atol, method="DOP853")
    check_surface(sol)
    return sol.t, sol.y[0], sol.y[1]


@dataclass
class WDSolution:
    p: Params
    rho_c: float
    xi1: float
    a_le: float
    M: float      # g (or kg in SI mode)
    R: float      # cm (or m)


def build_wd(p, rho_c, si=False, **le_kwargs):
    """Solve + scale one white dwarf of given central density."""
    positive(rho_c=rho_c)
    G, _, K, rho0 = consts(si)
    a2 = (N_POLY + 1.0) * K * rho_c ** (1.0 / N_POLY - 1.0) / (4.0 * np.pi * G)
    a_le = np.sqrt(a2)
    xi, theta, dtheta = lane_emden(p, rho_c, rho0, **le_kwargs)
    xi1 = xi[-1]
    I = simpson(xi**2 * np.clip(theta, 0.0, None) ** N_POLY, x=xi)
    M = 4.0 * np.pi * a_le**3 * rho_c * I
    R = a_le * xi1
    positive(mass=M, radius=R, integral=I)
    return WDSolution(p, rho_c, xi1, a_le, M, R)


def mr_pair(p, rho_c, si=False, **le_kwargs):
    """(standard, modified) solutions at the same central density."""
    std = build_wd(Params(0.0, p.alpha), rho_c, si=si, **le_kwargs)
    mod = build_wd(p, rho_c, si=si, **le_kwargs)
    return std, mod
