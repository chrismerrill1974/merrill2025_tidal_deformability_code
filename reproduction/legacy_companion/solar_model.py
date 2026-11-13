"""
SolarCore V3 — RDT-modified Lane-Emden model with density-law dimensionality.

Theory (see V3_REVISIONS.md sections B, G2):
  Opening fraction (series-canonical law, shared with Papers 2-3):
      Omega(rho)  = 1 - A * x/(1+x),   x = (rho/rho0)^alpha
      d_eff(rho)  = 3 * Omega(rho)
      delta(rho)  = d_eff - 3 = -3A * x/(1+x)   (<= 0)

  Only the gravitational divergence operator carries d_eff; physical space
  stays 3D for volume elements, mass integrals, and reaction kernels.
  Reference scale r0 = alpha_LE (the Lane-Emden length), so the operator is
  xi^(2+delta) and all logarithms are of the dimensionless xi -> results are
  unit-invariant by construction.

  Exact modified Lane-Emden equation (n = 3), no log expansion:
      theta'' + [ (2+delta)/xi + (d delta/d xi) * ln(xi) ] theta' + theta^3 = 0
  with delta = delta(rho_c * theta^3) so d delta/d xi = (d delta/d theta) theta'.

  Series IC near the center:  theta ~ 1 - xi^2 / (2*(3+delta_0)).

  Calibration to (M_sun, R_sun) is self-consistent in rho_c (delta depends on
  rho_c), solved by fixed-point iteration. Central temperature from ideal gas
  + polytropic calibration: T_c proportional to G * alpha_LE^2 * rho_c
  (NOT rho_c^(1/3) at fixed K — K differs between calibrated models).
"""

from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp
from scipy.integrate import simpson

N_POLY = 3

# Surface regularization: for alpha < 1/3 the density law has
# d(delta)/d(theta) ~ theta^(3*alpha - 1) -> integrable divergence at the
# surface, where delta itself is ~1e-5 and irrelevant. Below this theta the
# delta-derivative term is dropped (validated: gate 3b checks xi_1 is
# insensitive to the floor at the 1e-7 level).
THETA_FLOOR = 1e-8

# Solar constants. CGS by default; pass si=True to constants() for the
# unit-invariance validation gate.
M_SUN_CGS = 1.98847e33   # g
R_SUN_CGS = 6.957e10     # cm
M_SUN_SI = 1.98847e30    # kg
R_SUN_SI = 6.957e8       # m
RHO0_CGS = 150.0         # g/cm^3, series reference density (Papers 2-3)
RHO0_SI = 150.0e3        # kg/m^3


@dataclass(frozen=True)
class Params:
    """RDT density-law parameters. A=0 recovers the standard polytrope."""
    A: float
    alpha: float


def units(si=False):
    if si:
        return M_SUN_SI, R_SUN_SI, RHO0_SI
    return M_SUN_CGS, R_SUN_CGS, RHO0_CGS


# ---------------------------------------------------------------------------
# Dimensionality law
# ---------------------------------------------------------------------------

def delta_of_theta(theta, rho_c, p, rho0):
    """delta = d_eff - 3 evaluated at rho = rho_c * theta^3."""
    th = np.clip(theta, 1e-30, None)
    x = (rho_c * th**3 / rho0) ** p.alpha
    return -3.0 * p.A * x / (1.0 + x)


def ddelta_dtheta(theta, rho_c, p, rho0):
    """d(delta)/d(theta) = -9 A alpha x / (theta (1+x)^2)."""
    th = np.clip(theta, 1e-30, None)
    x = (rho_c * th**3 / rho0) ** p.alpha
    return -9.0 * p.A * p.alpha * x / (th * (1.0 + x) ** 2)


# ---------------------------------------------------------------------------
# Modified Lane-Emden solver
# ---------------------------------------------------------------------------

def lane_emden(p, rho_c, rho0, xi0=1e-6, xi_max=20.0, rtol=1e-10, atol=1e-12,
               max_step=0.01, log_offset=0.0):
    """
    Solve the exact modified Lane-Emden equation for given rho_c (which sets
    the density scale inside delta). Returns (xi, theta, dtheta) arrays ending
    at the surface xi_1 (theta = 0).

    log_offset L implements the reference-scale choice r0 in
    r^delta = (r/r0)^delta: the operator log becomes ln(xi) - L.
    L = 0 -> r0 = a (Lane-Emden scale, default); L = ln(xi_1) -> r0 = R_sun.
    """

    def rhs(xi, y):
        theta, thetap = y
        th = max(theta, 0.0)
        d = delta_of_theta(th, rho_c, p, rho0)
        dd_dxi = (ddelta_dtheta(th, rho_c, p, rho0) * thetap
                  if th > THETA_FLOOR else 0.0)
        coeff = (2.0 + d) / xi + dd_dxi * (np.log(xi) - log_offset)
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
    if sol.t_events[0].size == 0:
        raise RuntimeError(f"no surface found for {p}")
    return sol.t, sol.y[0], sol.y[1]


# ---------------------------------------------------------------------------
# Solar calibration (self-consistent in rho_c)
# ---------------------------------------------------------------------------

@dataclass
class Solution:
    p: Params
    xi: np.ndarray
    theta: np.ndarray
    dtheta: np.ndarray
    xi1: float
    a_le: float      # Lane-Emden length scale alpha_LE = R/xi1
    rho_c: float
    I_mass: float    # integral xi^2 theta^3 dxi
    Tc_unit: float   # T_c in units of (G mu m_u / k_B): equals a_le^2 * rho_c * pi... see note

    def kernel_integral(self, gamma):
        """integral xi^2 theta^gamma dxi over the solution."""
        th = np.clip(self.theta, 0.0, None)
        return simpson(self.xi**2 * th**gamma, x=self.xi)


def calibrate(p, si=False, tol=1e-12, max_iter=60, r0_mode="a", **le_kwargs):
    """
    Solve + calibrate to (M_sun, R_sun) with fixed-point iteration on rho_c.
    Standard model (A=0) converges in one pass (delta independent of rho_c).
    Extra kwargs are passed to lane_emden (solver tolerances for scans).

    r0_mode: "a" (default, reference scale = Lane-Emden length) or "Rsun"
    (reference scale = solar radius; log_offset = ln(xi_1), iterated to
    self-consistency along with rho_c).
    """
    M, R, rho0 = units(si)
    rho_c = 80.0 * (rho0 / 150.0)  # initial guess ~ standard solar polytrope
    log_offset = np.log(6.897) if r0_mode == "Rsun" else 0.0
    last = None
    for _ in range(max_iter):
        xi, theta, dtheta = lane_emden(p, rho_c, rho0,
                                       log_offset=log_offset, **le_kwargs)
        if r0_mode == "Rsun":
            log_offset = np.log(xi[-1])
        xi1 = xi[-1]
        a_le = R / xi1
        I_mass = simpson(xi**2 * np.clip(theta, 0.0, None) ** 3, x=xi)
        rho_c_new = M / (4.0 * np.pi * a_le**3 * I_mass)
        if last is not None and abs(rho_c_new - rho_c) / rho_c < tol:
            rho_c = rho_c_new
            break
        last = rho_c
        rho_c = rho_c_new
    # T_c = (mu m_u / k_B) * K * rho_c^{1/3}, with K = pi G a_le^2 rho_c^{2/3}
    # => T_c proportional to G * a_le^2 * rho_c. Store the unit-free part.
    Tc_unit = a_le**2 * rho_c
    return Solution(p, xi, theta, dtheta, xi1, a_le, rho_c, I_mass, Tc_unit)


# ---------------------------------------------------------------------------
# Observables: flux ratios, sound speed
# ---------------------------------------------------------------------------

# Temperature exponents beta_i (Bahcall-style, evaluated near solar T_c) and
# kernel exponents gamma_i = 6 + beta_i (rho^2 T^beta, two-body, rho ~ theta^3).
BETA = {"pp": 4.0, "pep": 2.0, "be7": 10.0, "b8": 20.0}

# Experimental flux ratios (data / SSM), verified 2026-06-10:
# Borexino 2018 (Nature 562, 505) Table 2 over B16-GS98 (Vinyoles+ 2017,
# ApJ 835, 202) Table 6. pep is HZ-conditioned. Sigmas are MEASUREMENT-ONLY;
# analyze.py adds SSM theory errors in quadrature for the headline statistic.
EXP_R = {"pp": 1.02, "pep": 0.88, "be7": 1.01, "b8": 1.04}
EXP_SIGMA = {"pp": 0.10, "pep": 0.15, "be7": 0.03, "b8": 0.08}


def flux_ratios(mod, std):
    """
    R_i = Phi_i(mod)/Phi_i(std) with full calibration prefactors:
        Phi_i ~ rho_c^2 * a_le^3 * T_c^beta_i * Integral(xi^2 theta^(6+beta_i))
    """
    out = {}
    pref_rho = (mod.rho_c / std.rho_c) ** 2
    pref_vol = (mod.a_le / std.a_le) ** 3
    tc = mod.Tc_unit / std.Tc_unit
    for br, beta in BETA.items():
        gamma = 6.0 + beta
        Ii = mod.kernel_integral(gamma) / std.kernel_integral(gamma)
        out[br] = pref_rho * pref_vol * tc**beta * Ii
    return out


def chi2_fluxes(ratios):
    chi2 = sum(((ratios[b] - EXP_R[b]) / EXP_SIGMA[b]) ** 2 for b in BETA)
    return chi2


def sound_speed_deviation(mod, std, x_lo=0.01, x_hi=0.98, n=400):
    """
    Fractional deviation of c_s vs fractional radius. c_s^2 ~ T ~ T_c theta.
    Returns (x, dcs_over_cs).
    """
    x = np.linspace(x_lo, x_hi, n)
    th_s = np.interp(x * std.xi1, std.xi, std.theta)
    th_m = np.interp(x * mod.xi1, mod.xi, mod.theta)
    ok = (th_s > 0) & (th_m > 0)
    tc = mod.Tc_unit / std.Tc_unit
    ratio = np.sqrt(tc * th_m[ok] / th_s[ok])
    return x[ok], ratio - 1.0
