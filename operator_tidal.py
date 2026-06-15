"""
Paper IV (V6) — tidal deformability on the CONTRACTED operator backgrounds.

Real Hinderer/PPL Love-number integration (reuses the validated tidal
pipeline of this repo: love_k2, cs2_of_P, constants from tidal.py), with the
canonical operator modification (Embedding 1): the field ratio
R_op = a_op/a_N multiplies BOTH the pressure gradient AND the metric
potential nu', consistently (hydrostatic balance dP/dr = -(eps+P) nu'/2).
Minimal coupling: GR perturbation equation on the modified background, as in
IV V5.

A=0 reduces to GR tidal (gate: Lambda(1.4) ~ 298, matching the V5 GR pipeline
and the literature ~297). Reports the COMPUTED Lambda reduction and the
k2-vs-R decomposition: does dLambda/Lambda exceed the R^5 estimate (5 dR/R)
because k2 drops too?

This is the V6 mechanism (operator form, contraction). The superseded V5
F-multiplier driver scripts are in archive_v5/.
"""
import os
import sys

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tidal import love_k2, cs2_of_P, G, C, MSUN, KM, KG   # noqa: E402
from eos import CompOSE_SLy4                               # noqa: E402
from tov_operator import delta, ddelta_drho, _drho_dP      # noqa: E402


def solve_op_tidal(eos, P_c, A, alpha, rstar=1.0e6, rtol=1e-8, max_step=2.0e4):
    eps_c = eos.eps_of_P(P_c)
    rho_c = eos.rho_of_P(P_c)
    r0 = 1.0e2
    m0 = 4.0 / 3.0 * np.pi * r0**3 * eps_c / C**2
    P0 = P_c - (2.0 * np.pi * G / 3.0) * (eps_c + P_c) \
        * (eps_c + 3.0 * P_c) / C**4 * r0**2
    d0 = delta(rho_c, A, alpha)
    aN0 = 4.0 * np.pi * G * rho_c / 3.0 * r0
    aop0 = 4.0 * np.pi * G * rho_c / (3.0 + d0) * r0
    P_surf = eos.P_min * 1.001

    def rhs(r, s):
        P, m, a_op, a_N, Y = s
        if P <= P_surf or a_N <= 0:
            return [0.0, 0.0, 0.0, 0.0, 0.0]
        eps = eos.eps_of_P(P); rho = eos.rho_of_P(P)
        Pg, eg = KG * P, KG * eps
        mg = G * m / C**2
        denom = 1.0 - 2.0 * mg / r
        if denom <= 0:
            return [0.0, 0.0, 0.0, 0.0, 0.0]
        elam = 1.0 / denom
        Rop = a_op / a_N if A != 0.0 else 1.0
        base = (mg + 4.0 * np.pi * r**3 * Pg) * elam / r**2
        dPg = -(eg + Pg) * base * Rop
        nu_p = 2.0 * base * Rop                      # modified metric potential
        cs2 = cs2_of_P(eos, P)
        Q = (4.0 * np.pi * elam * (5.0 * eg + 9.0 * Pg + (eg + Pg) / cs2)
             - 6.0 * elam / r**2 - nu_p**2)
        dY = -(Y**2 + Y * elam * (1.0 + 4.0 * np.pi * r**2 * (Pg - eg))
               + r**2 * Q) / r
        dP_phys = dPg * C**4 / G
        dm = 4.0 * np.pi * r**2 * eps / C**2
        src = 4.0 * np.pi * G * rho
        d = delta(rho, A, alpha)
        dprime = ddelta_drho(rho, A, alpha) * _drho_dP(eos, P) * dP_phys
        da_N = src - (2.0 / r) * a_N
        da_op = src - ((2.0 + d) / r + dprime * np.log(r / rstar)) * a_op
        return [dP_phys, dm, da_op, da_N, dY]

    def surface(r, s):
        return s[0] - P_surf
    surface.terminal = True; surface.direction = -1.0

    sol = solve_ivp(rhs, (r0, 3.0e6), [P0, m0, aop0, aN0, 2.0], events=surface,
                    rtol=rtol, atol=[P_surf * 1e-3, 1e25, 1e-3, 1e-3, 1e-8],
                    max_step=max_step, method="DOP853")
    if sol.t_events[0].size == 0:
        raise RuntimeError(f"no surface for P_c={P_c:.3e}")
    R = sol.t_events[0][0]; ev = sol.y_events[0][0]
    M, yR = ev[1], ev[4]
    eps_s = eos.eps_of_P(P_surf * 1.01)
    yR -= 4.0 * np.pi * R**3 * (KG * eps_s) / (G * M / C**2)
    Cc = G * M / (R * C**2)
    k2 = love_k2(Cc, yR)
    Lam = (2.0 / 3.0) * k2 / Cc**5
    return M / MSUN, R / KM, Cc, k2, Lam


def at_mass(eos, Mt, A, alpha):
    f = lambda lp: solve_op_tidal(eos, 10.0**lp, A, alpha)[0] - Mt
    lp = brentq(f, 34.5, 35.8, xtol=1e-7)
    return solve_op_tidal(eos, 10.0**lp, A, alpha)


if __name__ == "__main__":
    eos = CompOSE_SLy4()
    A, AL = 0.0087, 0.122
    print(f"{'M':>4} {'model':>5} {'R km':>8} {'k2':>8} {'Lambda':>9}")
    for Mt in (1.2, 1.4, 1.6):
        Mg, Rg, Cg, k2g, Lg = at_mass(eos, Mt, 0.0, AL)
        Mo, Ro, Co, k2o, Lo = at_mass(eos, Mt, A, AL)
        dR = Ro / Rg - 1; dk2 = k2o / k2g - 1; dL = Lo / Lg - 1
        print(f"{Mt:>4} {'GR':>5} {Rg:8.3f} {k2g:8.4f} {Lg:9.1f}")
        print(f"{Mt:>4} {'RDT':>5} {Ro:8.3f} {k2o:8.4f} {Lo:9.1f}"
              f"   dL/L={100*dL:+.2f}%  5dR/R={100*5*dR:+.2f}%  dk2/k2={100*dk2:+.2f}%")
    print("\n(reviewer prediction: dL/L ~ -4% to -5%, larger than R^5's ~-2%,"
          " because k2 drops with the contracted profile)")
