#!/usr/bin/env python3
"""
tov_solver.py - TOV Equation Solver for GR and RDT

Solves the Tolman-Oppenheimer-Volkoff equations for neutron star structure
in both standard General Relativity and Recursive Dimensionality Theory.

Author: Christopher Merrill
Computational collaboration: Claude (Anthropic)
"""

import numpy as np
from scipy.integrate import solve_ivp
from eos_sly4 import pressure_sly4, energy_density_sly4, RHO_NUC, C_CGS, G_CGS, MSUN_CGS

# Geometrized unit conversions
KM_TO_CM = 1e5
MSUN_KM = G_CGS * MSUN_CGS / C_CGS**2 / KM_TO_CM  # ~1.477 km


def F_rdt(rho, alpha):
    """
    RDT geometric correction factor.
    
    F(ρ) = (d_eff - 1) / 2
    
    where d_eff = 3 + α * log10(ρ/ρ_0) for ρ > ρ_0
    
    Parameters
    ----------
    rho : float
        Mass density [g/cm³]
    alpha : float
        RDT coupling parameter
    
    Returns
    -------
    F : float
        Geometric correction factor (F=1 for standard GR)
    """
    if rho <= RHO_NUC:
        return 1.0
    else:
        d_eff = 3.0 + alpha * np.log10(rho / RHO_NUC)
        return (d_eff - 1.0) / 2.0


def tov_equations_gr(r, y, eos_P, eos_eps):
    """
    Standard GR TOV equations.
    
    dm/dr = 4π r² ε
    dP/dr = -(ε + P)(m + 4π r³ P) / [r(r - 2m)]
    
    Parameters
    ----------
    r : float
        Radial coordinate [cm]
    y : array
        State vector [m, P] where m is in grams, P in dyn/cm²
    eos_P : callable
        Pressure as function of density
    eos_eps : callable
        Energy density as function of density
    
    Returns
    -------
    dy : array
        Derivatives [dm/dr, dP/dr]
    """
    m, P = y
    
    if P <= 0 or r <= 0:
        return [0.0, 0.0]
    
    # Invert EOS to get density from pressure (simple Newton iteration)
    rho = invert_eos(P, eos_P)
    epsilon = eos_eps(rho)
    
    # Convert to geometrized units for TOV
    r_geom = r / KM_TO_CM  # km
    m_geom = m / MSUN_CGS * MSUN_KM  # km
    P_geom = P * G_CGS / C_CGS**4 * KM_TO_CM**2  # 1/km²
    eps_geom = epsilon * G_CGS / C_CGS**2 * KM_TO_CM**2  # 1/km²
    
    # TOV equations in geometrized units
    factor = r_geom - 2 * m_geom
    if factor <= 0:
        return [0.0, 0.0]
    
    dm_dr_geom = 4 * np.pi * r_geom**2 * eps_geom
    dP_dr_geom = -(eps_geom + P_geom) * (m_geom + 4*np.pi * r_geom**3 * P_geom) / (r_geom * factor)
    
    # Convert back to CGS
    dm_dr = dm_dr_geom * MSUN_CGS / MSUN_KM * KM_TO_CM  # g/cm
    dP_dr = dP_dr_geom * C_CGS**4 / G_CGS / KM_TO_CM  # dyn/cm²/cm
    
    return [dm_dr, dP_dr]


def tov_equations_rdt(r, y, eos_P, eos_eps, alpha):
    """
    RDT-modified TOV equations.
    
    Same as GR but with F(ρ) correction to pressure gradient:
    dP/dr = -(ε + P)(m + 4π r³ P) / [r(r - 2m)] × F(ρ)
    
    Parameters
    ----------
    r : float
        Radial coordinate [cm]
    y : array
        State vector [m, P]
    eos_P : callable
        Pressure as function of density
    eos_eps : callable
        Energy density as function of density
    alpha : float
        RDT coupling parameter
    
    Returns
    -------
    dy : array
        Derivatives [dm/dr, dP/dr]
    """
    m, P = y
    
    if P <= 0 or r <= 0:
        return [0.0, 0.0]
    
    # Invert EOS to get density
    rho = invert_eos(P, eos_P)
    epsilon = eos_eps(rho)
    
    # RDT correction factor
    F = F_rdt(rho, alpha)
    
    # Convert to geometrized units
    r_geom = r / KM_TO_CM
    m_geom = m / MSUN_CGS * MSUN_KM
    P_geom = P * G_CGS / C_CGS**4 * KM_TO_CM**2
    eps_geom = epsilon * G_CGS / C_CGS**2 * KM_TO_CM**2
    
    factor = r_geom - 2 * m_geom
    if factor <= 0:
        return [0.0, 0.0]
    
    dm_dr_geom = 4 * np.pi * r_geom**2 * eps_geom
    dP_dr_geom = -(eps_geom + P_geom) * (m_geom + 4*np.pi * r_geom**3 * P_geom) / (r_geom * factor)
    
    # Apply RDT correction
    dP_dr_geom *= F
    
    # Convert back to CGS
    dm_dr = dm_dr_geom * MSUN_CGS / MSUN_KM * KM_TO_CM
    dP_dr = dP_dr_geom * C_CGS**4 / G_CGS / KM_TO_CM
    
    return [dm_dr, dP_dr]


def invert_eos(P, eos_P, rho_min=1e6, rho_max=1e16, tol=1e-8):
    """
    Invert EOS to get density from pressure using bisection.
    
    Parameters
    ----------
    P : float
        Pressure [dyn/cm²]
    eos_P : callable
        Pressure as function of density
    
    Returns
    -------
    rho : float
        Mass density [g/cm³]
    """
    if P <= 0:
        return rho_min
    
    # Bisection search
    rho_lo, rho_hi = rho_min, rho_max
    
    for _ in range(100):
        rho_mid = np.sqrt(rho_lo * rho_hi)  # Geometric mean for log-spaced
        P_mid = eos_P(rho_mid)
        
        if abs(P_mid - P) / P < tol:
            return rho_mid
        
        if P_mid < P:
            rho_lo = rho_mid
        else:
            rho_hi = rho_mid
    
    return rho_mid


def solve_tov(rho_c, alpha=None, r_max=50e5, n_points=10000):
    """
    Solve TOV equations for given central density.
    
    Parameters
    ----------
    rho_c : float
        Central density [g/cm³]
    alpha : float or None
        RDT coupling parameter. If None, solve standard GR.
    r_max : float
        Maximum radius [cm]
    n_points : int
        Number of output points
    
    Returns
    -------
    result : dict
        Dictionary containing:
        - r : radial coordinates [cm]
        - m : enclosed mass [g]
        - P : pressure [dyn/cm²]
        - rho : density [g/cm³]
        - R : stellar radius [km]
        - M : total mass [M_sun]
    """
    # Initial conditions
    P_c = pressure_sly4(rho_c)
    r0 = 100.0  # Start at small radius to avoid r=0 singularity
    m0 = 4/3 * np.pi * r0**3 * rho_c  # Small initial mass
    
    y0 = [m0, P_c]
    
    # Surface pressure threshold
    P_surf = 1e10  # dyn/cm²
    
    # Event function to detect surface
    def surface_event(r, y):
        return y[1] - P_surf
    surface_event.terminal = True
    surface_event.direction = -1
    
    # Choose equation set
    if alpha is None:
        tov_func = lambda r, y: tov_equations_gr(r, y, pressure_sly4, energy_density_sly4)
    else:
        tov_func = lambda r, y: tov_equations_rdt(r, y, pressure_sly4, energy_density_sly4, alpha)
    
    # Integrate
    r_span = (r0, r_max)
    r_eval = np.linspace(r0, r_max, n_points)
    
    sol = solve_ivp(tov_func, r_span, y0, method='RK45', 
                    t_eval=r_eval, events=surface_event,
                    max_step=1e4, rtol=1e-8, atol=1e-10)
    
    # Extract results
    r = sol.t
    m = sol.y[0]
    P = sol.y[1]
    
    # Compute density profile
    rho = np.array([invert_eos(p, pressure_sly4) if p > 0 else 0 for p in P])
    
    # Find surface
    surface_idx = np.argmax(P < P_surf) if np.any(P < P_surf) else -1
    if surface_idx > 0:
        R = r[surface_idx] / KM_TO_CM  # km
        M = m[surface_idx] / MSUN_CGS  # solar masses
    else:
        R = r[-1] / KM_TO_CM
        M = m[-1] / MSUN_CGS
    
    return {
        'r': r[:surface_idx] if surface_idx > 0 else r,
        'm': m[:surface_idx] if surface_idx > 0 else m,
        'P': P[:surface_idx] if surface_idx > 0 else P,
        'rho': rho[:surface_idx] if surface_idx > 0 else rho,
        'R': R,
        'M': M,
        'rho_c': rho_c
    }


def find_mass(M_target, alpha=None, rho_min=1e14, rho_max=3e15, tol=1e-4):
    """
    Find central density that produces target mass using bisection.
    
    Parameters
    ----------
    M_target : float
        Target mass [solar masses]
    alpha : float or None
        RDT coupling parameter
    
    Returns
    -------
    result : dict
        TOV solution for target mass
    """
    for _ in range(50):
        rho_mid = np.sqrt(rho_min * rho_max)
        result = solve_tov(rho_mid, alpha)
        M = result['M']
        
        if abs(M - M_target) / M_target < tol:
            return result
        
        if M < M_target:
            rho_min = rho_mid
        else:
            rho_max = rho_mid
    
    return result


if __name__ == "__main__":
    print("TOV Solver Test")
    print("=" * 60)
    
    # Test for 1.4 solar mass star
    M_target = 1.4
    
    print(f"\nFinding {M_target} M_sun neutron star...")
    
    # GR solution
    result_gr = find_mass(M_target, alpha=None)
    print(f"\nGR Solution:")
    print(f"  M = {result_gr['M']:.4f} M_sun")
    print(f"  R = {result_gr['R']:.3f} km")
    print(f"  ρ_c = {result_gr['rho_c']:.3e} g/cm³")
    
    # RDT solution
    result_rdt = find_mass(M_target, alpha=0.30)
    print(f"\nRDT Solution (α = 0.30):")
    print(f"  M = {result_rdt['M']:.4f} M_sun")
    print(f"  R = {result_rdt['R']:.3f} km")
    print(f"  ρ_c = {result_rdt['rho_c']:.3e} g/cm³")
    
    print(f"\nRDT Enhancement:")
    print(f"  ΔR/R = {100*(result_rdt['R']/result_gr['R'] - 1):+.2f}%")
