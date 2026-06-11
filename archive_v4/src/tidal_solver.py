#!/usr/bin/env python3
"""
tidal_solver.py - Relativistic Tidal Perturbation Solver

Solves the Hinderer (2008) equations for tidal perturbations
to compute the Love number k₂ and tidal deformability Λ.

Reference: Hinderer, T. (2008). ApJ 677, 1216.

Author: Christopher Merrill
Computational collaboration: Claude (Anthropic)
"""

import numpy as np
from scipy.integrate import solve_ivp
from eos_sly4 import energy_density_sly4, adiabatic_index_sly4, C_CGS, G_CGS

# Unit conversions
KM_TO_CM = 1e5


def tidal_y_equation(r, y_val, r_array, m_array, P_array, rho_array):
    """
    The Hinderer (2008) y(r) equation for tidal perturbations.
    
    r dy/dr + y² + y*F₁(r) + F₂(r) = 0
    
    which we rewrite as:
    dy/dr = -(y² + y*F₁ + F₂) / r
    
    Parameters
    ----------
    r : float
        Radial coordinate [cm]
    y_val : float
        Current value of y
    r_array, m_array, P_array, rho_array : arrays
        Interpolation data from TOV solution
    
    Returns
    -------
    dy_dr : float
        Derivative dy/dr
    """
    if r <= 0:
        return 0.0
    
    # Interpolate background quantities
    m = np.interp(r, r_array, m_array)
    P = np.interp(r, r_array, P_array)
    rho = np.interp(r, r_array, rho_array)
    
    if P <= 0 or rho <= 0:
        return 0.0
    
    # Energy density
    epsilon = energy_density_sly4(rho)
    
    # Convert to geometrized units (km)
    r_km = r / KM_TO_CM
    m_km = m / 1.989e33 * 1.477  # M_sun to km
    
    # Pressure and energy density in geometric units (1/km²)
    P_geom = P * G_CGS / C_CGS**4 * KM_TO_CM**2
    eps_geom = epsilon * G_CGS / C_CGS**2 * KM_TO_CM**2
    
    # Sound speed squared (dP/dε)
    Gamma = adiabatic_index_sly4(rho)
    cs2 = Gamma * P_geom / eps_geom if eps_geom > 0 else 0.1
    cs2 = min(cs2, 0.99)  # Causality bound
    
    # Metric function
    factor = r_km - 2 * m_km
    if factor <= 0:
        return 0.0
    
    # F₁ coefficient (Eq. 12 in Hinderer 2008)
    F1 = (r_km - 4 * np.pi * r_km**3 * (eps_geom - P_geom)) / factor
    
    # F₂ coefficient (Eq. 13 in Hinderer 2008)
    Q = 4 * np.pi * r_km**2 * (5*eps_geom + 9*P_geom + (eps_geom + P_geom)/cs2) - 6
    Q /= factor
    
    term2 = 4 * m_km**2 / (r_km**2 * factor**2)
    term2 *= (1 + 4*np.pi * r_km**3 * P_geom / m_km)**2 if m_km > 0 else 1
    
    F2 = Q - term2
    
    # The ODE: r dy/dr = -(y² + y*F₁ + F₂)
    dy_dr = -(y_val**2 + y_val * F1 + F2) / r_km * KM_TO_CM
    
    return dy_dr


def solve_tidal(r_array, m_array, P_array, rho_array):
    """
    Solve the tidal perturbation equation on a TOV background.
    
    Parameters
    ----------
    r_array : array
        Radial coordinates [cm]
    m_array : array
        Enclosed mass [g]
    P_array : array
        Pressure [dyn/cm²]
    rho_array : array
        Density [g/cm³]
    
    Returns
    -------
    y_R : float
        Value of y at the stellar surface
    y_array : array
        Full y(r) profile
    """
    # Starting point and boundary condition
    r0 = r_array[1]  # Small but non-zero radius
    y0 = 2.0  # Central boundary condition
    
    # Integration range
    R = r_array[-1]
    
    # Solve ODE
    def ode_func(r, y):
        return [tidal_y_equation(r, y[0], r_array, m_array, P_array, rho_array)]
    
    sol = solve_ivp(ode_func, (r0, R), [y0], 
                    method='RK45', 
                    t_eval=r_array[1:],
                    rtol=1e-8, atol=1e-10)
    
    y_array = np.concatenate([[y0], sol.y[0]])
    y_R = sol.y[0][-1]
    
    return y_R, y_array


def compute_k2(M, R, y_R):
    """
    Compute the tidal Love number k₂ from surface values.
    
    Uses Eq. 23-24 from Hinderer (2008).
    
    Parameters
    ----------
    M : float
        Stellar mass [solar masses]
    R : float
        Stellar radius [km]
    y_R : float
        Value of y at stellar surface
    
    Returns
    -------
    k2 : float
        Tidal Love number
    """
    # Compactness
    C = M * 1.477 / R  # M in km / R in km
    
    # Hinderer (2008) Eq. 23-24
    num = 8 * C**5 / 5 * (1 - 2*C)**2 * (2 + 2*C*(y_R - 1) - y_R)
    
    denom = (2*C * (6 - 3*y_R + 3*C*(5*y_R - 8)) +
             4*C**3 * (13 - 11*y_R + C*(3*y_R - 2) + 2*C**2*(1 + y_R)) +
             3*(1 - 2*C)**2 * (2 - y_R + 2*C*(y_R - 1)) * np.log(1 - 2*C))
    
    k2 = num / denom
    
    return k2


def compute_Lambda(k2, M, R):
    """
    Compute dimensionless tidal deformability Λ.
    
    Λ = (2/3) k₂ (R/M)⁵
    
    Parameters
    ----------
    k2 : float
        Love number
    M : float
        Mass [solar masses]
    R : float
        Radius [km]
    
    Returns
    -------
    Lambda : float
        Dimensionless tidal deformability
    """
    # Convert M to km
    M_km = M * 1.477
    
    Lambda = (2/3) * k2 * (R / M_km)**5
    
    return Lambda


def compute_Lambda_tilde(m1, m2, Lambda1, Lambda2):
    """
    Compute combined tidal deformability for binary system.
    
    Eq. (5) from Abbott et al. (2017).
    
    Parameters
    ----------
    m1, m2 : float
        Component masses [solar masses]
    Lambda1, Lambda2 : float
        Component tidal deformabilities
    
    Returns
    -------
    Lambda_tilde : float
        Combined tidal deformability
    """
    M_tot = m1 + m2
    
    term1 = (m1 + 12*m2) * m1**4 * Lambda1
    term2 = (m2 + 12*m1) * m2**4 * Lambda2
    
    Lambda_tilde = 16/13 * (term1 + term2) / M_tot**5
    
    return Lambda_tilde


if __name__ == "__main__":
    from tov_solver import find_mass
    
    print("Tidal Solver Test")
    print("=" * 60)
    
    # Get TOV solution for 1.4 solar mass star
    M_target = 1.4
    
    print(f"\nComputing tidal deformability for {M_target} M_sun star...")
    
    # GR
    result_gr = find_mass(M_target, alpha=None)
    y_R_gr, _ = solve_tidal(result_gr['r'], result_gr['m'], 
                            result_gr['P'], result_gr['rho'])
    k2_gr = compute_k2(result_gr['M'], result_gr['R'], y_R_gr)
    Lambda_gr = compute_Lambda(k2_gr, result_gr['M'], result_gr['R'])
    
    print(f"\nGR Results:")
    print(f"  R = {result_gr['R']:.3f} km")
    print(f"  C = {result_gr['M']*1.477/result_gr['R']:.4f}")
    print(f"  y_R = {y_R_gr:.4f}")
    print(f"  k₂ = {k2_gr:.4f}")
    print(f"  Λ = {Lambda_gr:.1f}")
    
    # RDT
    result_rdt = find_mass(M_target, alpha=0.30)
    y_R_rdt, _ = solve_tidal(result_rdt['r'], result_rdt['m'],
                             result_rdt['P'], result_rdt['rho'])
    k2_rdt = compute_k2(result_rdt['M'], result_rdt['R'], y_R_rdt)
    Lambda_rdt = compute_Lambda(k2_rdt, result_rdt['M'], result_rdt['R'])
    
    print(f"\nRDT Results (α = 0.30):")
    print(f"  R = {result_rdt['R']:.3f} km")
    print(f"  C = {result_rdt['M']*1.477/result_rdt['R']:.4f}")
    print(f"  y_R = {y_R_rdt:.4f}")
    print(f"  k₂ = {k2_rdt:.4f}")
    print(f"  Λ = {Lambda_rdt:.1f}")
    
    print(f"\nRDT Enhancement:")
    print(f"  ΔR/R = {100*(result_rdt['R']/result_gr['R'] - 1):+.2f}%")
    print(f"  Δk₂/k₂ = {100*(k2_rdt/k2_gr - 1):+.2f}%")
    print(f"  ΔΛ/Λ = {100*(Lambda_rdt/Lambda_gr - 1):+.2f}%")
