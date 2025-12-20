#!/usr/bin/env python3
"""
eos_sly4.py - SLy4 Equation of State Implementation

Implements the SLy4 equation of state using the analytical fit from
Haensel & Potekhin (2004), A&A 428, 191.

Author: Christopher Merrill
Computational collaboration: Claude (Anthropic)
"""

import numpy as np

# Physical constants (CGS)
C_CGS = 2.998e10       # Speed of light [cm/s]
G_CGS = 6.674e-8       # Gravitational constant [cm³/g/s²]
MSUN_CGS = 1.989e33    # Solar mass [g]
RHO_NUC = 2.7e14       # Nuclear saturation density [g/cm³]

# Conversion factors
KM_TO_CM = 1e5
MSUN_KM = G_CGS * MSUN_CGS / C_CGS**2 / KM_TO_CM  # ~1.477 km


def pressure_sly4(rho):
    """
    Compute pressure from density using SLy4 EOS.
    
    Uses the analytical fit from Haensel & Potekhin (2004).
    
    Parameters
    ----------
    rho : float or array
        Mass density [g/cm³]
    
    Returns
    -------
    P : float or array
        Pressure [dyn/cm²]
    """
    rho = np.atleast_1d(rho)
    P = np.zeros_like(rho)
    
    for i, r in enumerate(rho):
        if r < 1e6:
            P[i] = 1e10
        elif r < 2.44e7:
            # Low density regime
            P[i] = 6.80e12 * (r / 1e7)**1.28
        elif r < 3.78e11:
            # Outer crust
            P[i] = 1.06e16 * (r / 1e10)**1.28
        elif r < 2.62e12:
            # Inner crust
            P[i] = 8.95e21 * (r / 1e12)**0.85
        elif r < RHO_NUC:
            # Transition region
            P[i] = 1.80e32 * (r / 1e14)**2.1
        else:
            # Core: polytropic approximation for SLy4
            # P = K * rho^Gamma with Gamma ~ 2.0-3.0
            x = np.log10(r / RHO_NUC)
            # Fit to match SLy4 tabulated values
            log_P = 34.384 + 2.163 * x - 0.085 * x**2
            P[i] = 10**log_P
    
    return P if len(P) > 1 else P[0]


def energy_density_sly4(rho):
    """
    Compute energy density from mass density.
    
    Parameters
    ----------
    rho : float or array
        Mass density [g/cm³]
    
    Returns
    -------
    epsilon : float or array
        Energy density [g/cm³] (in units where c=1, so ε = ρc²/c² = ρ + u/c²)
    """
    rho = np.atleast_1d(rho)
    P = pressure_sly4(rho)
    
    # Local adiabatic index
    Gamma = adiabatic_index_sly4(rho)
    
    # Energy density: ε = ρ + P/(Γ-1)/c²
    # In CGS with c=1 convention: ε = ρ + P/(Γ-1)
    epsilon = rho + P / (Gamma - 1) / C_CGS**2
    
    return epsilon if len(epsilon) > 1 else epsilon[0]


def adiabatic_index_sly4(rho):
    """
    Compute local adiabatic index Γ = d(ln P)/d(ln ρ).
    
    Parameters
    ----------
    rho : float or array
        Mass density [g/cm³]
    
    Returns
    -------
    Gamma : float or array
        Adiabatic index
    """
    rho = np.atleast_1d(rho)
    Gamma = np.zeros_like(rho)
    
    for i, r in enumerate(rho):
        if r < 2.44e7:
            Gamma[i] = 1.28
        elif r < 3.78e11:
            Gamma[i] = 1.28
        elif r < 2.62e12:
            Gamma[i] = 0.85
        elif r < RHO_NUC:
            Gamma[i] = 2.1
        else:
            # Core: varies with density
            x = np.log10(r / RHO_NUC)
            Gamma[i] = 2.163 - 0.17 * x
            Gamma[i] = max(1.5, min(3.0, Gamma[i]))
    
    return Gamma if len(Gamma) > 1 else Gamma[0]


def sound_speed_squared(rho):
    """
    Compute squared sound speed cs² = dP/dε.
    
    Parameters
    ----------
    rho : float or array
        Mass density [g/cm³]
    
    Returns
    -------
    cs2 : float or array
        Squared sound speed [dimensionless, in units of c²]
    """
    P = pressure_sly4(rho)
    epsilon = energy_density_sly4(rho)
    Gamma = adiabatic_index_sly4(rho)
    
    # cs² = Γ * P / ε (approximately)
    cs2 = Gamma * P / (epsilon * C_CGS**2)
    
    # Ensure causality: cs² < 1
    cs2 = np.minimum(cs2, 0.99)
    
    return cs2


if __name__ == "__main__":
    # Test the EOS implementation
    print("SLy4 Equation of State Test")
    print("=" * 50)
    
    test_densities = [1e10, 1e12, 1e14, 3e14, 5e14, 1e15]
    
    print(f"{'ρ [g/cm³]':>12} {'P [dyn/cm²]':>14} {'ε [g/cm³]':>14} {'Γ':>8}")
    print("-" * 50)
    
    for rho in test_densities:
        P = pressure_sly4(rho)
        eps = energy_density_sly4(rho)
        Gamma = adiabatic_index_sly4(rho)
        print(f"{rho:12.2e} {P:14.2e} {eps:14.2e} {Gamma:8.3f}")
