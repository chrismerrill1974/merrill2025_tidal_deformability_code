#!/usr/bin/env python3
"""
mass_sequence.py - Generate Mass-Radius-Lambda Sequences

Computes full mass sequences for neutron stars in GR and RDT,
outputting M, R, k₂, Λ, and I for publication tables.

Author: Christopher Merrill
Computational collaboration: Claude (Anthropic)
"""

import numpy as np
import argparse
from tov_solver import find_mass
from tidal_solver import solve_tidal, compute_k2, compute_Lambda


def lattimer_schutz_I(M, R):
    """
    Approximate moment of inertia using Lattimer & Schutz (2005) formula.
    
    Parameters
    ----------
    M : float
        Mass [solar masses]
    R : float
        Radius [km]
    
    Returns
    -------
    I_bar : float
        Dimensionless moment of inertia I/M³
    """
    C = M * 1.477 / R  # Compactness
    I_bar = 0.237 / C**2 * (1 + 4.2*C + 90*C**4)
    return I_bar


def compute_sequence(masses, alpha=None, verbose=True):
    """
    Compute tidal sequence for given masses.
    
    Parameters
    ----------
    masses : list
        Target masses [solar masses]
    alpha : float or None
        RDT coupling parameter (None for GR)
    verbose : bool
        Print progress
    
    Returns
    -------
    results : list of dict
        Results for each mass point
    """
    theory = "RDT" if alpha else "GR"
    if verbose:
        print(f"\nComputing {theory} sequence...")
        print(f"{'M':>6} {'R':>8} {'C':>8} {'k2':>8} {'Lambda':>10} {'I_bar':>8}")
        print("-" * 56)
    
    results = []
    
    for M_target in masses:
        # Solve TOV
        tov = find_mass(M_target, alpha=alpha)
        
        # Solve tidal equations
        y_R, _ = solve_tidal(tov['r'], tov['m'], tov['P'], tov['rho'])
        
        # Compute Love number and Lambda
        k2 = compute_k2(tov['M'], tov['R'], y_R)
        Lambda = compute_Lambda(k2, tov['M'], tov['R'])
        
        # Compute moment of inertia
        I_bar = lattimer_schutz_I(tov['M'], tov['R'])
        
        # Compactness
        C = tov['M'] * 1.477 / tov['R']
        
        result = {
            'M': tov['M'],
            'R': tov['R'],
            'C': C,
            'k2': k2,
            'Lambda': Lambda,
            'I_bar': I_bar,
            'rho_c': tov['rho_c'],
            'y_R': y_R
        }
        results.append(result)
        
        if verbose:
            print(f"{result['M']:6.2f} {result['R']:8.3f} {result['C']:8.4f} "
                  f"{result['k2']:8.4f} {result['Lambda']:10.1f} {result['I_bar']:8.2f}")
    
    return results


def save_sequence(results, filename, alpha=None):
    """
    Save sequence to data file.
    """
    theory = f"RDT (alpha={alpha})" if alpha else "GR"
    
    with open(filename, 'w') as f:
        f.write(f"# Neutron Star Tidal Sequence - {theory}\n")
        f.write("# SLy4 Equation of State\n")
        f.write("# M[Msun]  R[km]    C        k2       Lambda    I_bar\n")
        
        for r in results:
            f.write(f"{r['M']:.3f}     {r['R']:.3f}   {r['C']:.5f}  "
                    f"{r['k2']:.5f}  {r['Lambda']:.1f}     {r['I_bar']:.3f}\n")
    
    print(f"Saved: {filename}")


def main():
    parser = argparse.ArgumentParser(description='Generate NS tidal sequences')
    parser.add_argument('--masses', nargs='+', type=float, 
                        default=[1.0, 1.2, 1.4, 1.6],
                        help='Target masses in solar masses')
    parser.add_argument('--alpha', type=float, default=None,
                        help='RDT coupling parameter (omit for GR)')
    parser.add_argument('--output', type=str, default='.',
                        help='Output directory')
    args = parser.parse_args()
    
    print("=" * 60)
    print("Neutron Star Tidal Deformability Sequence Generator")
    print("=" * 60)
    
    # Compute GR sequence
    results_gr = compute_sequence(args.masses, alpha=None)
    save_sequence(results_gr, f"{args.output}/gr_lambda_sequence.dat", alpha=None)
    
    # Compute RDT sequence if requested
    if args.alpha:
        results_rdt = compute_sequence(args.masses, alpha=args.alpha)
        save_sequence(results_rdt, f"{args.output}/rdt_lambda_sequence.dat", alpha=args.alpha)
        
        # Print comparison
        print("\n" + "=" * 60)
        print("RDT vs GR Comparison")
        print("=" * 60)
        print(f"{'M':>6} {'ΔR/R':>10} {'Δk2/k2':>10} {'ΔΛ/Λ':>10}")
        print("-" * 40)
        
        for gr, rdt in zip(results_gr, results_rdt):
            dR = 100 * (rdt['R'] / gr['R'] - 1)
            dk2 = 100 * (rdt['k2'] / gr['k2'] - 1)
            dL = 100 * (rdt['Lambda'] / gr['Lambda'] - 1)
            print(f"{gr['M']:6.2f} {dR:+10.2f}% {dk2:+10.2f}% {dL:+10.2f}%")


if __name__ == "__main__":
    main()
