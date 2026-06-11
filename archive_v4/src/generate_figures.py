#!/usr/bin/env python3
"""
generate_figures.py - Generate Publication Figures

Creates the four main figures for Paper IV:
  1. Lambda(M) comparison with GW170817 constraint
  2. GW170817 prediction space (Lambda1-Lambda2 plane)
  3. I-Love universal relation
  4. R^5 geometric scaling verification

Author: Christopher Merrill
Computational collaboration: Claude (Anthropic)
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy.interpolate import interp1d
import argparse
import os

# Publication plot settings
plt.rcParams.update({
    'font.size': 11,
    'axes.labelsize': 12,
    'axes.titlesize': 12,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 10,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'font.family': 'serif',
    'mathtext.fontset': 'cm'
})


# Data from validated calculations
M_vals = np.array([1.0, 1.2, 1.4, 1.6])

GR = {
    'R': np.array([13.452, 13.502, 13.373, 13.020]),
    'k2': np.array([0.2449, 0.1750, 0.1304, 0.0979]),
    'Lambda': np.array([10241, 2995, 984.7, 331.5]),
    'I_bar': np.array([25.86, 18.97, 14.39, 11.06])
}

RDT = {
    'R': np.array([13.475, 13.603, 13.639, 13.604]),
    'k2': np.array([0.2484, 0.1842, 0.1466, 0.1224]),
    'Lambda': np.array([10483, 3274, 1221, 516.4]),
    'I_bar': np.array([25.94, 19.22, 14.87, 11.87])
}

# Derived quantities
delta_Lambda_pct = (RDT['Lambda'] - GR['Lambda']) / GR['Lambda'] * 100
delta_R_pct = (RDT['R'] - GR['R']) / GR['R'] * 100
delta_k2_pct = (RDT['k2'] - GR['k2']) / GR['k2'] * 100


def fig1_lambda_comparison(output_dir):
    """Figure 1: Lambda(M) comparison with GW170817 constraint."""
    
    fig, ax = plt.subplots(figsize=(7, 5.5))
    
    # Main panel
    ax.semilogy(M_vals, GR['Lambda'], 'bo-', markersize=10, linewidth=2,
                markeredgecolor='black', markeredgewidth=1, label='GR (SLy4)')
    ax.semilogy(M_vals, RDT['Lambda'], 'r^-', markersize=10, linewidth=2,
                markeredgecolor='black', markeredgewidth=1, label=r'RDT ($\alpha=0.30$)')
    
    # GW170817 constraint
    ax.axhspan(10, 720, alpha=0.15, color='green', label='GW170817 90% CI')
    ax.axhline(720, color='green', linestyle='--', linewidth=1.5, alpha=0.7)
    
    ax.set_xlabel(r'Mass $M$ [$M_\odot$]')
    ax.set_ylabel(r'Tidal Deformability $\Lambda$')
    ax.set_xlim(0.9, 1.7)
    ax.set_ylim(100, 20000)
    ax.legend(loc='upper right')
    ax.grid(True, alpha=0.3, which='both')
    
    # Inset: fractional enhancement
    ax_inset = inset_axes(ax, width="40%", height="35%", loc='lower left',
                          bbox_to_anchor=(0.12, 0.12, 1, 1), bbox_transform=ax.transAxes)
    ax_inset.plot(M_vals, delta_Lambda_pct, 'ko-', markersize=6, linewidth=1.5)
    ax_inset.axhline(0, color='gray', linestyle='--', linewidth=0.8)
    ax_inset.set_xlabel(r'$M$ [$M_\odot$]', fontsize=9)
    ax_inset.set_ylabel(r'$\Delta\Lambda/\Lambda$ [%]', fontsize=9)
    ax_inset.tick_params(labelsize=8)
    ax_inset.set_xlim(0.9, 1.7)
    ax_inset.set_ylim(-5, 65)
    ax_inset.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/fig1_lambda.pdf', bbox_inches='tight')
    plt.savefig(f'{output_dir}/fig1_lambda.png', bbox_inches='tight')
    plt.close()
    print("Created: fig1_lambda.pdf")


def fig2_gw170817_predictions(output_dir):
    """Figure 2: GW170817 prediction space."""
    
    fig, ax = plt.subplots(figsize=(7, 6))
    
    # Binary parameters
    M_tot = 2.72
    q_range = np.linspace(0.73, 1.0, 50)
    
    # Interpolate Lambda(M)
    Lambda_interp_GR = interp1d(M_vals, GR['Lambda'], kind='cubic', fill_value='extrapolate')
    Lambda_interp_RDT = interp1d(M_vals, RDT['Lambda'], kind='cubic', fill_value='extrapolate')
    
    def get_component_masses(M_tot, q):
        m1 = M_tot / (1 + q)
        m2 = q * m1
        return m1, m2
    
    # Generate prediction curves
    Lambda1_GR, Lambda2_GR = [], []
    Lambda1_RDT, Lambda2_RDT = [], []
    
    for q in q_range:
        m1, m2 = get_component_masses(M_tot, q)
        Lambda1_GR.append(float(Lambda_interp_GR(m1)))
        Lambda2_GR.append(float(Lambda_interp_GR(m2)))
        Lambda1_RDT.append(float(Lambda_interp_RDT(m1)))
        Lambda2_RDT.append(float(Lambda_interp_RDT(m2)))
    
    # GW170817 constraint region
    ax.fill_between([0, 800], [0, 0], [800, 800], alpha=0.15, color='green',
                    label='GW170817 90% CI')
    
    # Constant Lambda_tilde contours
    Lambda1_grid = np.linspace(100, 2500, 100)
    for Lt in [300, 500, 720, 1000, 1500]:
        Lambda2_contour = 2*Lt - Lambda1_grid
        valid = Lambda2_contour > 0
        ax.plot(Lambda1_grid[valid], Lambda2_contour[valid], 'gray',
                linestyle=':', alpha=0.5, linewidth=0.8)
        idx = len(Lambda1_grid)//3
        if Lambda2_contour[idx] > 0 and Lambda2_contour[idx] < 2500:
            ax.text(Lambda1_grid[idx], Lambda2_contour[idx], f'{Lt}',
                    fontsize=8, color='gray', ha='center', va='bottom')
    
    # Prediction curves
    ax.plot(Lambda1_GR, Lambda2_GR, 'b-', linewidth=2.5, label='GR (SLy4)')
    ax.plot(Lambda1_RDT, Lambda2_RDT, 'r-', linewidth=2.5, label=r'RDT ($\alpha=0.30$)')
    
    # Mark q values
    for q_mark in [0.75, 0.85, 0.95]:
        m1, m2 = get_component_masses(M_tot, q_mark)
        L1_gr = float(Lambda_interp_GR(m1))
        L2_gr = float(Lambda_interp_GR(m2))
        L1_rdt = float(Lambda_interp_RDT(m1))
        L2_rdt = float(Lambda_interp_RDT(m2))
        ax.plot(L1_gr, L2_gr, 'bo', markersize=8, markeredgecolor='black')
        ax.plot(L1_rdt, L2_rdt, 'r^', markersize=8, markeredgecolor='black')
        ax.annotate(f'q={q_mark}', xy=(L1_gr, L2_gr), xytext=(5, 5),
                    textcoords='offset points', fontsize=8, color='blue')
    
    # Equal mass line
    ax.plot([0, 2500], [0, 2500], 'k--', linewidth=1, alpha=0.5)
    
    ax.set_xlabel(r'$\Lambda_1$ (heavier component)')
    ax.set_ylabel(r'$\Lambda_2$ (lighter component)')
    ax.set_xlim(0, 2500)
    ax.set_ylim(0, 2500)
    ax.legend(loc='upper left')
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)
    
    ax.text(0.97, 0.03, f'$M_{{\\rm tot}} = {M_tot}\\,M_\\odot$',
            transform=ax.transAxes, ha='right', va='bottom', fontsize=10,
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/fig2_gw170817.pdf', bbox_inches='tight')
    plt.savefig(f'{output_dir}/fig2_gw170817.png', bbox_inches='tight')
    plt.close()
    print("Created: fig2_gw170817.pdf")


def fig3_ilove_relation(output_dir):
    """Figure 3: I-Love universal relation."""
    
    fig, ax = plt.subplots(figsize=(7, 5.5))
    
    # Yagi-Yunes universal relation
    YY = {'a': 1.47, 'b': 0.0817, 'c': 0.0149, 'd': 2.87e-4, 'e': -3.64e-5}
    
    def YY_lnI(ln_lambda):
        return (YY['a'] + YY['b']*ln_lambda + YY['c']*ln_lambda**2 +
                YY['d']*ln_lambda**3 + YY['e']*ln_lambda**4)
    
    ln_lambda_range = np.linspace(5, 10, 100)
    ln_I_YY = YY_lnI(ln_lambda_range)
    
    # Universal relation curve
    ax.plot(ln_lambda_range, ln_I_YY, 'k-', linewidth=2.5,
            label='Yagi & Yunes (2013)')
    
    # Data points
    ax.scatter(np.log(GR['Lambda']), np.log(GR['I_bar']),
               s=120, c='royalblue', marker='o', edgecolors='black',
               linewidths=1.2, label='GR (SLy4)', zorder=10)
    ax.scatter(np.log(RDT['Lambda']), np.log(RDT['I_bar']),
               s=120, c='crimson', marker='^', edgecolors='black',
               linewidths=1.2, label=r'RDT ($\alpha=0.30$)', zorder=10)
    
    # Mass labels
    for i, M in enumerate(M_vals):
        ax.annotate(f'{M:.1f}',
                    xy=(np.log(GR['Lambda'][i]), np.log(GR['I_bar'][i])),
                    xytext=(-8, -12), textcoords='offset points',
                    fontsize=9, color='royalblue', ha='center')
    
    ax.set_xlabel(r'$\ln(\bar{\lambda})$')
    ax.set_ylabel(r'$\ln(\bar{I})$')
    ax.legend(loc='lower right')
    ax.grid(True, alpha=0.3)
    ax.set_xlim(5.5, 9.5)
    ax.set_ylim(2.3, 3.5)
    
    # Inset: residuals
    ax_inset = inset_axes(ax, width="40%", height="35%", loc='upper left',
                          bbox_to_anchor=(0.08, 0.55, 1, 1), bbox_transform=ax.transAxes)
    
    def get_residual(Lambda, I_bar):
        ln_lam = np.log(Lambda)
        ln_I_pred = YY_lnI(ln_lam)
        return (np.log(I_bar) - ln_I_pred) * 100
    
    res_GR = get_residual(GR['Lambda'], GR['I_bar'])
    res_RDT = get_residual(RDT['Lambda'], RDT['I_bar'])
    
    ax_inset.axhline(0, color='gray', linestyle='--', linewidth=0.8)
    ax_inset.axhspan(-10, 10, alpha=0.15, color='gray')
    ax_inset.plot(M_vals, res_GR, 'bo-', markersize=5, linewidth=1)
    ax_inset.plot(M_vals, res_RDT, 'r^-', markersize=5, linewidth=1)
    ax_inset.set_xlabel(r'$M$ [$M_\odot$]', fontsize=9)
    ax_inset.set_ylabel('Residual [%]', fontsize=9)
    ax_inset.tick_params(labelsize=8)
    ax_inset.set_xlim(0.9, 1.7)
    ax_inset.set_ylim(-15, 15)
    ax_inset.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/fig3_ilove.pdf', bbox_inches='tight')
    plt.savefig(f'{output_dir}/fig3_ilove.png', bbox_inches='tight')
    plt.close()
    print("Created: fig3_ilove.pdf")


def fig4_scaling_verification(output_dir):
    """Figure 4: R^5 scaling verification."""
    
    fig, ax = plt.subplots(figsize=(6, 5))
    
    # Expected from R^5 scaling
    geometric_contribution = 5 * delta_R_pct
    k2_contribution = delta_k2_pct
    total_expected = geometric_contribution + k2_contribution
    
    # Plot
    ax.scatter(total_expected, delta_Lambda_pct, s=150, c='royalblue',
               marker='o', edgecolors='black', linewidths=1.5, zorder=10)
    
    # Mass labels
    for i, M in enumerate(M_vals):
        ax.annotate(f'{M:.1f} $M_\\odot$',
                    xy=(total_expected[i], delta_Lambda_pct[i]),
                    xytext=(8, -5), textcoords='offset points', fontsize=10)
    
    # Perfect correlation line
    x_line = np.linspace(0, 70, 100)
    ax.plot(x_line, x_line, 'k--', linewidth=2, label='Perfect scaling')
    
    # Fit line
    slope = np.sum(total_expected * delta_Lambda_pct) / np.sum(total_expected**2)
    ax.plot(x_line, slope * x_line, 'r-', linewidth=1.5, alpha=0.7,
            label=f'Fit: slope = {slope:.3f}')
    
    ax.set_xlabel(r'$\Delta k_2/k_2 + 5 \times \Delta R/R$ [%]')
    ax.set_ylabel(r'$\Delta\Lambda/\Lambda$ [%]')
    ax.legend(loc='upper left')
    ax.grid(True, alpha=0.3)
    ax.set_xlim(-5, 70)
    ax.set_ylim(-5, 70)
    ax.set_aspect('equal')
    
    # Formula box
    textstr = r'$\Lambda = \frac{2}{3}k_2\left(\frac{R}{M}\right)^5$'
    ax.text(0.95, 0.05, textstr, transform=ax.transAxes, fontsize=11,
            ha='right', va='bottom',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.9))
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/fig4_scaling.pdf', bbox_inches='tight')
    plt.savefig(f'{output_dir}/fig4_scaling.png', bbox_inches='tight')
    plt.close()
    print("Created: fig4_scaling.pdf")


def main():
    parser = argparse.ArgumentParser(description='Generate publication figures')
    parser.add_argument('--output', type=str, default='.',
                        help='Output directory for figures')
    args = parser.parse_args()
    
    os.makedirs(args.output, exist_ok=True)
    
    print("=" * 60)
    print("Generating Publication Figures")
    print("=" * 60)
    
    fig1_lambda_comparison(args.output)
    fig2_gw170817_predictions(args.output)
    fig3_ilove_relation(args.output)
    fig4_scaling_verification(args.output)
    
    print("\nAll figures generated successfully!")


if __name__ == "__main__":
    main()
