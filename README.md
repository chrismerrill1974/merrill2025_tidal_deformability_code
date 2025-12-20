# RDT Paper IV: Tidal Deformability of Neutron Stars

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

This repository contains the computational code and data for **Paper IV** of the Recursive Dimensionality Theory (RDT) series:

> **"Tidal Deformability of Neutron Stars in Recursive Dimensionality Theory: Gravitational Wave Signatures and Universal Relations"**
> 
> Christopher Merrill, with computational collaboration from Claude (Anthropic)

## Overview

We compute tidal deformabilities of neutron stars within Recursive Dimensionality Theory (RDT), a framework where effective spatial dimensionality varies with matter density:

$$d_{\rm eff}(\rho) = 3 + \alpha \log_{10}\left(\frac{\rho}{\rho_0}\right), \quad \rho > \rho_0$$

where $\rho_0 = 2.7 \times 10^{14}$ g/cm³ is nuclear saturation density and $\alpha$ is the dimensionless coupling constant.

### Key Results

| Mass (M☉) | Λ_GR | Λ_RDT | Enhancement |
|-----------|------|-------|-------------|
| 1.0 | 10241 | 10483 | +2.4% |
| 1.2 | 2995 | 3274 | +9.3% |
| 1.4 | 985 | 1221 | **+24.0%** |
| 1.6 | 332 | 516 | +55.8% |

- **GW170817 predictions**: Λ̃_GR = 1276, Λ̃_RDT = 1519 (19% enhancement)
- **I-Love universality**: Preserved in RDT
- **EOS independence**: Fractional shifts identical for SLy4 and APR4

## Repository Structure

```
rdt_paper4_release/
├── src/
│   ├── tov_solver.py          # TOV equation solver (GR and RDT)
│   ├── tidal_solver.py        # Tidal perturbation equations
│   ├── eos_sly4.py            # SLy4 equation of state
│   ├── love_number.py         # k₂ and Λ extraction
│   ├── mass_sequence.py       # Generate M-R-Λ sequences
│   ├── gw170817_analysis.py   # Binary tidal deformability
│   ├── ilove_relations.py     # I-Love universal relations
│   └── generate_figures.py    # Publication figure generation
├── data/
│   ├── profiles/              # Stellar structure profiles
│   ├── tables/                # Publication tables (LaTeX)
│   ├── gr_lambda_vs_mass.dat  # GR tidal sequence
│   ├── rdt_lambda_vs_mass.dat # RDT tidal sequence
│   └── ilove_data.dat         # I-Love relation data
├── figures/
│   ├── fig1_lambda.pdf        # Λ(M) comparison
│   ├── fig2_gw170817.pdf      # GW170817 prediction space
│   ├── fig3_ilove.pdf         # I-Love universal relation
│   └── fig4_scaling.pdf       # R⁵ geometric scaling
├── manuscript/
│   └── paper4_complete.tex    # Full manuscript (RevTeX4)
├── README.md
├── requirements.txt
└── LICENSE
```

## Installation

```bash
git clone https://github.com/USERNAME/rdt-paper4-tidal.git
cd rdt-paper4-tidal
pip install -r requirements.txt
```

### Requirements

- Python 3.8+
- NumPy ≥ 1.20
- SciPy ≥ 1.7
- Matplotlib ≥ 3.5

## Usage

### Basic TOV + Tidal Calculation

```python
from src.tov_solver import solve_tov_gr, solve_tov_rdt
from src.tidal_solver import solve_tidal
from src.love_number import compute_k2, compute_Lambda

# Solve TOV for 1.4 solar mass star
M_target = 1.4  # solar masses
alpha = 0.30    # RDT coupling

# GR solution
r_gr, m_gr, P_gr, rho_gr = solve_tov_gr(M_target)
y_R_gr = solve_tidal(r_gr, m_gr, P_gr, rho_gr)
k2_gr = compute_k2(m_gr[-1], r_gr[-1], y_R_gr)
Lambda_gr = compute_Lambda(k2_gr, m_gr[-1], r_gr[-1])

# RDT solution
r_rdt, m_rdt, P_rdt, rho_rdt = solve_tov_rdt(M_target, alpha)
y_R_rdt = solve_tidal(r_rdt, m_rdt, P_rdt, rho_rdt)
k2_rdt = compute_k2(m_rdt[-1], r_rdt[-1], y_R_rdt)
Lambda_rdt = compute_Lambda(k2_rdt, m_rdt[-1], r_rdt[-1])

print(f"GR:  Λ = {Lambda_gr:.1f}")
print(f"RDT: Λ = {Lambda_rdt:.1f} ({100*(Lambda_rdt/Lambda_gr-1):+.1f}%)")
```

### Generate Full Mass Sequence

```bash
python src/mass_sequence.py --masses 1.0 1.2 1.4 1.6 --alpha 0.30 --output data/
```

### Reproduce Publication Figures

```bash
python src/generate_figures.py --output figures/
```

## Physical Background

### RDT-Modified TOV Equations

The standard TOV pressure equation is modified by a geometric factor:

$$\frac{dP}{dr} = -\frac{(\epsilon + P)(m + 4\pi r^3 P)}{r(r - 2m)} \times F(\rho)$$

where:

$$F(\rho) = \frac{d_{\rm eff}(\rho) - 1}{2} = \begin{cases} 1 & \rho \leq \rho_0 \\ 1 + \frac{\alpha}{2}\log_{10}(\rho/\rho_0) & \rho > \rho_0 \end{cases}$$

### Tidal Perturbation Equations

We solve the Hinderer (2008) $y(r)$ equation on the RDT-modified background:

$$r \frac{dy}{dr} + y^2 + yF_1(r) + F_2(r) = 0$$

with boundary condition $y(0) = 2$.

### Love Number Extraction

The tidal Love number $k_2$ is computed from:

$$k_2 = \frac{8C^5}{5}(1-2C)^2[2+2C(y_R-1)-y_R] \times [\cdots]^{-1}$$

and the dimensionless tidal deformability:

$$\Lambda = \frac{2}{3} k_2 \left(\frac{R}{M}\right)^5$$

## Validation

Our implementation is validated against:

1. **Hinderer (2008)** polytrope benchmarks: $k_2 = 0.260$ for n=1 polytrope ✓
2. **Published SLy4 values**: $\Lambda_{1.4} \approx 985$ ✓
3. **Geometric scaling**: $\Delta\Lambda/\Lambda \approx \Delta k_2/k_2 + 5\Delta R/R$ ✓

## Citation

If you use this code in your research, please cite:

```bibtex
@article{Merrill2024_RDT_IV,
  author  = {Merrill, Christopher},
  title   = {Tidal Deformability of Neutron Stars in Recursive 
             Dimensionality Theory: Gravitational Wave Signatures 
             and Universal Relations},
  journal = {Zenodo},
  year    = {2024},
  doi     = {10.5281/zenodo.XXXXXXX}
}
```

Also cite the foundational papers:

- **Paper I** (Solar neutrinos): [Zenodo link]
- **Paper II** (White dwarfs): [Zenodo link]  
- **Paper III** (Neutron star structure): [Zenodo link]

## Related Publications

- Hinderer, T. (2008). "Tidal Love numbers of neutron stars." ApJ 677, 1216.
- Yagi, K. & Yunes, N. (2013). "I-Love-Q relations in neutron stars." Science 341, 365.
- Abbott, B.P. et al. (2017). "GW170817: Observation of gravitational waves from a binary neutron star inspiral." PRL 119, 161101.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- Claude (Anthropic) for computational collaboration
- The LIGO-Virgo-KAGRA Collaboration for GW170817 data
- The nuclear physics community for equation of state tabulations
