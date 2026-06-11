# Tidal Deformability Code — Radial Dimensionality Theory Paper IV (V5)

Code for **Radial Dimensionality Theory IV: Tidal Deformability Ceilings
under Solar-Constrained Density-Dependent Dimensionality**
(C. K. Merrill and Claude, 2026). Supersedes the V4 code, now in
[`archive_v4/`](archive_v4/).

Relativistic tidal deformability (Hinderer 2008 / Postnikov et al. 2010)
on TOV backgrounds validated to ~0.1% (shared with the Paper III repo,
`merrill2025_neutronstar_code`). RDT enters as the series saturating law
bounded by Paper I's solar constraint (A < 0.0087, 95% CL). Headline
ceiling: ΔΛ(1.4)/Λ ≤ +5.0% (Λ = 298 → 313), split ~evenly between R⁵
scaling and the Love number; GW170817: Λ̃ = 359 (GR) → ≤ 377 (RDT), both
inside the measured band — reversing the V4 conclusion.

## Requirements

Python ≥ 3.10 with `numpy`, `scipy`, `matplotlib`.

## Files

| File | Purpose |
|---|---|
| `eos.py`, `tov.py` | SLy EOS (CompOSE crust-to-core table, ε taken directly) and TOV solver — same modules as the Paper III repo. |
| `tidal.py` | Hinderer y-equation integrated alongside TOV (minimal coupling: GR perturbation equations + standard Einstein ν′ on the F-modified background), Love number k₂, Λ, binary Λ̃; n=1 polytrope class for the benchmark gate. |
| `validate.py` | Gates. **Run first.** (1) n=1 polytrope Newtonian limit k₂(C→0) = 0.2594 vs analytic 0.2599 — note: V4's claimed benchmark "k₂ = 0.26 at C = 0.15" was this Newtonian value mislabeled (relativistic value ≈ 0.071); (2) background consistency vs tov.py; (3) published-SLy regression Λ(1.4) = 298.0 (literature ≈ 297) — the gate V4 failed by ~3×; (4) tolerance stability. |
| `predict.py` | Λ(M), k₂/R⁵ decomposition, and GW170817 Λ̃ along Paper I's exclusion boundary A₉₅(α). Prints the paper-numbers block; writes `results_tidal.csv`. |
| `figures.py` | Regenerates the paper figure (`figures/fig_lambda_v5.pdf`). |
| `data/` | CompOSE RG(SLY4) table (`eos.thermo`, `eos.nb`, reference curve `eos.mr`) and the Douchin & Haensel table. |
| `archive_v4/` | Superseded V4 code, data, figures, and manuscript, kept for the record (its EOS module carries the energy-density reconstruction bug documented in Papers III–IV). |

## Usage

```sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
       VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1

python3 validate.py   # all 4 gates must PASS first
python3 predict.py    # paper-numbers block + results_tidal.csv (~2 min)
python3 figures.py    # paper figure (~1 min)
```

## License

MIT — see [LICENSE](LICENSE).
