# Tidal Deformability Code — Radial Dimensionality Theory Paper IV (V6)

Code for **Radial Dimensionality Theory IV: Tidal Deformability Ceilings
under Solar-Constrained Density-Dependent Dimensionality**
(C. K. Merrill and Claude, 2026). V6 supersedes V5 (its force-rescaling
scripts are archived in [`archive_v5/`](archive_v5/); the earlier V4 code is
in [`archive_v4/`](archive_v4/)).

**What changed in V6 — the sign.** V5 computed Λ on backgrounds modified by
a force-rescaling factor F(ρ) = (d_eff − 1)/2 multiplying the TOV pressure
gradient, which weakens gravity and *expands* the star, giving a tidal
**enhancement** (ΔΛ(1.4) = +5%, Λ = 298 → 313). That is **not** the
mechanism of Papers I–II, which modify the gravitational **divergence
operator**. The faithful relativistic analog of that operator — the same
mechanism adopted in Paper III V8 (`merrill2025_neutronstar_code`) —
**contracts** the star, so the tidal deformability is **reduced**:

- ΔΛ(1.4) = **−4.0%** (Λ = 298 → 286), envelope −3.7% (1.2 M⊙) to −4.5%
  (1.6 M⊙), bounded by Paper I's solar constraint A < 0.0087 (95% CL).
- Decomposition: 5 ΔR/R ≈ −2.1% plus Δk₂/k₂ ≈ −1.9% — the contracted
  density profile lowers k₂ as well as R, so the full Love-number result is
  ~2× the R⁵-only estimate (decomposition closes,
  ΔΛ/Λ ≈ 5 ΔR/R + Δk₂/k₂).
- GW170817: Λ̃ = 359 (GR) → 345 (RDT ceiling), both inside the measured
  Λ̃ = 300⁺⁴²⁰₋₂₃₀ — the contraction nudges Λ̃ toward the lower edge, not
  the upper (V5 gave 377).

The shift is **computed**, not scaled: A = 0 reduces the operator pipeline
to GR (Λ(1.4) = 298, literature ≈ 297).

## Requirements

Python ≥ 3.10 with `numpy`, `scipy`, `matplotlib`.

## Files (V6 — operator form)

| File | Purpose |
|---|---|
| `eos.py`, `tov.py` | SLy EOS (CompOSE crust-to-core table, ε taken directly) and baseline TOV solver/constants — the same modules as the Paper III repo. (`tov.py` also still contains the **superseded** V5 F-multiplier `F_rdt`; V6 does not use it.) |
| `tidal.py` | Hinderer y-equation Love-number pipeline: `love_k2`, `cs2_of_P`, Λ, binary Λ̃, and the n=1 polytrope class for the benchmark gate. Its `solve_tidal` still carries the V5 F-multiplier path (exercised by `validate.py` only at A = 0, i.e. GR); the V6 operator form is in `operator_tidal.py`. |
| `tov_operator.py` | The divergence-operator modified Gauss law (δ = d_eff − 3 ≤ 0): `delta`, `ddelta_drho`, `_drho_dP`. Shared with the Paper III repo. |
| `operator_tidal.py` | **The V6 mechanism.** Hinderer integration on the operator (Embedding 1) backgrounds: the field ratio R_op = a_op/a_N multiplies **both** dP/dr and the metric potential ν′ — which is why k₂ moves, not the radius alone. A = 0 → GR (built-in gate: Λ(1.4) = 298). Prints Λ(M) and the k₂-vs-R decomposition at the global solar ceiling (A = 0.0087, α = 0.122). |
| `validate.py` | GR / Newtonian gates. **Run first.** (1) n=1 polytrope Newtonian limit k₂(C→0) = 0.2594 vs analytic 0.2599 — note: V4's claimed benchmark "k₂ = 0.26 at C = 0.15" was this Newtonian value mislabeled (relativistic value ≈ 0.071); (2) background consistency vs `tov.py`; (3) published-SLy regression Λ(1.4) = 298.0 (literature ≈ 297) — the gate V4 failed by ~3×; (4) tolerance stability. |
| `figures_v6.py` | Regenerates the paper figure (`figures/fig_lambda_v6.pdf`): Λ(M) with the RDT curve *below* GR (reduction), and the fixed-mass reduction envelope along Paper I's A₉₅(α) boundary. |
| `data/` | CompOSE RG(SLY4) table (`eos.thermo`, `eos.nb`, reference curve `eos.mr`) and the Douchin & Haensel table. |
| `archive_v5/` | Superseded V5 force-rescaling scripts (`predict.py`, `figures.py`, `results_tidal.csv`) — the F-multiplier that carried the wrong sign (+5% enhancement). Retained for the record. |
| `archive_v4/` | Superseded V4 code, data, figures, and manuscript (its EOS module carries the energy-density reconstruction bug documented in Papers III–IV). |

## Usage

```sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
       VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1

python3 validate.py        # GR / Newtonian gates — run first (all 4 must PASS)
python3 operator_tidal.py  # V6 operator form: Lambda(M) + k2-vs-R decomposition at the ceiling
python3 figures_v6.py      # paper figure -> figures/fig_lambda_v6.pdf
```

The operator backgrounds (`tov_operator.py`) and the EOS/TOV base
(`eos.py`, `tov.py`) are shared with the Paper III repo
(`merrill2025_neutronstar_code`); this repo adds the relativistic tidal
(Love-number) layer on top of them.

## Note on the sign reversal

V6 reverses a published headline: V5 reported RDT *raises* Λ (a +5%
enhancement, Λ̃ = 377 for GW170817); the canonical operator form shows RDT
*lowers* Λ (−4.0%, Λ̃ = 345). This is the downstream tidal counterpart of
Paper III's V8 sign reversal — III shows the operator form contracts the
star (ΔM_max = −0.71%, ΔR(1.4) = −0.42%), and the contracted density
profile both shrinks R and lowers k₂, so Λ falls by ~2× the R⁵-only
estimate. The reduction matches the radius deficit of Papers I–II
(everything contracts). The remaining theoretical freedom is the
~0.15-percentage-point non-uniqueness of the relativistic embedding (a
magnitude, no longer a sign, ambiguity), documented in the Paper III repo's
`SIGN_VERDICT.md`.

## License

MIT — see [LICENSE](LICENSE).
