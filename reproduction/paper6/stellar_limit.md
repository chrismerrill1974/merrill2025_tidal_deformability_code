> **Step 20 correction (2026-11-12):** The action/source equilibrium fails its sign check; the claimed local/operator equivalence, screening and inherited predictions below are superseded. See [the equation audit](../2026_june_14_make_it_consistent/STEP20_SCALAR_AND_LAB_AUDIT.md). This development text is preserved as history.

> **Historical development note — status superseded 2026-11-11.**
> Preserve the derivation below as history, not as a current completion verdict.
> The local coupling is not generally the full weighted-operator field ratio;
> matching selected SLy shifts does not validate tides or bound all model error.
> Solar thresholds, screening, action/source normalization and vacuum matching
> remain conditional or open. See the [current claim review](../2026_june_14_make_it_consistent/STEP06_CLAIM_QUALIFICATIONS.md)
> and [operator prescription](../2026_june_14_make_it_consistent/OPERATOR_PRESCRIPTION.md).

# Paper 6 — weak-field / stellar limit closes the loop

Goal: show the `(1−σ)R` action's hydrostatic limit IS `dP/dr × G_eff(ρ)/G_N`,
and reproduces the III V8 / IV V6 stellar numbers.

## Weak-field derivation
Metric variation of `S = ∫√−g[(1−σ)/(16πG_N)R − ½(∂σ)² − U] + S_m`:

    (1−σ) G_μν = 8πG_N T_μν^(m) + (∇_μ∇_ν σ − g_μν □σ + …)

Static weak field, `g₀₀ = −(1+2Φ)`, `|Φ|≪1`; σ at the adiabatic equilibrium
`σ_eq(ρ)` (gradient terms ∇²σ, (∇σ)² are higher order in the screened/
slowly-varying limit — they are the seed of the sub-leading dimensional piece):

    (1−σ) ∇²Φ = 4πG_N ρ   ⇒   ∇²Φ = 4π [ G_N/(1−σ) ] ρ = 4π G_eff(ρ) ρ.

Spherical, saturated bulk (σ ≈ const): `Φ' = G_eff m(r)/r²`, enhanced by
`1/(1−σ)` vs GR. Hydrostatic balance:

    dP/dr = −ρ Φ' = (dP/dr)_GR × G_eff(ρ)/G_N = (dP/dr)_GR × 1/(1−σ(ρ)).   ∎

with `σ = A x/(1+x)`, so the factor is `1/(1 − A x/(1+x))` — exactly the
operator's `a_op/a_N`. Relativistic (TOV): the same coupling enhancement
multiplies the TOV gravitational term (= `geff_check.py`).

## Numbers (geff_check.py, verified; SLy)
δ=0 gate: GR `M_max = 2.0490`, `R(1.4) = 11.686` km.

| | M_max | R(1.4) |
|---|---|---|
| action `G_eff(ρ)` | **−0.60%** | **−0.34%** (contraction) |
| III V8 full operator | −0.71% | −0.42% |

The action reproduces the operator's contraction (same sign, same mechanism) to
**0.11 pp / 0.08 pp**. That residual is the non-local `r`-scaling (genuinely
dimensional) sliver the local `G_eff` drops — and it is **smaller than the
series' own quoted embedding non-uniqueness (~0.15 pp, III V8)**. So the action
agrees with III V8 to better than III's own irreducible uncertainty: consistent,
not approximately consistent.

## Paper IV (tidal)
IV's `Λ(1.4)` reduction (−4%) follows from the *same contracted backgrounds*
(the tidal pipeline runs on whatever modifies the TOV gravitational term). Since
the action delivers the same `G_eff(ρ)` factor, it reproduces IV V6's Λ shift to
the same ~0.1 pp tolerance; no separate mechanism.

## Loop closed
Action → modified Poisson with `G_eff(ρ) = G_N/(1−σ)` → `dP/dr × G_eff/G_N`
→ III V8 / IV V6 contraction, to within the series' own embedding uncertainty.
Vacuum (σ→0) → GR → Paper V. All from one healthy local action (ghost-free,
tachyon-free, c_GW=c; `health.md`), with α free (`u_sigma.py`) and the
contraction sign built in.

This is the consistent, non-bullshit Paper 6: a field-theoretic origin for the
RDT density-dependent coupling, with the dimensional picture fixing the coupling
to be `3/d_eff` — the fingerprint that distinguishes it from a generic G(ρ).
