> **Step 20 correction (2026-11-12):** The action/source equilibrium fails its sign check; the claimed local/operator equivalence, screening and inherited predictions below are superseded. See [the equation audit](../2026_june_14_make_it_consistent/STEP20_SCALAR_AND_LAB_AUDIT.md). This development text is preserved as history.

> **Historical development note — status superseded 2026-11-11.**
> Preserve the derivation below as history, not as a current completion verdict.
> The local coupling is not generally the full weighted-operator field ratio;
> matching selected SLy shifts does not validate tides or bound all model error.
> Solar thresholds, screening, action/source normalization and vacuum matching
> remain conditional or open. See the [current claim review](../2026_june_14_make_it_consistent/STEP06_CLAIM_QUALIFICATIONS.md)
> and [operator prescription](../2026_june_14_make_it_consistent/OPERATOR_PRESCRIPTION.md).

# Paper 6 — Einstein-frame health of the G_eff(σ) action

Jordan-frame action:

    S = ∫ d⁴x √−g [ (1−σ)/(16πG_N) · R − ½(∂σ)² − U(σ) ] + S_m[ψ, g_μν]

A non-minimally-coupled f(σ)R scalar-tensor theory: coefficient of R/2 is
`F(σ) = (1−σ)/(8πG_N)`, Jordan kinetic coefficient `Z = 1` (correct sign),
matter minimally coupled to g_μν.

## Conformal transform to the Einstein frame
`g̃_μν = Ω² g_μν`, `Ω² = 8πG_N F = (1−σ) > 0` (since σ < A ≪ 1). Standard
scalar-tensor result (Fujii–Maeda / Faraoni):

    S = ∫ d⁴x √−g̃ [ R̃/(16πG_N) − ½ K(σ)(∂̃σ)² − U(σ)/(1−σ)² ] + S_m[ψ, (1−σ)⁻¹ g̃]

    K(σ) = Z/F + (3/2)(F'/F)²  ∝  1/(1−σ) + (3/2)/(1−σ)².

## Results
1. **Ghost-free.**  K = [Z/F > 0] + [(3/2)(F'/F)² ≥ 0] > 0 for every `F = (1−σ) > 0`,
   i.e. every σ < 1. Here σ ∈ [0, A), A < 0.0087 ⇒ F > 0.991, so **robustly no
   ghost**. (Even safer than f(R) gravity, which has Z = 0 and is still healthy
   via the (3/2)(F'/F)² term.)
2. **Tachyon-free / screening.**  σ does not sit at a bare U-minimum; it sits at
   the density-driven equilibrium σ_eq(ρ) of U_eff = U − σ J(ρ). The effective
   mass² = U''(σ_eq) > 0 (U convex — verified, `u_sigma.py`) and **grows toward
   the σ→A wall** ⇒ heavier / more screened at higher density: the chameleon
   mechanism. No tachyon.
3. **c_GW = c.**  Only an f(σ)R curvature coupling — no X-dependent (derivative)
   curvature term — so tensor modes are canonical. Manifest in the Einstein
   frame: gravity is standard Einstein–Hilbert `R̃/(16πG_N)` ⇒ gravitons at c.
   (And Ω² = 1−σ differs from 1 by < 0.9%, so the two frames nearly coincide.)
4. **Vacuum → GR.**  σ → 0 ⇒ F → 1/(8πG_N), Ω² → 1, U_E → U(0): exactly GR.
   Paper V holds by construction.

## Caveat (stated honestly)
- **Screening is required (not optional).** The Einstein-frame matter–scalar
  coupling is `β = ½ d ln(1−σ)⁻¹/dσ ≈ ½` — gravitational strength. An UNSCREENED
  light scalar with β ~ ½ violates Cassini (needs Brans–Dicke ω ≳ 4×10⁴). So the
  scalar must be heavy/short-range in the Solar System — the "stellar-relevant"
  screened branch the companion already identified — set by the σ kinetic
  normalization (an implicit mass scale M, free). A constraint, not a pathology;
  consistent with the series' solar bound.
- σ is dimensionless, so −½(∂σ)² carries an implicit scale M (i.e. −½M²(∂σ)²);
  M sets the scalar mass/range. It does not affect the ghost/tachyon signs.

## Verdict
**Healthy:** ghost-free, tachyon-free, c_GW = c, vacuum → GR — given F = 1−σ > 0
(always true) and the screened (heavy-scalar) regime. The theory is barely-
modified GR with a density-dependent coupling G_eff = G_N/(1−σ) = 3G_N/d_eff.
