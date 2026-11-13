> **Step 20 correction (2026-11-12):** The action/source equilibrium fails its sign check; the claimed local/operator equivalence, screening and inherited predictions below are superseded. See [the equation audit](../2026_june_14_make_it_consistent/STEP20_SCALAR_AND_LAB_AUDIT.md). This development text is preserved as history.

> **Historical development note — status superseded 2026-11-11.**
> Preserve the derivation below as history, not as a current completion verdict.
> The local coupling is not generally the full weighted-operator field ratio;
> matching selected SLy shifts does not validate tides or bound all model error.
> Solar thresholds, screening, action/source normalization and vacuum matching
> remain conditional or open. See the [current claim review](../2026_june_14_make_it_consistent/STEP06_CLAIM_QUALIFICATIONS.md)
> and [operator prescription](../2026_june_14_make_it_consistent/OPERATOR_PRESCRIPTION.md).

# Paper 6 — the G_eff(ρ) action for RDT

**Decision (Chris):** build RDT's completion as a density-dependent effective
gravitational coupling. The *dimensional picture* is what fixes the form of the
coupling — it is not decoration.

## The physical input (the I–V claim)
At density ρ, matter behaves as if in slightly fewer than 3 spatial dimensions:
`d_eff(ρ) = 3 − ε(ρ)`, `ε(ρ) = 3A·x/(1+x)`, `x = (ρ/ρ₀)^α`, bounded `A < 0.0087`.
Gravitational flux that spreads into `d_eff < 3` dimensions is *concentrated*, so
the field is stronger. The central field ratio of a `d_eff`-dimensional Gauss law
is `3/d_eff`. **This dictates the coupling — it is not a free choice:**

    G_eff(ρ) = G · 3/d_eff(ρ) = G / (1 − A·x/(1+x))   > G   ⇒ contraction.

(Numerically verified to reproduce the operator's contraction; `derivation.md`.)
This is the *dimensional reduction* enhancement; a generic G(ρ) modification
would not have the `3/d_eff` form. That form is RDT's fingerprint.

## The field
Introduce the **dimension-deficit field** `σ(x)`, with `d_eff = 3(1 − σ)`, so
`σ = ε/3 = A·x/(1+x) ∈ [0, A)`. Gravity carries the dimensional enhancement
through a field-dependent coupling:

    S = ∫ d⁴x √−g [ (1 − σ)/(16π G_N) · R  −  ½(∂σ)²  −  U(σ) ]  +  S_m[ψ, g_μν]

- Gravitational sector: `1/(16π G_eff) = (1 − σ)/(16π G_N)` ⇒ `G_eff = G_N/(1 − σ)`.
  σ > 0 (dense) ⇒ stronger gravity ⇒ contraction. σ → 0 (vacuum) ⇒ `G_eff → G_N`,
  exactly GR — **Paper V falls out by construction.**
- Matter is minimally coupled to `g_μν` (no separate fifth-force coupling); σ is
  sourced by curvature/density through the non-minimal `(1−σ)R` term.

## σ tracks the density (the saturating law)
In the adiabatic / heavy-field (short-range) limit, σ sits at the minimum of an
effective potential `U_eff(σ;ρ) = U(σ) − σ·J(ρ)`, with the source `J(ρ) ∝ ρ` from
the `(1−σ)R` coupling (Newtonian trace `R ∼ 8πG ρ`). Requiring the equilibrium
`σ_eq(ρ) = A·x/(1+x)` and inverting,

    x = σ/(A − σ)  ⇒  ρ(σ) = ρ₀ · [ σ/(A − σ) ]^{1/α}  ⇒  U′(σ) = J(ρ(σ)) ∝ [σ/(A−σ)]^{1/α}.

**This is the key point: α enters only as the power `1/α` in `U(σ)`** — it is a
*free shape parameter of the potential*, recovered anywhere in `[0.02, 0.5]`. It
is **not forced to α = 1**. (The companion's quadratic-potential chameleon forced
`φ_eq ∝ ρ`, i.e. `α = 1` — that rigidity is gone because here the saturating form
is built into `U(σ)` directly.)

## Status of the bar
- **Contraction sign:** ✓ `G_eff > G` for σ > 0 (the F-multiplier / companion sign
  error is structurally avoided — the `3/d_eff` form is > 1 by construction).
- **Vacuum → GR:** ✓ σ → 0, `G_eff → G_N`.
- **α free in [0.02, 0.5]:** ✓ set by `U(σ)` shape, not forced.
- **Solar bound, downstream ceilings:** unchanged — same `G_eff(ρ)` the series
  already computes with.
- **`c_GW = c`, ghost/tachyon-free:** expected (non-minimal `f(σ)R` scalar-tensor
  with no curvature-derivative coupling ⇒ tensor speed c; Brans-Dicke-class
  kinetic structure) — **TO VERIFY explicitly** (Einstein-frame kinetic sign).

## Honest framing
"Effective dimension" is the *organizing physical picture* — it fixes the coupling
to be `3/d_eff`, which is what makes this RDT and not a generic G(ρ). The
operational, observable content is the density-dependent coupling `G_eff(ρ)`. In
RDT's regime (dense matter; vacuum = GR) the two are indistinguishable; the
sub-leading literal-dimension (force-falloff) piece is the negligible ~0.1 pp of
`derivation.md`.

## U(σ) — DONE (`u_sigma.py`, verified)
Inverting `U'(σ_eq)=J(ρ)∝ρ` with target `σ_eq=A x/(1+x)`:

    U'(σ) = λ [σ/(A−σ)]^{1/α},   U(σ) = λ ∫₀^σ [s/(A−s)]^{1/α} ds.

Numerically `σ(potential) == A x/(1+x)` to machine precision for α = 0.12, 0.30,
0.50, 1.0 (outer crust → NS core). Convex (U''>0 ⇒ tachyon-free in range). Steep
wall as σ→A ⇒ floor `d_eff>2.974` automatic. Solar anchor: at ρ=ρ₀=150 g/cm³,
σ=A/2, d_eff=2.987 (α-independent, since x=1). **α enters only as 1/α ⇒ FREE in
[0.02,0.5]; α=1 is just the companion's special case.** Drift fixed.

## Next
1. [x] U(σ) reproduces the saturating law; α free. (done)
2. [ ] Einstein-frame transform: confirm ghost-/tachyon-free and `c_GW = c`.
3. [ ] Show the weak-field stellar limit returns `dP/dr × G_eff(ρ)/G` (→ III/IV
       numbers) and vacuum → Schwarzschild/Kerr (→ Paper V).
