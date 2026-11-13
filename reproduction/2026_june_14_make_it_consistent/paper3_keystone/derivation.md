# Paper III keystone — relativistic divergence-operator modification

Goal: derive the GR (TOV-level) version of the **Papers I–II divergence-operator**
modification, integrate it at the canonical parameter point, and read off the
**sign** of the NS effect (ΔM_max, ΔR, ΔΛ). No tuning: the sign is whatever the
faithful GR analog of I–II produces. We test robustness with ≥2 inequivalent
embeddings; if they agree we trust the sign, if they disagree the GR embedding is
genuinely ambiguous (→ audit branch 1b).

## 1. The canonical Papers I–II content (what must be reproduced)
Modification lives in the gravitational **field/Gauss law**, not in matter:
- Mass integral stays 3D: `dm/dr = 4π r² ρ` (Newtonian) / `= 4π r² ε/c²` (GR).
- The radial gravitational acceleration `a(r)` felt by matter obeys a
  **(3+δ)-dimensional Gauss law**, δ = d_eff − 3 = 3(Ω−1) ≤ 0:
  ```
  (1/r^{2+δ}) d/dr( r^{2+δ} a ) = 4πGρ          (★ field form)
  ⇔  a' + [ (2+δ)/r + δ' ln(r/r⋆) ] a = 4πGρ
  ```
  with reference scale r⋆ = a_LE (so the log is ln ξ; unit-invariant).
- Hydrostatic balance: `dP/dr = −ρ a`.
- δ=0 ⇒ a' + 2a/r = 4πGρ ⇒ a = Gm/r² (standard). ✓

**Established result (computed last session, Papers I–II):** δ<0 in (★) →
the modified field reorganizes so the calibrated star is **smaller** (ξ₁
decreases for the Sun; white dwarfs show a fixed-mass radius **deficit**
|ΔR/R| ≤ 9.5e−3). So the canonical sign is **contraction**.

Contrast — the III/IV F-multiplier: `dP/dr → F·dP/dr`, F=(d_eff−1)/2=1+δ/2<1,
weakens the pressure gradient → **larger** star (expansion). Opposite sign.
This is the inconsistency to resolve.

## 2. Standard relativistic anchor (δ=0, must reproduce V7)
```
dP/dr = −(ε+P) · 𝒜_TOV(r),   𝒜_TOV = (m + 4πr³P)/(r(r−2m))     [G=c=1]
dm/dr = 4π r² ε
```
𝒜_TOV is the relativistic gravitational acceleration (= ½ ν'). Newtonian limit
𝒜_TOV → m/r² = a_N. V7 baseline: M_max=2.0495, R(1.4)=11.686.

## 3. The embeddings (the only ambiguity is here)
The invariant content (★) fixes the **weak-field** acceleration. The relativistic
promotion — how the pressure source 4πr³P and the metric factor (r−2m) ride along
with the r^{2+δ} operator — is not unique. Two defensible, inequivalent forms,
each satisfying both limits (δ=0 → TOV; weak field → ★):

### Embedding 1 — operator correction factor (primary)
Solve (★) for `a_op(r)` and its δ=0 counterpart `a_N(r)` on the *same*
relativistic density profile ρ(r); define the dimensionless ratio
`R_op(r) = a_op(r)/a_N(r)` (the pure operator effect on the field). Apply it
multiplicatively to the relativistic acceleration:
```
dP/dr = −(ε+P) · 𝒜_TOV(r) · R_op(r)
```
- δ=0 ⇒ R_op=1 ⇒ exact TOV. ✓
- weak field ⇒ 𝒜_TOV→a_N ⇒ dP/dr → −ρ a_N R_op = −ρ a_op (= ★). ✓
- Isolates the operator sign as R_op vs 1. Contraction ⟺ R_op>1 (stronger field).

### Embedding 2 — effective-dimension mass (alternative)
Promote the geometric measure in the acceleration directly, paralleling the
Newtonian closed form `a = 4πG M_{3+δ}/r^{2+δ}`, `M_{3+δ}=∫ρ r'^{2+δ}dr'`:
```
𝒜_mod = ( m_δ + 4π r³ (r/r⋆)^δ P ) / ( r²(r/r⋆)^δ · (1 − 2m/r) )
m_δ(r) = ∫₀^r 4π r'^2 (r'/r⋆)^δ ε dr'     (3D mass with the (r/r⋆)^δ operator weight)
dP/dr = −(ε+P) 𝒜_mod
```
- δ=0 ⇒ (r/r⋆)^δ=1 ⇒ m_δ=m, 𝒜_mod=𝒜_TOV. ✓
- weak field ⇒ reduces to the integral form of (★). ✓
- Genuinely different from Emb. 1: modifies via a non-local weighted mass + the
  geometric falloff, not a local field ratio.

(Reference-scale handling `r^{2+δ} = r²(r/r⋆)^δ`, r⋆=a_LE, closes the units gap
flagged in `paper6_units_note.md`. Sensitivity to r⋆=R vs a_LE is a secondary
robustness check, as in Paper I.)

## 4. Gates (refuse to report a sign unless all pass)
1. **δ=0 → V7**: both embeddings reproduce M_max=2.0495, R(1.4)=11.686 to 0.1%.
2. **Weak-field → I–II**: at low compactness the fractional ΔR matches the
   Newtonian operator result (sign and rough magnitude of the WD deficit).
3. **Tolerance/units invariance**: stable to 10× tighter integrator tol; cgs/SI
   (or geometric) invariant; r⋆ sensitivity bounded.

## 5. Sign verdict logic (→ SIGN_VERDICT.md, committed before any number rewrite)
At the canonical A₉₅(α) boundary points (same as III V7; NOT α=1):
- Both embeddings **contract** (smaller M_max, smaller R, smaller Λ) → audit
  branch **1a**: the operator form is the canonical relativistic mechanism;
  rewrite III's numbers with the contraction sign; IV re-runs downstream.
- Embeddings **disagree in sign** → the GR embedding is genuinely ambiguous →
  audit branch **1b**: F-multiplier stays as an explicitly labeled non-canonical
  ansatz, +0.9% demoted to "ansatz illustration", one-paragraph open problem.

**Prior (to be tested, not assumed):** contraction (branch 1a). The Newtonian
operator already contracts; relativistic corrections amplify gravity without
flipping the source→field relation that sets the sign. Magnitude expected
comparable to the WD/NS bound scale (sub-percent to ~1%).

## 6. Build plan
- `tov_operator.py`: TOV with δ=0 anchor + Embeddings 1 & 2, reusing ns_v7
  eos.py (EOS) and surface detection. Single-thread.
- `validate_operator.py`: the three gates.
- `mr_operator.py`: M–R sweep + M_max + R(1.4) for GR / Emb1 / Emb2.
- Read sign → SIGN_VERDICT.md → check in with Chris before rewriting III.
