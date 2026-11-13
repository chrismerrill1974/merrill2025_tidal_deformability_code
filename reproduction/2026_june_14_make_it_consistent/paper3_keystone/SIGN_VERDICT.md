# SIGN_VERDICT — Paper III keystone (committed before any number rewrite)

**Verdict: audit branch 1a. The divergence-operator form CONTRACTS.**
RDT makes neutron stars SMALLER and LOWERS the maximum mass — matching the
Papers I–II standard, and OPPOSITE to the III/IV F-multiplier (expansion).

## Numbers (canonical global ceiling A=0.0087, α=0.12; CompOSE SLy; G=c baseline V7)
| model | M_max [M⊙] | ΔM_max | R(1.4) [km] | ΔR(1.4) | Λ(1.4) |
|---|---|---|---|---|---|
| GR (V7) | 2.0495 | — | 11.686 | — | baseline |
| Emb1 (operator ratio) | 2.0349 | **−0.71%** | 11.637 | **−0.42%** | reduced |
| Emb2 (eff-dim mass) | 2.0378 | **−0.57%** | 11.639 | **−0.40%** | reduced |
| F-multiplier (III V7) | 2.0677 | +0.89% | 11.745 | +0.50% | +5% |

Λ(1.4) sign follows R (Λ ~ k₂R⁵): ΔR<0 ⇒ **Λ reduced** (full value = Paper IV
re-run, item 2). Representative headline shift: **ΔM_max ≈ −0.6%, ΔR(1.4) ≈
−0.4%**; the two embeddings bracket it (−0.57% to −0.71%).

## Why it's trustworthy (not tuned)
- **Gate 1** δ=0 → V7 exactly: M_max 2.0495 (0.00%), R(1.4) 11.686 (0.00%).
- **Gate 3** tolerance: ΔM=5e−9, ΔR=6e−6 km under 10× tighter tol.
- **Two inequivalent embeddings agree in sign** and to ~0.15 pp in magnitude.
  The relativistic-dressing ambiguity does NOT flip the sign.
- **Emb2 r⋆ sensitivity flat**: ΔM_max=−0.57% for r⋆=5/10/20 km — because δ
  is saturated (≈−3A=const) through the NS bulk, so the operator is just
  "(3−3A)-dimensional gravity"; the ln(r/r⋆) term is negligible.
- **Physical mechanism**: operator with δ<0 gives R_op = a_op/a_N > 1
  (stronger field; R_op→3/(3+δ)=1/(1−A) at the center) → more compact star,
  lower M_max. The F-multiplier did the opposite by weakening dP/dr.

## Consequence for the series (per audit)
- **Paper III:** rewrite M_max/R/Λ with the contraction sign. RDT does NOT
  raise M_max (so it does NOT help heavy pulsars — it slightly lowers the
  ceiling). The +0.89% "raise" of V7 is withdrawn; correct sign is ≈ −0.6%.
- **Paper IV:** Λ enhancement (+5%) flips to a **reduction** (~−2%), via the
  same backgrounds; re-run pipeline unchanged (item 2).
- Mechanism now matches I–II everywhere: every paper contracts.

## Open methodological note
Two embeddings agree, so the sign is robust. They differ ~0.15 pp in
magnitude (the genuine relativistic-dressing freedom). For the III rewrite,
quote Emb1 (primary, cleanest) as the value and Emb2 as the robustness
bracket, OR quote the range. Decision deferred to Chris at check-in.

Files: tov_operator.py, run_keystone.py, derivation.md (this folder).
Reproduce: `python3 run_keystone.py` (single-thread; ~30 s).

## Post-review verification (Checks A/B/C) — 2026-06-14
- **Check A (gate independence) PASS.** The modified emb1/emb2 RHS at A=1e-15
  (forced past the A==0 V7 fallback) reproduces the V7 2-ODE solver to ~7e-8 in
  M and ~2e-5 in R (integrator level, from 4-ODE vs 2-ODE adaptive stepping —
  10³–10⁴× below the effect). Linear-scaling: (ΔM/M)/A is constant at
  −0.966→−0.971 across A=1e-5..1e-2 ⇒ ΔM/M ∝ A ⇒ the path → V7 as δ→0. The
  contraction at A=0.0087 is produced by the SAME code path, only A varied.
  (The earlier "ratio 12" at A=1e-9 was integrator noise below the floor, not
  a failure.)  File: gate_independence.py.
- **Check B (embedding non-uniqueness, relabeled).** The Emb1–Emb2 spread
  (~0.15 pp) is NOT robustness/noise — it is the genuine non-uniqueness of the
  relativistic embedding of the I–II operator (the residue of the series'
  covariant-formulation gap, shrunk from a SIGN ambiguity to a sub-0.1%
  MAGNITUDE ambiguity). Quote Emb1 (primary), report Emb2 as the spread.
- **Check C (A95(α) ceiling, Paper II envelope format).** boundary_sweep.py:
  | α | A95 | ΔM_max (Emb1/Emb2) | ΔR(1.4) (Emb1/Emb2) |
  |---|---|---|---|
  | 0.12  | 0.0087 | −0.71% / −0.57% | −0.42% / −0.40% |
  | 0.15  | 0.0083 | −0.69% / −0.55% | −0.41% / −0.38% |
  | 0.334 | 0.0061 | −0.51% / −0.41% | −0.30% / −0.29% |
  **CEILING (max contraction, at α=0.12, Emb1): |ΔM_max| ≤ 0.71%,
  |ΔR(1.4)| ≤ 0.42%.** Λ(1.4) reduced (full value from Paper IV re-run).
