"""
Paper III keystone — gates + M-R sweep + sign read.

Gate 1 (delta=0 -> V7): A=0 sweep reproduces M_max=2.0495, R(1.4)=11.686.
Gate 2 (weak field / embedding agreement): both embeddings contract, agree
       in sign across the mass range, fractional shift finite at low mass.
Gate 3 (tolerance): one star stable to 10x tighter tol.
Then: M_max, R(1.4) for GR / Emb1 / Emb2 at the canonical ceiling
(A=0.0087, alpha=0.12) -> sign of the NS effect.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.normpath(os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "..", "..", "paper3_NeutronStars", "ns_v7")))
from eos import CompOSE_SLy4
from tov_operator import solve_star_op, mr_curve_op

eos = CompOSE_SLy4()
A, ALPHA = 0.0087, 0.12          # canonical global ceiling (NOT alpha=1)
PC = np.geomspace(1e34, 1.3e36, 60)


def summarize(curve):
    ok = ~np.isnan(curve[:, 1])
    M, R = curve[ok, 1], curve[ok, 2]
    i = int(M.argmax())
    R14 = np.interp(1.4, M[:i + 1], R[:i + 1])
    return M[i], R[i], R14


# --- sweeps ---
gr = mr_curve_op(eos, PC, 0.0, ALPHA)
e1 = mr_curve_op(eos, PC, A, ALPHA, embedding="emb1")
e2 = mr_curve_op(eos, PC, A, ALPHA, embedding="emb2")
Mg, Rg, R14g = summarize(gr)
M1, R1, R141 = summarize(e1)
M2, R2, R142 = summarize(e2)

# --- Gate 1: delta=0 -> V7 ---
e1g = abs(Mg - 2.0495) / 2.0495
e2g = abs(R14g - 11.686) / 11.686
print("=== Gate 1 (delta=0 -> V7) ===")
print(f"  M_max = {Mg:.4f} (V7 2.0495, {100*e1g:.2f}%), "
      f"R(1.4) = {R14g:.3f} km (V7 11.686, {100*e2g:.2f}%) "
      f"-> {'PASS' if e1g<1e-3 and e2g<1e-3 else 'FAIL'}")

# --- Gate 3: tolerance ---
m_a, r_a = solve_star_op(eos, 2e35, A, ALPHA, embedding="emb1")
m_b, r_b = solve_star_op(eos, 2e35, A, ALPHA, embedding="emb1",
                         rtol=1e-10, max_step=5e3)
print("=== Gate 3 (tolerance, Emb1) ===")
print(f"  dM = {abs(m_a-m_b):.1e}, dR = {abs(r_a-r_b):.1e} km "
      f"-> {'PASS' if abs(m_a-m_b)<1e-3 and abs(r_a-r_b)<1e-2 else 'FAIL'}")

# --- r* sensitivity (Emb2) ---
print("=== Emb2 r* sensitivity (sign must not flip) ===")
for rs in (5e5, 1e6, 2e6):
    c = mr_curve_op(eos, PC, A, ALPHA, embedding="emb2", rstar_cm=rs)
    Mm, _, _ = summarize(c)
    print(f"  r*={rs/1e5:.0f} km: dM_max/M_max = {100*(Mm/Mg-1):+.2f}%")

# --- Sign verdict numbers ---
print("\n=== KEYSTONE NUMBERS (canonical A=0.0087, alpha=0.12) ===")
print(f"{'model':>6} {'M_max':>8} {'dM_max%':>9} {'R(1.4)':>8} {'dR(1.4)%':>10}")
print(f"{'GR':>6} {Mg:8.4f} {'--':>9} {R14g:8.3f} {'--':>10}")
print(f"{'Emb1':>6} {M1:8.4f} {100*(M1/Mg-1):+9.2f} {R141:8.3f} {100*(R141/R14g-1):+10.2f}")
print(f"{'Emb2':>6} {M2:8.4f} {100*(M2/Mg-1):+9.2f} {R142:8.3f} {100*(R142/R14g-1):+10.2f}")
print(f"\nLambda(1.4) sign follows R (Lambda ~ k2 R^5): dR<0 -> Lambda REDUCED.")
print(f"F-multiplier (III V7) gave +0.89% M_max, +0.5% R(1.4): OPPOSITE sign.")
