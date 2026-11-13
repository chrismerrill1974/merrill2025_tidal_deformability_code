"""
g_regression_v2.py  --  VERIFIED-dataset regression of RDT's terrestrial
prediction against the modern precision-G determinations.

Supersedes g_regression.py (which used best-recall G values). Every G value,
uncertainty, and source-mass material here is pinned to a primary or
authoritative source:

  * Values/uncertainties: Schlamminger, Gundlach & Newman, "Recent measurements
    of the gravitational constant as a function of time" (2015), Table II
    [written by three of the actual G experimentalists; matches the CODATA
    input set]. NIST/BIPM-2026 from Schlamminger et al., Metrologia 63 (2026)
    025012 (open access), Eq. (64) / Table 16.
  * Source-mass materials (the RDT variable): the experiment papers and the
    Rothleitner-Schlamminger 2017 review --
      NIST-82  : two 10.5 kg tungsten balls          (Luther & Towler 1982)
      LANL-97  : tungsten                             (Bagley & Luther 1997)
      UWash-00 : eight 316 stainless-steel spheres    (Gundlach & Merkowitz 2000)
      BIPM-01  : copper source masses                 (Quinn et al. 2001)
      UWup-02  : two ~576 kg brass masses             (Kleinevoss 2002)
      UZur-06  : mercury source on copper test masses (Schlamminger et al. 2006)
      JILA-10  : four 120 kg tungsten cylinders       (Parks & Faller 2010)
      BIPM-14  : copper source masses                 (Quinn et al. 2013/14)
      LENS-14  : tungsten alloy INERMET IT180, ~18.0  (Rosi et al. 2014)
      HUST-18  : stainless-steel source spheres       (Li et al. 2018)
      NIST-26  : copper + sapphire SM sets, same Mark II BIPM balance
                 (Schlamminger et al. 2026); consensus 6.67387 taken at copper
                 density (copper-free dominates the consensus). The internal
                 copper-vs-sapphire swap is analysed SEPARATELY below.

CHANGES vs g_regression.py:
  - UWash-00 corrected 6.674230 -> 6.674255 (post-correction final value)
  - NIST/BIPM-26 corrected 6.67366(20) -> 6.67387(38) (Metrologia 2026)
  - LENS-14 source density 19.3 -> 18.0 (tungsten heavy alloy, not pure W)
  - HUST-18 ADDED as two points (TOS 6.674184, AAF 6.674484; steel, rho 8.0)
  - The "controlled swap" is now the NIST-2026 internal copper-vs-sapphire 2x2
    (same apparatus, same time) -- NOT the confounded BIPM-14(France,Cu) ->
    NIST(sapphire) comparison used in v1, which mixed a ~250 ppm apparatus/
    analysis revision with the density change.

CAVEATS (load-bearing, unchanged): the G scatter is ~5x the quoted
uncertainties (Birge ratio ~5-6), widely attributed to per-experiment
SYSTEMATICS. We inflate the fit uncertainty by the Birge ratio. Small N,
apparatus-correlated points (the two BIPM, the two HUST, the NIST 2x2 share
groups/hardware) -> the cross-experiment A-bound is INDICATIVE, not rigorous.
This is a screen + a proposal, not a measurement of A.
"""
import numpy as np

RHO0 = 150.0
A_SOLAR = 0.0087

# name, year, G[1e-11], u[1e-11], source-mass material, rho_src[g/cc]
DATA = [
    ("NIST-82",   1982, 6.67260,  0.00050,  "tungsten",         19.30),
    ("LANL-97",   1997, 6.67400,  0.00070,  "tungsten",         19.30),
    ("UWash-00",  2000, 6.674255, 0.000092, "stainless steel",   8.00),
    ("BIPM-01",   2001, 6.67559,  0.00027,  "copper",            8.96),
    ("UWup-02",   2002, 6.67422,  0.00098,  "brass",             8.50),
    ("UZur-06",   2006, 6.67425,  0.00012,  "mercury",          13.53),
    ("JILA-10",   2010, 6.67234,  0.00014,  "tungsten",         19.30),
    ("BIPM-14",   2014, 6.67554,  0.00016,  "copper",            8.96),
    ("LENS-14",   2014, 6.67191,  0.00099,  "tungsten alloy",   18.00),
    ("HUST-18TOS",2018, 6.674184, 0.000078, "stainless steel",   8.00),
    ("HUST-18AAF",2018, 6.674484, 0.000078, "stainless steel",   8.00),
    ("NIST-26",   2026, 6.67387,  0.00038,  "copper(cons.)",     8.96),
]

# NIST-2026 internal 2x2 (Schlamminger et al. 2026, Table 16): the one clean,
# same-apparatus, same-time source-density swap. G[1e-11], u_r[ppm], SM material.
NIST26_2x2 = {
    ("copper",   "servo"): (6.673642, 23, 8.96),
    ("copper",   "free"):  (6.674021, 30, 8.96),
    ("sapphire", "servo"): (6.672637, 38, 3.98),
    ("sapphire", "free"):  (6.673636, 94, 3.98),
}


def s_shape(rho, alpha):
    x = (rho / RHO0) ** alpha
    return x / (1.0 + x)          # sigma = A * s_shape


def wls(x, y, u):
    """Weighted least squares y = a + b x; returns a, b, sigma_b, chi2/dof."""
    w = 1.0 / u**2
    Sw, Sx = w.sum(), (w*x).sum()
    Sy, Sxx = (w*y).sum(), (w*x*x).sum()
    Sxy = (w*x*y).sum()
    D = Sw*Sxx - Sx*Sx
    b = (Sw*Sxy - Sx*Sy) / D
    a = (Sy - b*Sx) / Sw
    var_b = Sw / D
    resid = y - (a + b*x)
    chi2 = (w*resid**2).sum()
    dof = len(x) - 2
    return a, b, np.sqrt(var_b), chi2/dof


if __name__ == "__main__":
    G = np.array([d[2] for d in DATA])
    u = np.array([d[3] for d in DATA])
    rho = np.array([d[5] for d in DATA])
    Gref = G.mean()

    print("VERIFIED modern-G dataset (source-mass density is the RDT variable):")
    print(f"  {'exp':>12} {'year':>5} {'G[1e-11]':>10} {'u[ppm]':>7} "
          f"{'source':>16} {'rho':>6}")
    for nm, yr, g, uu, mat, r in DATA:
        print(f"  {nm:>12} {yr:>5} {g:>10.5f} {1e6*uu/g:>7.0f} {mat:>16} {r:>6.2f}")

    # over-scatter (Birge ratio under a constant-G model)
    w = 1/u**2
    Gw = (w*G).sum()/w.sum()
    birge = np.sqrt(((w*(G-Gw)**2).sum())/(len(G)-1))
    print(f"\nBirge ratio (scatter/quoted unc, constant-G model) = {birge:.1f}  "
          f"[N={len(G)}; scatter is systematics-dominated]")

    # ---- model-free: fixed-density spans (scatter is NOT density-driven) ----
    def span(mask, label):
        gg = G[mask]
        print(f"  {label:>34}: {gg.min():.5f}..{gg.max():.5f}"
              f"  = {1e6*(gg.max()-gg.min())/Gref:>4.0f} ppm  (N={mask.sum()})")
    print("\nFixed-density spans (same source density, different experiments):")
    span(np.isclose(rho, 19.30), "pure tungsten (rho=19.3)")
    span(np.isin(rho, [19.30, 18.00]), "tungsten+alloy (rho=18-19.3)")
    span(np.isclose(rho, 8.96), "copper (rho=8.96)")
    span(np.isclose(rho, 8.00), "stainless steel (rho=8.0)")
    print("  RDT predicts <~10 ppm across each of these narrow density ranges,")
    print("  so the fixed-density spans (45-313 ppm) are NOT density-driven.")

    # ---- RDT fit  G = G_N(1 + A*s(rho;alpha)) ----
    print("\nRDT fit  G = G_N(1 + A*s(rho;alpha))  [A>0 = RDT; A<0 unphysical]:")
    print(f"  errors inflated by Birge ratio ({birge:.1f}x) for honest bounds\n")
    print(f"  {'alpha':>6} {'A_fit':>10} {'+/-(infl)':>11} {'95% UL on A':>12} "
          f"{'vs solar 0.0087':>16}")
    for alpha in (0.02, 0.12, 0.334, 0.50):
        s = s_shape(rho, alpha)
        a, b, sb, chi2dof = wls(s, G, u)
        A_fit = b / a
        sA = (sb / a) * birge
        A_ul = max(A_fit, 0.0) + 1.645*sA
        tag = "TIGHTER" if A_ul < A_SOLAR else "looser"
        print(f"  {alpha:>6} {A_fit:>+10.4f} {sA:>11.4f} {A_ul:>12.4f} "
              f"{tag:>16}")

    # ---- the clean controlled swap: NIST-2026 internal Cu vs sapphire 2x2 ----
    print("\n" + "="*70)
    print("CONTROLLED SWAP -- NIST 2026, SAME Mark II BIPM balance, same run")
    print("(Schlamminger et al., Metrologia 63 (2026) 025012, Table 16):")
    print("="*70)
    (cs, _, _) = NIST26_2x2[("copper", "servo")]
    (cf, _, _) = NIST26_2x2[("copper", "free")]
    (ss, _, _) = NIST26_2x2[("sapphire", "servo")]
    (sf, _, _) = NIST26_2x2[("sapphire", "free")]
    g0 = 6.6737
    print(f"                       servo        free-deflection")
    print(f"   copper  (8.96)   {cs:.6f}      {cf:.6f}")
    print(f"   sapphire(3.98)   {ss:.6f}      {sf:.6f}")
    print(f"\n  Same-method copper-minus-sapphire (the RDT density observable):")
    print(f"    servo: {1e6*(cs-ss)/g0:+.0f} ppm     free: {1e6*(cf-sf)/g0:+.0f} ppm"
          f"   (copper higher = RDT direction, both methods)")
    print(f"  Fixed-density method effect (free-minus-servo, NOT density):")
    print(f"    copper: {1e6*(cf-cs)/g0:+.0f} ppm    sapphire: {1e6*(sf-ss)/g0:+.0f} ppm"
          f"   (comparable to the density effect -> entangled)")
    print(f"\n  A implied by the same-method copper-sapphire difference:")
    for alpha in (0.12, 0.334, 0.50):
        ds = s_shape(8.96, alpha) - s_shape(3.98, alpha)
        print(f"    alpha={alpha:<5}  from servo: A={ (cs-ss)/g0/ds:.4f}"
              f"   from free: A={ (cf-sf)/g0/ds:.4f}"
              f"   (RDT@ceiling predicts {1e6*A_SOLAR*ds:+.0f} ppm)")
    print("  -> copper (denser) is higher than sapphire in BOTH methods (RDT")
    print("     direction), implied A ~0.001-0.006 (at/below the solar ceiling),")
    print("     BUT the fixed-density servo-vs-free systematic is just as large,")
    print("     so density cannot be cleanly separated from method. The authors")
    print("     model it as correlated systematics + 'dark uncertainty' and")
    print("     attribute the shift to instrument, not composition.")

    # ---- second within-group, fixed-density disagreement: HUST-18 ----
    tos = [d[2] for d in DATA if d[0] == "HUST-18TOS"][0]
    aaf = [d[2] for d in DATA if d[0] == "HUST-18AAF"][0]
    print(f"\nHUST-18 (same group, same steel SM rho=8.0, two methods):")
    print(f"  TOS {tos:.6f}  vs  AAF {aaf:.6f} : {1e6*(aaf-tos)/g0:+.0f} ppm at "
          f"FIXED density -> a second within-group method systematic.")
