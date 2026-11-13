"""Legacy nominal scoring of a complete, provenance-checked solar scan.
Requires --input and a new --figures directory. The helioseismic term is a
chosen deterministic proxy, not an observed Gaussian residual. Source-precision
scoring and full interval extraction remain step 11; no coverage is claimed.
"""

import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))


SIGMA_HELIO = 2.0e-3

# Experimental ratios (Borexino 2018 / B16-GS98, verified 2026-06-10).
# pep is the HZ-conditioned value. Two error models:
#   "meas":   experimental uncertainties only
#   "theory": legacy selective addition of 6% and 12% for Be7 and B8
EXP_R = {"pp": 1.02, "pep": 0.88, "be7": 1.01, "b8": 1.04}
SIG = {
    "meas":   {"pp": 0.10, "pep": 0.15, "be7": 0.03, "b8": 0.08},
    "theory": {"pp": 0.10, "pep": 0.15,
               "be7": np.hypot(0.03, 0.06), "b8": np.hypot(0.08, 0.12)},
}
ERROR_MODEL = "theory"   # used for the headline; both are reported

def main(csv_path, figure_dir):
    from scan import load_complete
    rows = load_complete(csv_path)
    if len({r["A"] for r in rows}) < 2 or len({r["alpha"] for r in rows}) < 2:
        raise ValueError("A map needs at least two points per axis")
    # No figures are touched until the full manifested grid is accepted.
    if os.path.exists(figure_dir):
        raise ValueError("Use a new figure directory; historical output is preserved")
    FIGDIR = figure_dir
    os.makedirs(FIGDIR)
    A = np.array([float(r["A"]) for r in rows])
    al = np.array([float(r["alpha"]) for r in rows])
    Tc_flat = np.array([float(r["Tc_ratio"]) for r in rows])
    dcs_flat = np.array([float(r["max_abs_dcs"]) for r in rows])


    def chi2_nu_with(sig):
        out = np.zeros(len(rows))
        for b in ("pp", "pep", "be7", "b8"):
            Rb = np.array([float(r[f"R_{b}"]) for r in rows])
            out += ((Rb - EXP_R[b]) / sig[b]) ** 2
        return out


    chi2_nu_flat = chi2_nu_with(SIG[ERROR_MODEL])
    chi2_nu_meas_flat = chi2_nu_with(SIG["meas"])

    A_vals = np.unique(A)
    al_vals = np.unique(al)
    idx = {(a, l): i for i, (a, l) in enumerate(zip(A, al))}


    def grid(arr):
        return np.array([[arr[idx[(a, l)]] for l in al_vals] for a in A_vals])


    CHI2_NU = grid(chi2_nu_flat)
    DCS = grid(dcs_flat)
    TC = grid(Tc_flat)
    CHI2_CO = CHI2_NU + (DCS / SIGMA_HELIO) ** 2


    def upper_limit_A(CHI2):
        """95% CL one-sided upper limit on A at each alpha (profile in A),
        Delta-chi^2 = 2.71 above the profile minimum. NaN = beyond grid."""
        out = []
        for j in range(len(al_vals)):
            prof = CHI2[:, j]
            jmin = int(np.argmin(prof))
            target = prof[jmin] + 2.71
            a95 = np.nan
            for i in range(jmin, len(A_vals) - 1):
                if prof[i] <= target <= prof[i + 1]:
                    f = (target - prof[i]) / (prof[i + 1] - prof[i])
                    a95 = A_vals[i] + f * (A_vals[i + 1] - A_vals[i])
                    break
            out.append(a95)
        return np.array(out)


    def report(name, CHI2):
        ib = np.unravel_index(np.argmin(CHI2), CHI2.shape)
        c_min = CHI2[ib]
        c_null = CHI2[A_vals == 0.0, :].min()
        A95 = upper_limit_A(CHI2)
        print(f"\n--- {name} ---")
        print(f"chi2(null) = {c_null:.3f}; chi2(best) = {c_min:.3f} at "
              f"A = {A_vals[ib[0]]:.4f}, alpha = {al_vals[ib[1]]:.4f}; "
              f"Delta chi2 = {c_null - c_min:.3f} (nominal score; no calibrated significance)")
        print("A_95(alpha):")
        for k in range(0, len(al_vals), max(1, len(al_vals)//8)):
            lim = f"A < {A95[k]:.4f}" if np.isfinite(A95[k]) else "> grid (0.04)"
            print(f"  alpha = {al_vals[k]:.4f}:  {lim}")
        print(f"Sampled nominal envelope: A < {np.nanmax(A95):.4f} at all alpha in grid"
              if np.all(np.isfinite(A95)) else
              "Nominal envelope: grid-limited at small alpha")
        # Papers 2-3 fiducial exclusion
        i_f = int(np.argmin(np.abs(A_vals - 0.0334)))
        for l_t in (0.15, 0.30):
            j = int(np.argmin(np.abs(al_vals - l_t)))
            d = CHI2[i_f, j] - c_min
            print(f"Fiducial A={A_vals[i_f]:.4f}, alpha={al_vals[j]:.4f}: "
                  f"Delta chi2 = {d:.1f} (score difference only)")
        return A95, ib


    CHI2_NU_MEAS = grid(chi2_nu_meas_flat)
    CHI2_CO_MEAS = CHI2_NU_MEAS + (DCS / SIGMA_HELIO) ** 2

    print("=== Paper-numbers block ===")
    print(f"(helioseismic proxy tolerance sigma_h = {SIGMA_HELIO:.1e}; "
          f"headline error model = {ERROR_MODEL})")
    A95_nu, ib_nu = report("Neutrinos only", CHI2_NU)
    A95_co, ib_co = report("Combined (nu + helio proxy)", CHI2_CO)
    A95_co_m, _ = report("Combined, measurement-only errors (comparison)",
                         CHI2_CO_MEAS)

    print(f"\nBest fit (combined): Tc_ratio = {TC[ib_co]:.5f}, "
          f"max|dcs/cs| = {DCS[ib_co]:.2e}")

    # --- Figures ---------------------------------------------------------------
    if len(al_vals) < 2 or len(A_vals) < 2:
        raise ValueError("A map needs at least two points per axis")
    AL, AA = np.meshgrid(al_vals, A_vals)

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), sharey=True)
    for ax, CHI2, A95, ib, title in [
            (axes[0], CHI2_NU, A95_nu, ib_nu, "Neutrino fluxes only"),
            (axes[1], CHI2_CO, A95_co, ib_co, "Combined: neutrinos + helioseismic proxy")]:
        D = CHI2 - CHI2.min()
        pc = ax.pcolormesh(AL, AA, np.clip(D, 0, 30), shading="auto",
                           cmap="viridis_r")
        cs = ax.contour(AL, AA, D, levels=[2.30, 6.18], colors=["w", "r"],
                        linewidths=1.5)
        ax.clabel(cs, fmt={2.30: r"$\Delta Q=2.30$", 6.18: r"$\Delta Q=6.18$"})
        ax.plot(al_vals[ib[1]], A_vals[ib[0]], "w*", ms=13, mec="k",
                label="best fit")
        ok = np.isfinite(A95)
        ax.plot(al_vals[ok], A95[ok], "k--", lw=1.5, label=r"$A_{95}(\alpha)$")
        ax.plot([0.15], [0.0334], "rx", ms=10, mew=2.5,
                label="Papers 2-3 fiducial")
        ax.set_xscale("log")
        ax.set_xlabel(r"$\alpha$")
        ax.set_title(title, fontsize=11)
    axes[0].set_ylabel(r"$A$")
    axes[0].legend(loc="upper left", fontsize=8)
    fig.colorbar(pc, ax=axes, label=r"$\Delta\chi^2$ (clipped at 30)")
    plt.savefig(os.path.join(FIGDIR, "fig_chi2_map.pdf"), bbox_inches="tight")
    plt.close()

    plt.figure(figsize=(7, 4.5))
    ok = np.isfinite(A95_nu)
    plt.plot(al_vals[ok], A95_nu[ok], "C0-", lw=2, label=r"neutrinos only")
    plt.plot(al_vals, A95_co, "k-", lw=2, label=r"combined ($+$ helio proxy)")
    plt.fill_between(al_vals, np.where(np.isfinite(A95_co), A95_co, 0.04), 0.04,
                     alpha=0.15, color="r")
    plt.plot([0.15], [0.0334], "rx", ms=10, mew=2.5, label="Papers 2-3 fiducial")
    plt.xscale("log")
    plt.xlabel(r"$\alpha$")
    plt.ylabel(r"$A_{95}$")
    plt.ylim(0, 0.04)
    plt.title(r"Legacy nominal score thresholds (coverage not calibrated)")
    plt.legend(fontsize=9)
    plt.tight_layout()
    plt.savefig(os.path.join(FIGDIR, "fig_A95.pdf"))
    plt.close()

    print(f"\nFigures written to {FIGDIR}/ (fig_chi2_map.pdf, fig_A95.pdf)")


if __name__ == "__main__":
    import argparse
    parser=argparse.ArgumentParser(description="Analyze a complete manifested scan with legacy nominal scoring")
    parser.add_argument("--input",required=True)
    parser.add_argument("--figures",required=True)
    args=parser.parse_args()
    main(args.input,args.figures)
