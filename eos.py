"""
NS V7 — SLy4 equation of state (Douchin & Haensel 2001 tables).

Loads the manually entered D&H table (037/eos_sly4_dh2001_v2.dat,
columns: rho [g/cm^3], P [dyne/cm^2], epsilon [erg/cm^3], n_b [1/cm^3])
with log-log interpolation in P.

Two epsilon conventions, selectable:
  "rho_c2": eps = rho * c^2.  In D&H 2001 the tabulated rho IS the total
            mass-energy density (rest mass + internal), so this should be
            the correct convention.
  "table":  the v2 file's 'recalculated' epsilon column — suspected of
            double-counting internal energy (V6's residual 5-6% radius
            error). Kept for the A/B test in validate.py.

Also loads the CompOSE RG(SLY4) reference mass-radius curve (eos.mr,
extracted from 025/eos.zip) as an independent cross-check; note RG(SLY4)
is the Gulminelli-Raduta functional with the SLy4 interaction — close to,
but not identical to, D&H SLy (tolerance ~1-2%).
"""

import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
TABLE = os.path.join(HERE, "data", "eos_sly4_dh2001_v2.dat")
MR_REF = os.path.join(HERE, "data", "eos.mr")

C_CGS = 2.99792458e10


class SLy4:
    def __init__(self, epsilon_mode="rho_c2"):
        d = np.loadtxt(TABLE)
        # ensure strictly increasing in P (required for interpolation)
        order = np.argsort(d[:, 1])
        d = d[order]
        keep = np.concatenate([[True], np.diff(d[:, 1]) > 0])
        d = d[keep]
        self.rho = d[:, 0]
        self.P = d[:, 1]
        self.eps = d[:, 0] * C_CGS**2 if epsilon_mode == "rho_c2" else d[:, 2]
        self.mode = epsilon_mode
        self._lp = np.log(self.P)
        self._lr = np.log(self.rho)
        self._le = np.log(self.eps)

    def rho_of_P(self, P):
        return np.exp(np.interp(np.log(P), self._lp, self._lr))

    def eps_of_P(self, P):
        return np.exp(np.interp(np.log(P), self._lp, self._le))

    def P_of_rho(self, rho):
        return np.exp(np.interp(np.log(rho), self._lr, self._lp))

    @property
    def P_min(self):
        return self.P[0]

    @property
    def P_max(self):
        return self.P[-1]


def compose_reference():
    """CompOSE RG(SLY4) (R [km], M [Msun]) reference curve."""
    d = np.loadtxt(MR_REF)
    return d[:, 0], d[:, 1]


if __name__ == "__main__":
    for mode in ("rho_c2", "table"):
        e = SLy4(mode)
        # sanity: monotone, causal sound speed (finite-difference)
        lp, le = np.log(e.P), np.log(e.eps)
        cs2_over_c2 = np.diff(e.P) / np.diff(e.eps)  # dP/deps, in c^2 units
        print(f"mode={mode:7s}: {len(e.P)} pts, "
              f"rho {e.rho[0]:.2e}..{e.rho[-1]:.2e}, "
              f"P monotone: {np.all(np.diff(e.P) > 0)}, "
              f"eps monotone: {np.all(np.diff(e.eps) > 0)}, "
              f"max cs^2/c^2 = {cs2_over_c2.max():.3f}, "
              f"eps/rho_c2 range: {(e.eps/(e.rho*C_CGS**2)).min():.4f}"
              f"..{(e.eps/(e.rho*C_CGS**2)).max():.4f}")
    R, M = compose_reference()
    i = np.argmin(np.abs(M - 1.4))
    print(f"CompOSE ref: {len(M)} pts, M_max = {M.max():.3f}, "
          f"R(~1.4) = {R[i]:.2f} km")


# CompOSE RG(SLY4) full-range table (includes outer crust below neutron
# drip, which the D&H 77-point table lacks — that gap costs ~0.3-0.4 km
# of radius). Format per CompOSE manual: eos.thermo rows are
# i_T i_nb i_yq Q1..Q7 [+extra], with Q1 = p/nb [MeV], Q7 = e/(nb m_n) - 1;
# eos.nb lists the nb grid [fm^-3].
MEV_FM3_TO_CGS = 1.602176634e33   # erg/cm^3 (or dyne/cm^2)
M_U = 1.66053906660e-24           # g


class CompOSE_SLy4:
    def __init__(self):
        nb_lines = open(os.path.join(HERE, "data", "eos.nb")).read().split()
        nb = np.array([float(x) for x in nb_lines[2:]])
        rows = []
        with open(os.path.join(HERE, "data", "eos.thermo")) as f:
            header = f.readline().split()
            m_n = float(header[0])  # MeV
            for line in f:
                w = line.split()
                if len(w) < 10:
                    continue
                i_nb = int(w[1])
                Q1, Q7 = float(w[3]), float(w[9])
                rows.append((i_nb, Q1, Q7))
        idx0 = min(r[0] for r in rows)
        P, eps, rho = [], [], []
        for i_nb, Q1, Q7 in rows:
            n = nb[i_nb - idx0]
            P.append(Q1 * n * MEV_FM3_TO_CGS)
            eps.append(n * m_n * (1.0 + Q7) * MEV_FM3_TO_CGS)
            rho.append(n * 1.0e39 * M_U)
        order = np.argsort(P)
        P = np.array(P)[order]; eps = np.array(eps)[order]; rho = np.array(rho)[order]
        keep = np.concatenate([[True], np.diff(P) > 0])
        self.P, self.eps, self.rho = P[keep], eps[keep], rho[keep]
        self.mode = "compose"
        self._lp, self._lr, self._le = np.log(self.P), np.log(self.rho), np.log(self.eps)

    rho_of_P = SLy4.rho_of_P
    eps_of_P = SLy4.eps_of_P
    P_of_rho = SLy4.P_of_rho
    P_min = SLy4.P_min
    P_max = SLy4.P_max
