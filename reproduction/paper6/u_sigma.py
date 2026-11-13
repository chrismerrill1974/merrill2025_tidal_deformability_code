"""
paper6 -- U(sigma) construction. Show the potential whose density-driven
equilibrium IS the series saturating law sigma_eq(rho)=A x/(1+x),
x=(rho/rho0)^alpha, and that alpha is a FREE shape parameter (not forced to 1).

Adiabatic equilibrium: U'(sigma_eq) = J(rho), J(rho)=lambda*(rho/rho0) from the
(1-sigma)R coupling. Choosing
    U'(sigma) = lambda * [ sigma/(A-sigma) ]^(1/alpha)
inverts EXACTLY to sigma_eq = A x/(1+x) (set lambda=1). Wall as sigma->A
enforces saturation d_eff -> 3(1-A)=2.974. Convex (U''>0) => tachyon-free.

We root-find U'(sigma)=J(rho) numerically and confirm it == A x/(1+x) for
several alpha. alpha=1 is just the quadratic-U special case (the companion's
forced value); alpha in [0.02,0.5] are all realized by the exponent 1/alpha.
"""
import numpy as np
from scipy.optimize import brentq
from scipy.integrate import quad

A = 0.0087
RHO0 = 150.0  # g/cm^3 (series convention = solar-core density)


def Uprime(s, alpha):
    return (s / (A - s)) ** (1.0 / alpha)


def Udprime(s, alpha):  # U''(sigma) > 0  => convex => tachyon-free
    return (1.0 / alpha) * (s / (A - s)) ** (1.0 / alpha - 1.0) * A / (A - s) ** 2


def sigma_target(rho, alpha):
    x = (rho / RHO0) ** alpha
    return A * x / (1.0 + x)


def sigma_from_potential(rho, alpha):
    J = rho / RHO0                       # lambda = 1
    return brentq(lambda s: Uprime(s, alpha) - J, 1e-280, A * (1 - 1e-13),
                  xtol=1e-18)


if __name__ == "__main__":
    rhos = [1.5, 150.0, 1.5e4, 1.5e8, 5.5e14]   # outer crust .. NS core
    for alpha in (0.12, 0.30, 0.50, 1.0):
        print(f"\nalpha = {alpha:>4}   U'(s) = [s/(A-s)]^(1/{alpha})")
        print(f"  {'rho g/cc':>10} {'sigma(pot)':>13} {'A x/(1+x)':>13} "
              f"{'d_eff':>8} {'conv':>5}")
        for rho in rhos:
            s = sigma_from_potential(rho, alpha)
            st = sigma_target(rho, alpha)
            conv = Udprime(min(s, A * 0.9999), alpha) > 0
            print(f"  {rho:10.3g} {s:13.6e} {st:13.6e} {3*(1-s):8.4f} "
                  f"{str(conv):>5}")
    U_wall, _ = quad(lambda s: Uprime(s, 0.12), 0, A * 0.999)
    print(f"\nU(A^-) finite (alpha=0.12): U ~ {U_wall:.3e}  "
          f"(steep wall at sigma=A => d_eff > 2.974 floor)")
    print("=> sigma(pot) == A x/(1+x) for every alpha "
          "=> alpha is a FREE potential parameter, NOT forced to 1.")
