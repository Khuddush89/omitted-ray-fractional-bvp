"""Check of Theorem 2 against exact solutions (Appendix A).

For h(s) = s^lambda the problem D^alpha u + h = 0, u(0) = 0, D^nu u(1) = int u dA
has the exact solution  u(t) = -I^alpha h(t) + c1 t^nu  with
    I^alpha h(t) = Gamma(lambda+1)/Gamma(lambda+alpha+1) t^(lambda+alpha),
    c1 = (1/Delta) [ 1/(lambda+1) - Gamma(lambda+1)/Gamma(lambda+alpha+1) int tau^(lambda+alpha) dA ].
We compare it with int_0^1 calG(t,s) s^lambda ds computed by adaptive quadrature
(scipy.integrate.quad, break points at s = t and at the atoms of dA).

Output: results/kernel_check.json
"""
import warnings
from math import gamma

import numpy as np
from scipy.integrate import quad, IntegrationWarning

from _setup import save_json
from fracbvp import Problem
from compute_constants import EX1

ALPHAS = (1.08, 1.2, 1.5, 1.75, 2.0)
LAMBDAS = (0.0, 0.5, 1.0, 2.0)
TS = (0.05, 0.2, 0.5, 0.75, 0.9, 1.0)


def exact(p, lam, t):
    a = p.alpha
    c = gamma(lam + 1) / gamma(lam + a + 1)
    c1 = (1.0 / (lam + 1) - c * np.sum(p.w_a * p.xi_a ** (lam + a))) / p.Delta
    return -c * t ** (lam + a) + c1 * t ** p.nu, c1


def main():
    warnings.simplefilter("ignore", IntegrationWarning)
    worst, worst_bc, rows = 0.0, 0.0, []
    for al in ALPHAS:
        p = Problem(al, EX1["xi"], EX1["w"])
        for lam in LAMBDAS:
            for t in TS:
                pts = sorted({min(max(t, 1e-12), 1 - 1e-12), *p.xi_a})
                val, err = quad(lambda s: float(p.calG(t, s)) * s ** lam, 0.0, 1.0,
                                points=pts, epsabs=1e-14, epsrel=1e-12, limit=400)
                ex, c1 = exact(p, lam, t)
                dev = abs(val - ex)
                worst = max(worst, dev)
                rows.append(dict(alpha=al, lam=lam, t=t, quad=val, exact=float(ex), deviation=dev))
            # boundary conditions of the exact solution
            ex_xi = np.array([exact(p, lam, x)[0] for x in p.xi_a])
            _, c1 = exact(p, lam, 1.0)
            lhs = -1.0 / (lam + 1) + c1 * p.Ga             # D^nu u(1) by (17)
            rhs = float(np.sum(p.w_a * ex_xi))             # int u dA
            worst_bc = max(worst_bc, abs(lhs - rhs))
    out = dict(alphas=ALPHAS, lambdas=LAMBDAS, ts=TS, max_deviation=worst,
               max_bc_residual=worst_bc, rows=rows)
    save_json("kernel_check.json", out)
    print("largest deviation of int calG(t,s) s^lambda ds from the exact solution:", worst)
    print("largest boundary-condition residual of the exact solutions:", worst_bc)


if __name__ == "__main__":
    main()
