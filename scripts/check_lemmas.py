"""Numerical checks of Lemma 5 (G1)-(G3) (Table 5), Lemma 8 (sup R, Table 4)
and Proposition 1 (no uniform Harnack constant, Table 6).

Output: results/lemmas.json, results/table_lemmaG.tex,
        results/table_theta.tex, results/table_harnack.tex
"""
import numpy as np

from _setup import save_json, save_text
from fracbvp import Problem, theta_star


def lemma_G(alphas=(1.05, 1.2, 1.5, 1.8, 1.95, 2.0), n=700):
    rows, tex = [], []
    t = np.linspace(0.0, 1.0, n)
    T, S = np.meshgrid(t, t, indexing="ij")
    for al in alphas:
        p = Problem(al, (0.5,), (0.0,))
        nu, Ga = p.nu, p.Ga
        G = p.G(T, S)
        c1 = float(np.max(G - np.minimum(T, S) ** nu / Ga))           # (G1): <= 0
        c2 = float(np.min(G - T ** nu * p.G(1.0, S)))                 # (G3): >= 0
        Q = G[1:, :] / T[1:, :] ** nu
        c3 = float(np.max(np.diff(Q, axis=0)))                        # (G2): <= 0
        c0 = float(np.min(G))                                         # (G1): >= 0
        rows.append(dict(alpha=al, max_G1=c1, min_G3=c2, max_diff_G2=c3, min_G=c0))
        fmt = lambda x: "0" if x == 0 else ("%.2e" % x).replace("e-", "\\times10^{-").replace("e+", "\\times10^{") + "}"
        tex.append(f"${al:.2f}$&${fmt(c1)}$&${fmt(c2)}$&${fmt(c3)}$\\\\")
    return rows, "\n".join(tex)


def sup_R(cases=((1.05, .75, .25), (1.05, .9, .1), (1.2, .75, .25), (1.5, .75, .25),
                 (1.5, 1.0, .5), (1.7, .75, .25), (1.9, .9, .1), (2.0, .75, .25)),
          N=2_000_000):
    """sup of R(s) = G(xi2,s)/G(xi1,s) over the grid s_k = k/N, 0 < k < N (it contains xi2)."""
    s = np.arange(1, N) / N
    rows, tex = [], []
    for al, x1, x2 in cases:
        p = Problem(al, (x1, x2), (1.0, -1.0))
        R = p.G(x2, s) / p.G(x1, s)
        k = int(np.argmax(R))
        th = theta_star(x1, x2, p.nu)
        # largest value away from xi2 (shows that the maximiser is isolated when nu < 1)
        away = np.abs(s - x2) >= 1e-3
        rows.append(dict(alpha=al, xi1=x1, xi2=x2, grid_max=float(R.max()), argmax=float(s[k]),
                         closed_form=th, difference=float(th - R.max()),
                         max_away_from_xi2=float(R[away].max())))
        arg = "%.2f" % s[k] if al < 2 else "---"
        tex.append(f"${al:.2f}$&${x1:.2f}$&${x2:.2f}$&${th:.9f}$&${R.max():.9f}$&${arg}$&"
                   f"${R[away].max():.6f}$\\\\")
    return rows, "\n".join(tex)


def harnack(alphas=(1.08, 1.2, 1.5, 1.8, 2.0), svals=(1e-1, 1e-2, 1e-3, 1e-4), t0=0.75):
    rows, tex = [], []
    for al in alphas:
        p = Problem(al, (0.5,), (0.0,))
        vals = []
        for s in svals:
            r = (s ** p.nu / p.Ga) / float(p.G(t0, s))
            vals.append(r)
        # the asymptotic rate t0^{1-nu} nu^{-1} s^{nu-1}
        rate = [t0 ** (1 - p.nu) / p.nu * s ** (p.nu - 1) for s in svals]
        rows.append(dict(alpha=al, s=list(svals), ratio=vals, asymptotic=rate))
        tex.append(f"${al:g}$&" + "&".join(f"${v:.2f}$" for v in vals) + "\\\\")
    return rows, "\n".join(tex)


def main():
    r1, t1 = lemma_G()
    r2, t2 = sup_R()
    r3, t3 = harnack()
    save_json("lemmas.json", dict(lemma_G=r1, sup_R=r2, harnack=r3))
    save_text("table_lemmaG.tex", t1 + "\n")
    save_text("table_theta.tex", t2 + "\n")
    save_text("table_harnack.tex", t3 + "\n")
    print(t1, t2, t3, sep="\n\n")
    for r in r2:
        print(r)


if __name__ == "__main__":
    main()
