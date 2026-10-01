"""Closed-form constants of the two examples (Section 7), the check of (E1)-(E4),
and Table 3 of the paper (constants as functions of alpha).

Output: results/constants.json and results/table_constants.tex
"""
import numpy as np

from _setup import save_json, save_text
from fracbvp import Problem, theta_star, F_MINUS, F_PLUS, F_LIP

XI = (0.75, 0.25)
EX1 = dict(alpha=1.5, xi=XI, w=(3 / 5, -9 / 50),
           t1=0.25, t2=0.75, b=1.2, c=0.8, A0=0.9, A1=1.0, a=0.6, d=0.7)
EX2 = dict(alpha=1.08, xi=XI, w=(3 / 5, -1 / 50),
           t1=0.25, t2=0.75, b=1.1, c=0.7, A0=0.9, A1=1.0, a=0.9, d=1.2)


def example_constants(ex):
    p = Problem(ex["alpha"], ex["xi"], ex["w"])
    fm, fp = F_MINUS, F_PLUS
    psi1, psi2 = float(p.Psi(ex["t1"])), float(p.Psi(ex["t2"]))
    b, c, A0, A1, a, d = (ex[k] for k in ("b", "c", "A0", "A1", "a", "d"))
    s = np.linspace(1e-9, 1 - 1e-9, 2_000_001)
    KD = p.KA(s) / p.Delta
    out = dict(
        alpha=p.alpha, nu=p.nu, xi=list(p.xi_a), w=list(p.w_a),
        Gamma_alpha=p.Ga, mu=p.mu, Delta=p.Delta,
        theta_star=p.theta_star(), ratio=float(p.w_a[0] / abs(p.w_a[1])),
        Kbar=p.Kbar, m_Psi=p.m_Psi, M_Psi=p.M_Psi,
        Psi_quarter=psi1, Psi_three_quarters=psi2,
        Psi_070=float(p.Psi(0.70)), Psi_1=float(p.Psi(1.0)),
        t_star=p.t_star, max_Psi=float(p.Psi(p.t_star)),
        KA1_over_Delta=float(p.KA(np.array([1.0]))[0] / p.Delta),
        mu_over_Delta_Gamma=p.mu / (p.Delta * p.Ga),
        min_KA_over_Delta_grid=float(KD.min()), argmin_grid=float(s[KD.argmin()]),
        grid="2000001 equally spaced points of [1e-9, 1-1e-9]",
        lipschitz_times_maxPsi=F_LIP * float(p.Psi(p.t_star)),
        q_070_quarter=float(p.Psi(0.70) / p.Psi(0.25)),
        q_quarter_three_quarters=psi1 / psi2,
    )
    E = dict(
        E1_A0psi2=A0 * psi2, E1_b=b, E1_A0psi1=A0 * psi1, E1_a=a,
        E2_maxdev=max(abs(fp - A0), abs(fm - A0)), E2_maxpsi=max(psi1, psi2),
        E2_left=max(abs(fp - A0), abs(fm - A0)) * max(psi1, psi2),
        E2_right=b - A0 * psi2,
        E3_A1psi2=A1 * psi2, E3_c=c, E3_A1psi1=A1 * psi1, E3_d=d,
        E4_fmpsi2=fm * psi2, E4_margin=fm * psi2 - c,
        forced_inward_fppsi2=fp * psi2,
    )
    E["E2_margin"] = E["E2_right"] - E["E2_left"]
    E["all_hold"] = bool(E["E1_A0psi2"] < b and E["E1_A0psi1"] >= a and E["E2_margin"] > 0
                         and E["E3_A1psi2"] > c and E["E3_A1psi1"] <= d and E["E4_margin"] > 0
                         and 0 < c < b)
    out["E"] = E
    out["parameters"] = {k: ex[k] for k in ("t1", "t2", "b", "c", "A0", "A1", "a", "d")}
    return out


def table_constants(alphas=(1.08, 1.2, 1.4, 1.5, 1.6, 1.8, 2.0)):
    rows, tex = [], []
    for al in alphas:
        p = Problem(al, XI, EX1["w"])
        th = p.theta_star()
        ok = (p.w_a[0] / abs(p.w_a[1])) >= th
        r = dict(alpha=al, Gamma=p.Ga, mu=p.mu, Delta=p.Delta, Kbar=p.Kbar,
                 theta_star=th, criterion=bool(ok), max_Psi=float(p.Psi(p.t_star)),
                 t_star=p.t_star)
        rows.append(r)
        tex.append(f"${al:g}$&${p.Ga:.6f}$&${p.mu:.6f}$&${p.Delta:.6f}$&${p.Kbar:.6f}$&"
                   f"${th:.6f}$&{'yes' if ok else 'no'}&${r['max_Psi']:.6f}$&${p.t_star:.6f}$\\\\")
    return rows, "\n".join(tex)


def ex1_weights_at_ex2_order():
    """The weights of Example 1 at alpha = 1.08: K_A and calG change sign."""
    p = Problem(EX2["alpha"], XI, EX1["w"])
    s = np.linspace(0, 1, 4001)
    KA = p.KA(s)
    T, S = np.meshgrid(np.linspace(0, 1, 801), np.linspace(0, 1, 801), indexing="ij")
    cg = p.calG(T, S)
    k = np.unravel_index(np.argmin(cg), cg.shape)
    return dict(alpha=p.alpha, theta_star=p.theta_star(), ratio=10 / 3, Delta=p.Delta,
                KA_at_xi2=float(p.KA(np.array([0.25]))[0]),
                min_KA=float(KA.min()), argmin_KA=float(s[KA.argmin()]),
                calG_1_quarter=float(p.calG(1.0, 0.25)),
                min_calG_grid=float(cg.min()), argmin_calG=[float(T[k]), float(S[k])])


def main():
    res = dict(example1=example_constants(EX1), example2=example_constants(EX2))
    rows, tex = table_constants()
    res["table_constants"] = rows
    res["example1_weights_at_alpha_1_08"] = ex1_weights_at_ex2_order()
    save_json("constants.json", res)
    save_text("table_constants.tex", tex + "\n")
    for key in ("example1", "example2"):
        e = res[key]
        print(f"--- {key}: alpha={e['alpha']}")
        for k in ("Gamma_alpha", "mu", "Delta", "theta_star", "ratio", "Kbar", "Psi_quarter",
                  "Psi_three_quarters", "Psi_070", "Psi_1", "t_star", "max_Psi",
                  "KA1_over_Delta", "min_KA_over_Delta_grid", "lipschitz_times_maxPsi",
                  "q_070_quarter", "q_quarter_three_quarters"):
            print(f"  {k:28s} {e[k]:.10f}")
        for k, v in e["E"].items():
            print(f"  E.{k:26s} {v}")
    print(tex)
    print(res["example1_weights_at_alpha_1_08"])


if __name__ == "__main__":
    main()
