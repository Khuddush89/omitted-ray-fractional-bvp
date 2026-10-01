"""Numerical solutions of the two examples (Section 7) and of Table 7 (several alpha).

Method: Picard iteration u_{k+1} = T_h u_k, u_0 = 0, where T_h is the product-
integration discretisation of T (src/fracbvp.py): f(s,u(s)) is replaced by its
piecewise linear interpolant on the graded mesh s_j = (j/n)^2 (plus the points
1/4, 1/2, 3/4 and the atoms of dA), and the integrals of the kernel against that
interpolant are computed exactly.

Output: results/solutions.json, results/table_solutions.tex,
        results/table_refinement.tex, results/table_example2.tex,
        results/solution_*.npz (profiles used by the figure scripts)
"""
import os

import numpy as np

from _setup import RESULTS, save_json, save_text
from fracbvp import Problem, f_example, F_MINUS, F_PLUS, solve_problem, locate_max
from compute_constants import EX1, EX2

N_STD, Q = 2400, 2.0


def node_value(D, u, x):
    k = int(np.argmin(np.abs(D.s - x)))
    assert abs(D.s[k] - x) < 1e-14, (x, D.s[k])
    return float(u[k])


def summarize(prob, D, u, hist):
    phi = f_example(D.s, u)
    tmax, umax = locate_max(D, u)
    lhs, rhs = D.bc_check(phi)
    ratio = u[1:] / prob.Psi(D.s[1:])
    diffs = np.array(hist)
    rates = diffs[1:] / diffs[:-1]
    return dict(
        nodes=int(D.s.size), iterations=len(hist), last_increment=float(hist[-1]),
        u_quarter=node_value(D, u, 0.25), u_half=node_value(D, u, 0.5),
        u_three_quarters=node_value(D, u, 0.75), u_1=node_value(D, u, 1.0),
        max_u=umax, argmax_u=tmax,
        bc_lhs=lhs, bc_rhs=rhs, bc_residual=abs(lhs - rhs),
        min_u_over_Psi=float(ratio.min()), max_u_over_Psi=float(ratio.max()),
        mean_contraction_rate=float(np.exp(np.mean(np.log(rates[2:-3])))) if rates.size > 6 else None,
        increments=[float(x) for x in hist],
    )


def run_example(ex, tag):
    prob = Problem(ex["alpha"], ex["xi"], ex["w"])
    out = {}
    for n in (N_STD // 2, N_STD, 2 * N_STD):
        D, u, hist = solve_problem(prob, n=n, q=Q)
        out[str(n)] = summarize(prob, D, u, hist)
        if n == N_STD:
            np.savez(os.path.join(RESULTS, f"solution_{tag}.npz"), s=D.s, u=u,
                     Psi=prob.Psi(D.s), increments=np.array(hist))
    std, fine, coarse = out[str(N_STD)], out[str(2 * N_STD)], out[str(N_STD // 2)]
    for key in ("u_quarter", "u_three_quarters", "u_1", "max_u"):
        e1 = abs(std[key] - coarse[key])
        e2 = abs(fine[key] - std[key])
        out.setdefault("refinement", {})[key] = dict(
            coarse=coarse[key], standard=std[key], fine=fine[key],
            diff_coarse_standard=e1, diff_standard_fine=e2,
            observed_order=float(np.log2(e1 / e2)) if e2 > 0 else None)
    return prob, out


def table_solutions(alphas=(1.2, 1.4, 1.5, 1.6, 1.8, 2.0)):
    rows, tex = [], []
    for al in alphas:
        prob = Problem(al, EX1["xi"], EX1["w"])
        D, u, hist = solve_problem(prob, n=N_STD, q=Q)
        tmax, umax = locate_max(D, u)
        ok = prob.w_a[0] / abs(prob.w_a[1]) >= prob.theta_star()
        psi = float(prob.Psi(0.75))
        r = dict(alpha=al, criterion=bool(ok), fmPsi=F_MINUS * psi, u34=node_value(D, u, 0.75),
                 fpPsi=F_PLUS * psi, u14=node_value(D, u, 0.25), max_u=umax, argmax=tmax,
                 iterations=len(hist))
        rows.append(r)
        mark = "" if ok else "^\\dagger"
        tex.append(f"${al:.1f}{mark}$&${r['fmPsi']:.6f}$&${r['u34']:.6f}$&${r['fpPsi']:.6f}$&"
                   f"${r['u14']:.6f}$&${umax:.6f}$&${tmax:.4f}$\\\\")
        np.savez(os.path.join(RESULTS, f"solution_alpha_{al:.1f}.npz"), s=D.s, u=u, tmax=tmax, umax=umax)
    return rows, "\n".join(tex)


def table_example2(prob, tag="ex2", points=(0.05, None, 0.25, 0.5, 0.75, 1.0)):
    data = np.load(os.path.join(RESULTS, f"solution_{tag}.npz"))
    s, u = data["s"], data["u"]
    from fracbvp import Discretisation
    D = Discretisation(prob, s)
    phi = f_example(s, u)
    rows, tex = [], []
    for t in points:
        if t is None:
            t = prob.t_star
        val = float(D.evaluate(phi, [t])[0])
        psi = float(prob.Psi(t))
        rows.append(dict(t=t, fmPsi=F_MINUS * psi, u=val, fpPsi=F_PLUS * psi))
        label = "t^\\ast=%.4f" % t if abs(t - prob.t_star) < 1e-12 else "%g" % t
        tex.append(f"${label}$&${F_MINUS * psi:.6f}$&${val:.6f}$&${F_PLUS * psi:.6f}$\\\\")
    return rows, "\n".join(tex)


def main():
    res = {}
    p1, res["example1"] = run_example(EX1, "ex1")
    p2, res["example2"] = run_example(EX2, "ex2")
    rows, tex = table_solutions()
    res["table_solutions"] = rows
    save_text("table_solutions.tex", tex + "\n")
    rows2, tex2 = table_example2(p2)
    res["table_example2"] = rows2
    save_text("table_example2.tex", tex2 + "\n")
    # refinement table
    lines = []
    for n in (N_STD // 2, N_STD, 2 * N_STD):
        a = res["example1"][str(n)]
        b = res["example2"][str(n)]
        lines.append(f"${a['nodes']}$&${a['u_three_quarters']:.10f}$&${a['max_u']:.10f}$&"
                     f"${b['u_three_quarters']:.10f}$&${b['max_u']:.10f}$\\\\")
    save_text("table_refinement.tex", "\n".join(lines) + "\n")
    save_json("solutions.json", res)
    for key in ("example1", "example2"):
        std = res[key][str(N_STD)]
        print(f"--- {key}")
        for k, v in std.items():
            if k != "increments":
                print(f"  {k:24s} {v}")
        for k, v in res[key]["refinement"].items():
            print(f"  refine {k:16s} {v}")
    print(tex)
    print(tex2)
    print("\n".join(lines))


if __name__ == "__main__":
    main()
