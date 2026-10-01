"""All figures of the paper (Figures 1-8, Section 7).

Run after solve_examples.py (Figures 5, 6 and 8 read the stored solutions).
Output: figures/fig*.pdf
"""
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle

from _setup import FIGURES, RESULTS
from fracbvp import Problem, theta_star, theta_star_r, F_MINUS, F_PLUS
from compute_constants import EX1, EX2

# ----------------------------------------------------------------------------
# style
# ----------------------------------------------------------------------------
TEAL, TEAL2 = "#0E6E74", "#0F7C82"
DEEP = "#123B54"
CORAL, CORAL2 = "#C13B4E", "#D1495B"
MOSS = "#2F7D5B"
AMB = "#C8801A"
PLUM = "#7A3E9D"
GREY = "#8A99A8"
EDGE = "#5B6B7B"

plt.rcParams.update({
    "pdf.fonttype": 42, "ps.fonttype": 42,
    "font.family": "serif", "font.serif": ["DejaVu Serif"],
    "mathtext.fontset": "dejavusans",
    "font.size": 9.5, "axes.titlesize": 10, "axes.labelsize": 10,
    "axes.edgecolor": EDGE, "axes.labelcolor": "#1d2b36",
    "xtick.color": EDGE, "ytick.color": EDGE,
    "axes.grid": True, "grid.color": "#dde4ea", "grid.linewidth": 0.7,
    "axes.spines.top": False, "axes.spines.right": False,
    "legend.frameon": True, "legend.framealpha": 0.92, "legend.edgecolor": "#c9d3dc",
    "savefig.bbox": "tight", "savefig.pad_inches": 0.03,
})
FIELD = LinearSegmentedColormap.from_list(
    "field", ["#ffffff", "#e6f1f1", "#a9d0cf", "#4f9d9d", TEAL, DEEP, "#231f3f"])
FIELD2 = LinearSegmentedColormap.from_list(
    "field2", ["#ffffff", "#d9ebea", "#8cc2c0", "#3f8f90", TEAL])
# ordered family of curves indexed by alpha: lightness decreases with alpha,
# so the order is readable without colour; every curve is also labelled
ALPHA_CMAP = LinearSegmentedColormap.from_list(
    "alpha", ["#C27612", "#C13B4E", "#7A3E9D", "#0B5E66", "#10283A"])


def alpha_color(al):
    return ALPHA_CMAP(float(np.clip((al - 1.0) / 1.0, 0.0, 1.0)))

P1 = Problem(EX1["alpha"], EX1["xi"], EX1["w"])
P2 = Problem(EX2["alpha"], EX2["xi"], EX2["w"])


def save(fig, name):
    path = os.path.join(FIGURES, name)
    fig.savefig(path)
    plt.close(fig)
    return path


# ----------------------------------------------------------------------------
def fig1_channels():
    p = P1
    t = np.linspace(0, 1, 601)
    T, S = np.meshgrid(t, t, indexing="xy")          # x: t, y: s
    fields = [p.G(T, S), T ** p.nu * p.KA(S) / p.Delta, p.calG(T, S)]
    titles = [r"local Green function $G(t,s)$",
              r"nonlocal part $\Delta^{-1}t^{\nu}\mathcal{K}_A(s)$",
              r"perturbed kernel $\mathcal{G}(t,s)$"]
    vmax = float(fields[2].max())
    fig, axs = plt.subplots(1, 3, figsize=(7.9, 2.75), sharey=True)
    for k, (ax, F, ti) in enumerate(zip(axs, fields, titles)):
        im = ax.imshow(F, origin="lower", extent=(0, 1, 0, 1), cmap=FIELD, vmin=0, vmax=vmax,
                       aspect="auto", interpolation="bilinear", rasterized=True)
        ax.contour(t, t, F, levels=np.linspace(0.1, vmax, 11), colors="white",
                   linewidths=0.6, alpha=0.9)
        ax.set_title(ti, fontsize=9.5)
        ax.set_xlabel(r"$t$")
        ax.grid(False)
        ax.set_xticks([0, 0.5, 1])
        if k == 0:
            ax.plot([0, 1], [0, 1], "--", color=CORAL2, lw=1.4)
            ax.text(0.62, 0.70, r"ridge $t=s$", rotation=45, color=CORAL2, fontsize=8)
            ax.set_ylabel(r"$s$")
        else:
            for xv, lab in ((0.75, r"$\xi_1$"), (0.25, r"$\xi_2$")):
                ax.axhline(xv, color=AMB, ls=":", lw=1.3)
                ax.text(0.02, xv + 0.025, lab, color=AMB, fontsize=9)
    cb = fig.colorbar(im, ax=axs, fraction=0.035, pad=0.02)
    cb.set_label("kernel value")
    cb.outline.set_edgecolor("#c9d3dc")
    return save(fig, "fig1_channels.pdf")


# ----------------------------------------------------------------------------
def fig2_positivity():
    p = P1
    fig, axs = plt.subplots(1, 3, figsize=(9.6, 3.15))
    # (a) R(s)
    ax = axs[0]
    s = np.linspace(1e-4, 1 - 1e-4, 4000)
    R = p.G(0.25, s) / p.G(0.75, s)
    th = p.theta_star()
    ymax = 4.25
    ax.fill_between(s, th, ymax, color=TEAL, alpha=0.13, lw=0)
    ax.text(0.5, 3.93, r"admissible ratios $w_1/|w_2|\geq\theta^\ast$", color=TEAL,
            fontsize=8.3, ha="center")
    ax.plot(s, R, color=TEAL, lw=2.4)
    ax.axhline(th, color=CORAL2, ls=(0, (6, 3)), lw=1.4)
    ax.axhline(10 / 3, color="#1d2b36", lw=1.3)
    ax.text(0.97, 10 / 3 + 0.05, r"$w_1/|w_2|=10/3$", ha="right", va="bottom", fontsize=8.3)
    ax.plot([0.25], [th], "o", color=CORAL2, ms=6, mec="white", mew=1.0, zorder=5)
    ax.axvline(0.25, ymax=(th - 0.35) / (ymax - 0.35), color=CORAL2, ls=":", lw=1.1)
    ax.annotate(r"$\theta^\ast=\sqrt{3}+\sqrt{2}$", xy=(0.25, th), xytext=(0.38, 2.75),
                color=CORAL2, fontsize=8.5,
                arrowprops=dict(arrowstyle="-", color=CORAL2, lw=0.8),
                bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.9))
    ax.text(0.55, 1.35, r"$R(s)$", color=TEAL, fontsize=9.5)
    ax.set_ylim(0.35, ymax)
    ax.set_xlim(0, 1)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1])
    ax.set_xticklabels(["0", r"$\xi_2$", "0.5", r"$\xi_1$", "1"])
    ax.set_xlabel(r"$s$")
    ax.set_ylabel(r"$R(s)=G(\xi_2,s)/G(\xi_1,s)$")
    ax.set_title(r"(a)  criterion at $\alpha=3/2$")
    # (b) threshold surface
    ax = axs[1]
    nu = np.linspace(0.05, 1.0, 400)
    r = np.linspace(0.05, 0.95, 400)
    NU, RR = np.meshgrid(nu, r)
    TH = theta_star_r(RR, NU)
    im = ax.imshow(np.minimum(TH, 12), origin="lower", extent=(nu[0], nu[-1], r[0], r[-1]),
                   aspect="auto", cmap=FIELD2, vmin=1, vmax=12, interpolation="bilinear",
                   rasterized=True)
    cs = ax.contour(NU, RR, TH, levels=[2, 5, 8], colors="#34495e", linewidths=0.7)
    ax.clabel(cs, fmt="%g", fontsize=7.5)
    c2 = ax.contour(NU, RR, TH, levels=[10 / 3], colors=CORAL2, linewidths=2.2)
    ax.clabel(c2, fmt={10 / 3: r"$\theta^\ast=10/3$"}, fontsize=8, manual=[(0.33, 0.62)])
    ax.plot([0.5, 0.5], [0.05, 1 / 3], ":", color="#1d2b36", lw=0.9)
    ax.plot([0.05, 0.5], [1 / 3, 1 / 3], ":", color="#1d2b36", lw=0.9)
    ax.plot([0.5], [1 / 3], "o", color="#1d2b36", ms=6.5, mec="white", mew=1.2, zorder=5)
    ax.annotate(r"$(\nu,r)=(1/2,\,1/3)$", xy=(0.5, 1 / 3), xytext=(0.14, 0.19), fontsize=8.5,
                arrowprops=dict(arrowstyle="-", color="#1d2b36", lw=0.8),
                bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.9))
    ax.set_xlabel(r"$\nu=\alpha-1$")
    ax.set_ylabel(r"$r=\xi_2/\xi_1$")
    ax.set_title(r"(b)  threshold surface $\theta^\ast(\nu,r)$")
    ax.grid(False)
    cb = fig.colorbar(im, ax=ax, extend="max", fraction=0.05, pad=0.03)
    cb.set_label(r"$\theta^\ast$")
    cb.set_ticks([1, 3, 6, 9, 12])
    cb.outline.set_edgecolor("#c9d3dc")
    # (c) K_A / Delta
    ax = axs[2]
    s = np.linspace(0, 1, 4001)
    K = p.KA(s) / p.Delta
    ax.fill_between(s, 0, K, color="#9fb2c2", alpha=0.35, lw=0)
    ax.plot(s, K, color=DEEP, lw=2.4)
    ax.axhline(0, color=EDGE, lw=0.8)
    ax.axvline(0.25, color=CORAL2, ls=":", lw=1.1)
    ax.axvline(0.75, color=MOSS, ls=":", lw=1.1)
    ax.text(0.23, 1.02, r"$-|w_2|$ at $\xi_2$", color=CORAL2, ha="right", fontsize=8.5)
    ax.text(0.73, 1.02, r"$+w_1$ at $\xi_1$", color=MOSS, ha="right", fontsize=8.5)
    ax.plot([0, 1], [K[0], K[-1]], "o", color="white", mec=DEEP, mew=1.3, ms=5, zorder=5)
    ax.annotate(r"$\mathcal{K}_A(0)=0$", xy=(0.0, 0.0), xytext=(0.08, 0.12), fontsize=8.5,
                color=DEEP, arrowprops=dict(arrowstyle="-", color=DEEP, lw=0.8))
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.05, 1.25)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1])
    ax.set_xticklabels(["0", r"$\xi_2$", "0.5", r"$\xi_1$", "1"])
    ax.set_xlabel(r"$s$")
    ax.set_ylabel(r"$\Delta^{-1}\mathcal{K}_A(s)$")
    ax.set_title(r"(c)  $\mathcal{K}_A\geq0$ for signed $\mathrm{d}A$")
    fig.tight_layout(w_pad=1.2)
    return save(fig, "fig2_positivity.pdf")


# ----------------------------------------------------------------------------
def fig3_psi_shift():
    p = P1
    t = np.linspace(0, 1, 1200)
    Psi, Psi0 = p.Psi(t), p.Psi0(t)
    ts, t0 = p.t_star, p.nu
    fig, ax = plt.subplots(figsize=(4.5, 3.2))
    ax.fill_between(t, Psi0, Psi, color="#F2DDB8", alpha=0.65, lw=0, label="nonlocal contribution")
    ax.plot(t, Psi0, ls=(0, (6, 3)), color=GREY, lw=2.0, label=r"$\mathrm{d}A=0$")
    ax.plot(t, Psi, color=TEAL, lw=2.6, label=r"$\mathrm{d}A$ sign-changing")
    y0, ys = float(p.Psi0(t0)), float(p.Psi(ts))
    ax.plot([t0], [y0], "o", color=GREY, mec="white", ms=6, zorder=5)
    ax.plot([ts], [ys], "o", color=TEAL, mec="white", ms=6, zorder=5)
    ya = 1.045
    ax.annotate("", xy=(ts, ya), xytext=(t0, ya),
                arrowprops=dict(arrowstyle="-|>", color=CORAL2, lw=2.2, mutation_scale=14))
    ax.text((t0 + ts) / 2, ya + 0.06, "maximiser shifts right", color=CORAL2, ha="center",
            fontsize=8.5)
    ax.plot([t0, t0], [y0, ya], ":", color=GREY, lw=0.9)
    ax.plot([ts, ts], [ys, ya], ":", color=TEAL, lw=1.1)
    ax.annotate(r"$\alpha-1=0.5$", xy=(t0, y0), xytext=(0.31, 0.43), color=GREY, fontsize=8.5,
                arrowprops=dict(arrowstyle="-", color=GREY, lw=0.8))
    ax.annotate(r"$t^\ast=%.6f$" % ts, xy=(ts, ys), xytext=(0.75, 0.75), color=TEAL, fontsize=8.5,
                arrowprops=dict(arrowstyle="-", color=TEAL, lw=0.8),
                bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.9))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1.2)
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$\Psi(t)$")
    ax.legend(loc="lower center", fontsize=8.5)
    return save(fig, "fig3_psi_shift.pdf")


# ----------------------------------------------------------------------------
def fig4_harnack():
    fig, ax = plt.subplots(figsize=(4.5, 3.3))
    s = np.logspace(-5, np.log10(0.4), 400)
    t0 = 0.75
    cols = {al: alpha_color(al) for al in (1.08, 1.2, 1.4, 1.6, 1.8, 2.0)}
    for al, c in cols.items():
        nu = al - 1
        ratio = s ** nu / (t0 ** nu - np.maximum(t0 - s, 0) ** nu)
        ax.loglog(s, ratio, color=c, lw=2.1)
        ax.text(8.5e-6, ratio[0], r"$\alpha=%g$" % al, color=c, fontsize=8.5, ha="right",
                va="center")
    # slope triangle on the top curve (alpha = 1.08), where there is free space above it
    nu = 0.08
    f = lambda x: x ** nu / (t0 ** nu - (t0 - x) ** nu)
    x1, x2 = 1e-3, 1e-2
    ax.plot([x1, x2], [f(x1), f(x1)], color="#1d2b36", lw=1.1)
    ax.plot([x2, x2], [f(x1), f(x2)], color="#1d2b36", lw=1.1)
    ax.text(np.sqrt(x1 * x2), f(x1) * 1.35, "1 decade", ha="center", fontsize=8)
    ax.text(x2 * 1.12, np.sqrt(f(x1) * f(x2)), r"$10^{\nu-1}$", ha="left", va="center", fontsize=8)
    ax.text(1.1e-4, f(1e-4) * 3.2, r"slope $\nu-1$", fontsize=8.5)
    ax.set_xlim(4e-7, 1)
    ax.set_ylim(0.5, 3e6)
    ax.set_xlabel(r"$s$")
    ax.set_ylabel(r"$\sup_t\,G(t,s)\,/\,G(3/4,s)$")
    return save(fig, "fig4_harnack.pdf")


# ----------------------------------------------------------------------------
def _band_plot(ax, p, sol, b, c, A0, A1, t2, lab_value):
    s, u = sol["s"], sol["u"]
    t = np.linspace(0, 1, 1500)
    Psi = p.Psi(t)
    sig = F_MINUS / F_PLUS
    lo = sig * Psi / float(p.Psi(t2)) * b
    hi = Psi / float(p.Psi(t2)) * b / sig
    ax.fill_between(t, lo, hi, color=TEAL, alpha=0.13, lw=0)
    ax.plot(t, lo, color=TEAL, lw=0.9, alpha=0.8)
    ax.plot(t, hi, color=TEAL, lw=0.9, alpha=0.8)
    ax.axhspan(c, b, color=TEAL, alpha=0.05, lw=0)
    ax.axhline(b, color=CORAL2, ls=":", lw=1.1)
    ax.axhline(c, color=MOSS, ls=":", lw=1.1)
    ax.plot(t, A0 * Psi, ls=(0, (7, 3)), color=AMB, lw=1.9, label=r"$x_0=%g\,\Psi$" % A0)
    ax.plot(t, A1 * Psi, ls="-.", color=PLUM, lw=1.9,
            label=r"$x_1=\Psi$" if A1 == 1 else r"$x_1=%g\,\Psi$" % A1)
    ax.plot(s, u, color="#15202b", lw=2.4, label=r"solution $u^\ast$")
    k = int(np.argmin(np.abs(s - t2)))
    ax.plot([t2, t2], [0, u[k]], ":", color=CORAL2, lw=1.0)
    ax.plot([0, t2], [u[k], u[k]], ":", color=CORAL2, lw=0.9)
    ax.plot([t2], [u[k]], "o", color=CORAL2, mec="white", ms=7, zorder=6)
    return u[k]


def fig5_localisation():
    p = P1
    sol = np.load(os.path.join(RESULTS, "solution_ex1.npz"))
    fig, ax = plt.subplots(figsize=(5.1, 3.5))
    val = _band_plot(ax, p, sol, EX1["b"], EX1["c"], EX1["A0"], EX1["A1"], EX1["t2"], None)
    ax.text(0.015, EX1["b"] + 0.03, r"$b=1.2$", color=CORAL2, fontsize=8.5)
    ax.text(0.015, EX1["c"] - 0.075, r"$c=0.8$", color=MOSS, fontsize=8.5)
    ax.text(0.30, 1.39, "admissible band", color=TEAL, fontsize=8.5)
    ax.annotate(r"$u^\ast(3/4)=%.6f$" % val, xy=(0.75, val), xytext=(0.81, 0.58), color=CORAL2,
                fontsize=8.5, arrowprops=dict(arrowstyle="-", color=CORAL2, lw=0.8),
                bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.9))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1.58)
    ax.set_xticks([0, 0.2, 0.4, 0.6, 0.75, 1.0])
    ax.set_xticklabels(["0.0", "0.2", "0.4", "0.6", "3/4", "1.0"])
    ax.set_xlabel(r"$t$")
    ax.set_ylabel("value")
    ax.legend(loc="lower right", fontsize=8.5)
    return save(fig, "fig5_localisation.pdf")


# ----------------------------------------------------------------------------
def fig6_solutions():
    fig, axs = plt.subplots(1, 2, figsize=(7.6, 2.9))
    ax = axs[0]
    cols = {al: alpha_color(al) for al in (1.2, 1.4, 1.5, 1.6, 1.8, 2.0)}
    offs = {1.2: (0.0, 0.035), 1.4: (0.0, 0.035), 1.5: (-0.01, 0.035), 1.6: (-0.03, 0.035),
            1.8: (-0.07, 0.03), 2.0: (-0.075, -0.075)}
    for al, c in cols.items():
        d = np.load(os.path.join(RESULTS, f"solution_alpha_{al:.1f}.npz"))
        ax.plot(d["s"], d["u"], color=c, lw=2.0, label=r"$\alpha=%.1f$" % al)
        tm, um = float(d["tmax"]), float(d["umax"])
        ax.plot([tm], [um], "o", color=c, ms=5, mec="white", mew=0.8, zorder=5)

    ax.set_xlim(0, 1)
    ax.set_ylim(-0.02, 1.06)
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$u^\ast(t)$")
    ax.set_title("solution profiles")
    ax.legend(ncol=2, fontsize=8, loc="lower right", frameon=False)
    ax = axs[1]
    inc = np.load(os.path.join(RESULTS, "solution_ex1.npz"))["increments"]
    k = np.arange(1, inc.size + 1)
    ax.semilogy(k, inc, "-o", color=TEAL, ms=3.5, mfc="white", lw=1.8)
    ax.fill_between(k, 1e-16, inc, color=TEAL, alpha=0.1, lw=0)
    ax.set_ylim(3e-17, 30)
    ax.set_xlabel(r"iteration $k$")
    ax.set_ylabel(r"$\|u_{k+1}-u_k\|_\infty$")
    ax.set_title(r"convergence, $\alpha=3/2$")
    fig.tight_layout(w_pad=2.0)
    return save(fig, "fig6_solutions.pdf")


# ----------------------------------------------------------------------------
def fig7_failureset():
    fig, axs = plt.subplots(1, 2, figsize=(7.4, 3.4), sharey=False)
    v = 1.0
    for ax, d2 in zip(axs, (0.5, -0.5)):
        u = np.linspace(-3, 3, 1201)
        d1 = np.linspace(-3, 3, 1201)
        U, D1 = np.meshgrid(u, d1)
        F1 = (U * D1 >= 0) & (np.abs(U) >= v) & (np.abs(D1) >= abs(d2))
        F2 = (d2 >= 0) & (v >= np.abs(U)) & (abs(d2) >= np.abs(D1))
        ax.contourf(U, D1, F1.astype(float), levels=[0.5, 1.5], colors=[CORAL2], alpha=0.3)
        if d2 >= 0:
            ax.contourf(U, D1, F2.astype(float), levels=[0.5, 1.5], colors=[PLUM], alpha=0.3)
        ax.contour(U, D1, (F1 | F2).astype(float), levels=[0.5], colors="#15202b", linewidths=1.3)
        ax.axvline(v, color=TEAL, ls=":", lw=1.0)
        ax.axvline(-v, color=TEAL, ls=":", lw=1.0)
        ax.axhline(abs(d2), color=AMB, ls=":", lw=1.0)
        ax.axhline(-abs(d2), color=AMB, ls=":", lw=1.0)
        ax.axhline(0, color=GREY, lw=0.6)
        ax.axvline(0, color=GREY, lw=0.6)
        ax.set_xlim(-3, 3)
        ax.set_ylim(-3, 3)
        ax.set_xlabel(r"$u$")
        ax.set_title(r"$v=1,\ \ d_2=%g$" % d2)
        ax.grid(False)
        if d2 > 0:
            ax.set_ylabel(r"$d_1$")
            ax.plot([2], [2], "s", color=CORAL2, mec="white", ms=8, zorder=6)
            ax.text(1.85, 2.35, r"$(2,2)\in$ (F1)", color=CORAL2, fontsize=8.5, ha="right")
            ax.plot([0.2], [-2], "o", color=TEAL, ms=6, zorder=6)
            ax.text(0.45, -2.05, r"$(0.2,-2)$: inequality holds", color=TEAL, fontsize=8.5,
                    va="center")
    from matplotlib.patches import Patch
    handles = [Patch(color=CORAL2, alpha=0.3, label="(F1)"), Patch(color=PLUM, alpha=0.3, label="(F2)"),
               Patch(fc="white", ec="#15202b", label=r"$L<\mathsf{R}$ holds outside")]
    axs[1].legend(handles=handles, loc="lower right", fontsize=8.5)
    fig.tight_layout(w_pad=2.5)
    return save(fig, "fig7_failureset.pdf")


# ----------------------------------------------------------------------------
def fig8_example2():
    fig, axs = plt.subplots(1, 2, figsize=(8.2, 3.2))
    # (a) K_A / Delta at alpha = 1.08 for the two weight choices
    ax = axs[0]
    s = np.linspace(0, 1, 4001)
    pa = Problem(EX2["alpha"], EX2["xi"], EX1["w"])     # weights of Example 1
    Ka = pa.KA(s) / pa.Delta
    Kb = P2.KA(s) / P2.Delta
    ax.fill_between(s, Ka, 0, where=Ka < 0, color=CORAL2, alpha=0.25, lw=0)
    ax.plot(s, Ka, color=CORAL2, lw=2.0, label=r"$w_1/|w_2|=10/3<\theta^\ast$")
    ax.plot(s, Kb, color=DEEP, lw=2.4, label=r"$w_1/|w_2|=30>\theta^\ast$")
    ax.axhline(0, color=EDGE, lw=0.8)
    ax.axvline(0.25, color=CORAL2, ls=":", lw=1.0)
    ax.axvline(0.75, color=MOSS, ls=":", lw=1.0)
    ax.text(0.53, 0.52, r"$\theta^\ast=%.3f$ at $\alpha=1.08$" % P2.theta_star(), ha="center",
            fontsize=8.5, color="#15202b",
            bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="#c9d3dc"))
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.4, 1.7)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1])
    ax.set_xticklabels(["0", r"$\xi_2$", "0.5", r"$\xi_1$", "1"])
    ax.set_xlabel(r"$s$")
    ax.set_ylabel(r"$\Delta^{-1}\mathcal{K}_A(s)$")
    ax.set_title(r"(a)  sign of $\mathcal{K}_A$ at $\alpha=1.08$")
    ax.legend(loc="upper left", fontsize=8)
    # (b) solution with envelope
    ax = axs[1]
    sol = np.load(os.path.join(RESULTS, "solution_ex2.npz"))
    t = np.linspace(0, 1, 2000)
    Psi = P2.Psi(t)
    ax.fill_between(t, F_MINUS * Psi, F_PLUS * Psi, color=TEAL, alpha=0.15, lw=0,
                    label=r"$[\,f_-\Psi,\ f_+\Psi\,]$")
    ax.plot(t, F_MINUS * Psi, color=TEAL, lw=0.9)
    ax.plot(t, F_PLUS * Psi, color=TEAL, lw=0.9)
    ax.plot(sol["s"], sol["u"], color="#15202b", lw=2.3, label=r"solution $u^\ast$")
    ax.plot(t, EX2["A0"] * Psi, ls=(0, (7, 3)), color=AMB, lw=1.7, label=r"$x_0=0.9\,\Psi$")
    ax.axhline(EX2["b"], color=CORAL2, ls=":", lw=1.1)
    ax.axhline(EX2["c"], color=MOSS, ls=":", lw=1.1)
    ax.text(0.62, EX2["b"] + 0.025, r"$b=1.1$", color=CORAL2, fontsize=8.5)
    ax.text(0.30, EX2["c"] - 0.075, r"$c=0.7$", color=MOSS, fontsize=8.5)
    s2, u2 = sol["s"], sol["u"]
    k = int(np.argmin(np.abs(s2 - 0.75)))
    ax.plot([0.75, 0.75], [0, u2[k]], ":", color=CORAL2, lw=1.0)
    ax.plot([0.75], [u2[k]], "o", color=CORAL2, mec="white", ms=7, zorder=6)
    ax.annotate(r"$u^\ast(3/4)=%.6f$" % u2[k], xy=(0.75, u2[k]), xytext=(0.72, 0.30),
                color=CORAL2, fontsize=8.5, arrowprops=dict(arrowstyle="-", color=CORAL2, lw=0.8),
                bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.9))
    ts = P2.t_star
    ax.axvline(ts, color=GREY, ls=":", lw=0.9)
    ax.text(ts + 0.015, 0.45, r"$t^\ast=%.4f$" % ts, color=GREY, fontsize=8.5)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1.5)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.set_xticklabels(["0", "1/4", "1/2", "3/4", "1"])
    ax.set_xlabel(r"$t$")
    ax.set_ylabel("value")
    ax.set_title(r"(b)  solution for $\alpha=1.08$")
    ax.legend(loc="lower center", bbox_to_anchor=(0.36, 0.02), fontsize=8)
    fig.tight_layout(w_pad=2.0)
    return save(fig, "fig8_example2.pdf")


def main():
    for f in (fig1_channels, fig2_positivity, fig3_psi_shift, fig4_harnack, fig5_localisation,
              fig6_solutions, fig7_failureset, fig8_example2):
        print(f())


if __name__ == "__main__":
    main()
