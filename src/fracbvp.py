"""
Closed-form quantities and a product-integration solver for the problem

    D^alpha u(t) + f(t, u(t)) = 0,   0 < t < 1,
    u(0) = 0,   D^{nu} u(1) = int_0^1 u(tau) dA(tau),     nu = alpha - 1,

with 1 < alpha <= 2 and a finite atomic Stieltjes measure
dA = sum_k w_k delta_{xi_k}.  The weights w_k may have either sign.

Notation follows the paper:
    G(t,s)      local Green function (the case dA = 0), eq. (7)
    K_A(s)      = int G(tau, s) dA(tau)
    mu, Delta   mu = int tau^nu dA,  Delta = Gamma(alpha) - mu
    calG(t,s)   perturbed kernel G(t,s) + t^nu K_A(s) / Delta   (Theorem 2, eq. (14))
    Psi(t)      = int_0^1 calG(t,s) ds                            (Lemma 9, eq. (34))
    Kbar        = Delta^{-1} int_0^1 K_A(s) ds
    t_star      maximiser of Psi                                  (eq. (36))
    theta_star  weight-ratio threshold of Theorem 3               (eq. (25))
"""
from __future__ import annotations

from dataclasses import dataclass
from math import gamma

import numpy as np


# ----------------------------------------------------------------------------
# closed-form quantities
# ----------------------------------------------------------------------------
def theta_star(xi1: float, xi2: float, nu: float) -> float:
    """Threshold theta* = xi2^nu / (xi1^nu - (xi1 - xi2)^nu), 0 < xi2 < xi1 <= 1."""
    return xi2 ** nu / (xi1 ** nu - (xi1 - xi2) ** nu)


def theta_star_r(r: np.ndarray | float, nu: np.ndarray | float):
    """Same threshold written with r = xi2/xi1 (Corollary 1): r^nu / (1 - (1-r)^nu)."""
    r = np.asarray(r, float)
    nu = np.asarray(nu, float)
    return r ** nu / (1.0 - (1.0 - r) ** nu)


@dataclass
class Problem:
    alpha: float
    xi: tuple
    w: tuple

    def __post_init__(self):
        if not (1.0 < self.alpha <= 2.0):
            raise ValueError("need 1 < alpha <= 2")
        self.nu = self.alpha - 1.0
        self.Ga = gamma(self.alpha)
        self.xi_a = np.asarray(self.xi, float)
        self.w_a = np.asarray(self.w, float)
        self.mu = float(np.sum(self.w_a * self.xi_a ** self.nu))
        self.Delta = self.Ga - self.mu

    # local Green function -------------------------------------------------
    def G(self, t, s):
        t = np.asarray(t, float)
        s = np.asarray(s, float)
        t, s = np.broadcast_arrays(t, s)
        d = np.clip(t - s, 0.0, None)
        return (t ** self.nu - d ** self.nu) / self.Ga

    def KA(self, s):
        s = np.asarray(s, float)
        out = np.zeros_like(s, dtype=float)
        for wk, xk in zip(self.w_a, self.xi_a):
            out = out + wk * self.G(xk, s)
        return out

    def calG(self, t, s):
        t = np.asarray(t, float)
        s = np.asarray(s, float)
        t, s = np.broadcast_arrays(t, s)
        return self.G(t, s) + t ** self.nu * self.KA(s) / self.Delta

    # profile Psi ----------------------------------------------------------
    @property
    def Kbar(self) -> float:
        a, nu = self.alpha, self.nu
        s = np.sum(self.w_a * (self.xi_a ** nu - self.xi_a ** a / a))
        return float(s / (self.Delta * self.Ga))

    def Psi(self, t):
        t = np.asarray(t, float)
        return t ** self.nu * ((1.0 - t / self.alpha) / self.Ga + self.Kbar)

    def Psi0(self, t):
        """Profile for dA = 0."""
        t = np.asarray(t, float)
        return t ** self.nu * (1.0 - t / self.alpha) / self.Ga

    @property
    def t_star(self) -> float:
        return min(1.0, self.nu * (1.0 + self.Ga * self.Kbar))

    @property
    def m_Psi(self) -> float:
        return (1.0 - 1.0 / self.alpha) / self.Ga + self.Kbar

    @property
    def M_Psi(self) -> float:
        return 1.0 / self.Ga + self.Kbar

    def theta_star(self) -> float:
        """Threshold for a two-atom measure with w2 < 0 < w1 (xi1 > xi2)."""
        i1 = int(np.argmax(self.w_a))
        i2 = int(np.argmin(self.w_a))
        return theta_star(self.xi_a[i1], self.xi_a[i2], self.nu)


# ----------------------------------------------------------------------------
# product integration
# ----------------------------------------------------------------------------
_GL_Z, _GL_W = np.polynomial.legendre.leggauss(12)
_GL_Z = 0.5 * (_GL_Z + 1.0)      # nodes on [0, 1]
_GL_W = 0.5 * _GL_W              # weights on [0, 1]


def make_mesh(n: int, q: float = 2.0, extra=()) -> np.ndarray:
    """Graded mesh s_j = (j/n)^q on [0,1], with the points in `extra` added."""
    s = (np.arange(n + 1) / n) ** q
    s = np.concatenate([s, np.asarray(extra, float)])
    s = np.unique(np.round(s, 15))
    s[0], s[-1] = 0.0, 1.0
    return s


def _interval_weights(a, h, nu):
    """Weights of int_a^{a+h} r^nu (r-a)/h dr and int_a^{a+h} r^nu (a+h-r)/h dr.

    Returned as (WL, WR): WL multiplies the value at the interval end that is
    farther from the target (distance b = a+h), WR the nearer one (distance a).
    Closed forms are used when a < 2h, Gauss-Legendre otherwise (no cancellation).
    """
    a = np.asarray(a, float)
    h = np.asarray(h, float)
    WL = np.empty_like(a)
    WR = np.empty_like(a)
    near = a < 2.0 * h
    # closed forms (stable when a is comparable with h)
    an, hn = a[near], h[near]
    bn = an + hn
    p1 = (bn ** (nu + 1) - an ** (nu + 1)) / (nu + 1)
    p2 = (bn ** (nu + 2) - an ** (nu + 2)) / (nu + 2)
    WL[near] = (p2 - an * p1) / hn          # coefficient (r - a)/h
    WR[near] = (bn * p1 - p2) / hn          # coefficient (b - r)/h
    # Gauss-Legendre for well separated intervals
    af, hf = a[~near], h[~near]
    accL = np.zeros_like(af)
    accR = np.zeros_like(af)
    for z, wq in zip(_GL_Z, _GL_W):
        val = (af + hf * z) ** nu
        accL += wq * val * z
        accR += wq * val * (1.0 - z)
    WL[~near] = hf * accL
    WR[~near] = hf * accR
    return WL, WR


def riemann_liouville_weights(s: np.ndarray, x: np.ndarray, nu: float, chunk=400):
    """Matrix W with  int_0^x (x - r)^nu phi(r) dr = W @ phi  exactly for every
    phi that is piecewise linear on the mesh s.  Targets x may be arbitrary points
    of [0,1]; a target inside an interval uses the partial interval [s_j, x].
    """
    s = np.asarray(s, float)
    x = np.atleast_1d(np.asarray(x, float))
    n = s.size - 1
    left, right = s[:-1], s[1:]
    h = right - left
    W = np.zeros((x.size, n + 1))
    for i0 in range(0, x.size, chunk):
        xs = x[i0:i0 + chunk]
        # full intervals: right <= x
        a = xs[:, None] - right[None, :]            # distance to right end
        full = a >= -1e-15
        a = np.where(full, np.maximum(a, 0.0), 0.0)
        hh = np.broadcast_to(h, a.shape)
        WL, WR = _interval_weights(a.ravel(), hh.ravel(), nu)
        WL = WL.reshape(a.shape) * full
        WR = WR.reshape(a.shape) * full
        # variable r = x - s: the right node s_{j+1} is the near end (distance a)
        # phi = phi_j (s_{j+1}-s)/h + phi_{j+1} (s-s_j)/h, and with r = x - s,
        # (s_{j+1}-s)/h = (r - a)/h  -> left node gets WL,
        # (s - s_j)/h   = (b - r)/h  -> right node gets WR.
        block = np.zeros((xs.size, n + 1))
        block[:, :-1] += WL
        block[:, 1:] += WR
        # partial interval containing x (if x is not a node)
        for k, xv in enumerate(xs):
            j = np.searchsorted(s, xv, side="right") - 1
            if j < 0 or j >= n or abs(xv - s[j]) < 1e-15:
                continue
            hj = s[j + 1] - s[j]
            L = xv - s[j]                          # partial length
            # int_{s_j}^{x} (x-s)^nu [phi_j (s_{j+1}-s)/hj + phi_{j+1} (s-s_j)/hj] ds
            # with r = x - s in [0, L]:  s - s_j = L - r,  s_{j+1} - s = hj - L + r
            i_r0 = L ** (nu + 1) / (nu + 1)
            i_r1 = L ** (nu + 2) / (nu + 2)
            block[k, j] += ((hj - L) * i_r0 + i_r1) / hj
            block[k, j + 1] += (L * i_r0 - i_r1) / hj
        W[i0:i0 + chunk] = block
    return W


def trapezoid_weights(s: np.ndarray) -> np.ndarray:
    h = np.diff(s)
    m = np.zeros_like(s)
    m[:-1] += h / 2
    m[1:] += h / 2
    return m


class Discretisation:
    """Nystrom/product-integration discretisation of (T u)(t) = int calG(t,s) f(s,u(s)) ds.

    phi = f(s,u(s)) is replaced by its piecewise linear interpolant on the mesh,
    and every integral of the kernel against that interpolant is computed exactly.
    """

    def __init__(self, prob: Problem, s: np.ndarray):
        self.p = prob
        self.s = np.asarray(s, float)
        nu, Ga = prob.nu, prob.Ga
        self.M = trapezoid_weights(self.s)                      # int phi
        self.W = riemann_liouville_weights(self.s, self.s, nu)   # int_0^{s_i}(s_i-r)^nu phi
        self.Wxi = riemann_liouville_weights(self.s, prob.xi_a, nu)
        # coefficient vector of t^nu in u:  c1 = M/Ga + C,  C = sum_k w_k(xi_k^nu M - Wxi_k)/(Delta Ga)
        Cvec = np.zeros_like(self.M)
        for wk, xk, row in zip(prob.w_a, prob.xi_a, self.Wxi):
            Cvec += wk * (xk ** nu * self.M - row)
        self.Cvec = Cvec / (prob.Delta * Ga)
        self.c1vec = self.M / Ga + self.Cvec
        tn = self.s ** nu
        self.B = (np.outer(tn, self.M) - self.W) / Ga + np.outer(tn, self.Cvec)

    def apply(self, phi: np.ndarray) -> np.ndarray:
        return self.B @ phi

    def evaluate(self, phi: np.ndarray, t) -> np.ndarray:
        """Value of int calG(t,s) phi_h(s) ds at arbitrary t (phi_h the interpolant)."""
        t = np.atleast_1d(np.asarray(t, float))
        Wt = riemann_liouville_weights(self.s, t, self.p.nu)
        return (t ** self.p.nu * (self.c1vec @ phi)) - (Wt @ phi) / self.p.Ga

    def solve(self, f, tol=1e-15, maxit=500, u0=None):
        """Picard iteration u_{k+1} = T_h u_k from u0 (default 0)."""
        u = np.zeros_like(self.s) if u0 is None else np.array(u0, float)
        hist = []
        for k in range(1, maxit + 1):
            un = self.apply(f(self.s, u))
            d = float(np.max(np.abs(un - u)))
            hist.append(d)
            u = un
            if d <= tol:
                break
        return u, hist

    def bc_check(self, phi: np.ndarray):
        """D^{nu} u(1) and int u dA for u = T_h phi (both from the same quadrature)."""
        p = self.p
        C = float(self.Cvec @ phi)
        lhs = C * p.Ga                                   # D^{nu}u(1) = -int phi + c1 Ga
        Mphi = float(self.M @ phi)
        uxi = (p.xi_a ** p.nu * Mphi - self.Wxi @ phi) / p.Ga + p.xi_a ** p.nu * C
        rhs = float(np.sum(p.w_a * uxi))
        return lhs, rhs


# ----------------------------------------------------------------------------
# the nonlinearity of the examples
# ----------------------------------------------------------------------------
def f_example(t, z):
    """f(t,z) = 11/10 + sin(2 pi z)/10, so f_- = 1, f_+ = 6/5, Lipschitz constant pi/5."""
    return 1.1 + 0.1 * np.sin(2.0 * np.pi * np.asarray(z, float))


F_MINUS, F_PLUS, F_LIP = 1.0, 1.2, np.pi / 5.0


def solve_problem(prob: Problem, n=2400, q=2.0, extra=(), f=f_example, tol=1e-15):
    pts = tuple(extra) + tuple(prob.xi_a) + (0.25, 0.5, 0.75, 1.0)
    s = make_mesh(n, q, pts)
    D = Discretisation(prob, s)
    u, hist = D.solve(f, tol=tol)
    return D, u, hist


def locate_max(D: Discretisation, u: np.ndarray, f=f_example):
    """Maximum of the natural interpolant of the discrete solution."""
    phi = f(D.s, u)
    i = int(np.argmax(u))
    lo = D.s[max(i - 2, 0)]
    hi = D.s[min(i + 2, D.s.size - 1)]
    for _ in range(4):
        tt = np.linspace(lo, hi, 201)
        vv = D.evaluate(phi, tt)
        k = int(np.argmax(vv))
        lo, hi = tt[max(k - 1, 0)], tt[min(k + 1, tt.size - 1)]
    tt = np.linspace(lo, hi, 201)
    vv = D.evaluate(phi, tt)
    k = int(np.argmax(vv))
    return float(tt[k]), float(vv[k])
