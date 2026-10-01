"""Search related to Question 2 of Section 8 (a relaxed version only).

For boundary points x of the wedge P_sigma(beta_1, b) (beta_1(x) = x(t2) = b) the
three numbers X1 = x(t1), Y1 = Tx(t1), Y2 = Tx(t2) satisfy, by (47) and (43),
    X1 in [q b / Theta, q b Theta],  Y1 in [f_- psi_1, f_+ psi_1],  Y2 in [f_- psi_2, f_+ psi_2],
with q = psi_1/psi_2 and Theta = f_+/f_-.  Here the three numbers are allowed to
vary independently over this box B, which is a larger set than the true set S of
Section 8.  After scaling we may take psi_2 = 1, f_- = 1, f_+ = Theta, psi_1 = q.

A choice of constants (q, Theta, b, x0(t1) = p1, x0(t2) = p2, a) is called
admissible if
  (1) (A1) holds:  a <= p1,  p2 < b,  and  q p2/Theta <= p1 <= q p2 Theta
      (the last two inequalities are (47) for x0 in the cone);
  (2) outward points exist in the box:  Theta > b  (i.e. f_+ psi_2 > b);
  (3) the omitted ray inequality (62) holds at every box point on the branch
      X1 >= a of (A2) and at every box point on the branch Y1 < a of (A3).
Failure of (62) is decided with Theorem 6, i.e. by testing (F1) or (F2), which
involves only non-strict inequalities and is safe in floating point.
The box is sampled on a grid that contains its corners and the values a, b.

Output: results/box_scan.json
"""
import itertools

import numpy as np

from _setup import save_json


def fails(u, v, d1, d2):
    f1 = (u * d1 >= 0) & (np.abs(u) >= v) & (np.abs(d1) >= np.abs(d2))
    f2 = (d2 >= 0) & (v >= np.abs(u)) & (np.abs(d2) >= np.abs(d1))
    return f1 | f2


def main(nbox=17):
    Qs = (0.3, 0.5, 0.8, 1.0, 1.25, 1.5, 2.0, 3.0, 5.0)
    Thetas = (1.05, 1.2, 1.5, 2.0, 3.0)
    tested = admissible = 0
    examples = []
    for q, Th in itertools.product(Qs, Thetas):
        sig = 1.0 / Th
        for b in np.linspace(0.1 * Th, 0.95 * Th, 10):          # Theta > b: outward points exist
            I = np.linspace(sig * q * b, q * b / sig, nbox)       # X1 = x(t1)
            J1 = np.linspace(q, Th * q, nbox)                     # Y1 = Tx(t1)
            J2 = np.linspace(1.0, Th, nbox)                       # Y2 = Tx(t2)
            for p2 in np.linspace(0.05 * b, 0.95 * b, 6):
                v = b - p2
                for p1 in np.linspace(sig * q * p2, q * p2 / sig, 6):
                    for a in np.linspace(0.0, p1, 6):
                        tested += 1
                        X1 = np.union1d(I, [a] if I[0] <= a <= I[-1] else [])
                        Y1 = np.union1d(J1, [a] if J1[0] <= a <= J1[-1] else [])
                        Y2 = np.union1d(J2, [b] if J2[0] <= b <= J2[-1] else [])
                        X, Y, Z = np.meshgrid(X1, Y1, Y2, indexing="ij")
                        u, d1, d2 = X - p1, Y - X, Z - b
                        bad = fails(u, v, d1, d2)
                        on_branch = (X >= a) | (Y < a)
                        if np.any(bad & on_branch):
                            continue
                        admissible += 1
                        if len(examples) < 20:
                            examples.append(dict(q=q, Theta=Th, b=float(b), p1=float(p1),
                                                 p2=float(p2), a=float(a)))
    out = dict(parameter_sets=tested, admissible=admissible, box_grid=nbox, examples=examples)
    save_json("box_scan.json", out)
    print(out)


if __name__ == "__main__":
    main()
