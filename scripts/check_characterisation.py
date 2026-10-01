"""Exact-arithmetic check of Theorem 6 (Appendix A).

For random rational quadruples (u, v, d1, d2) with v > 0 we compare
    failure of  max{|u+d1|, |v+d2|} < max{|u|, v} + max{|d1|, |d2|}
with the disjunction (F1) or (F2).  Since the left side never exceeds the right
side, failure means exact equality, so floating point cannot be used here.
The values are drawn from a coarse set of rationals so that the boundary cases
(|u| = v, |d1| = |d2|, u d1 = 0, d2 = 0) occur often.

Output: results/characterisation_check.json
"""
import random
from fractions import Fraction as Fr

from _setup import save_json

N = 300_000
SEED = 20260929


def fails(u, v, d1, d2):
    L = max(abs(u + d1), abs(v + d2))
    R = max(abs(u), v) + max(abs(d1), abs(d2))
    assert L <= R
    return L == R


def F1(u, v, d1, d2):
    return u * d1 >= 0 and abs(u) >= v and abs(d1) >= abs(d2)


def F2(u, v, d1, d2):
    return d2 >= 0 and v >= abs(u) and abs(d2) >= abs(d1)


def main():
    rng = random.Random(SEED)
    vals = [Fr(k, 4) for k in range(-12, 13)]
    pos = [Fr(k, 4) for k in range(1, 13)]
    mismatches, n_fail, n_f1, n_f2 = 0, 0, 0, 0
    for _ in range(N):
        u, d1, d2 = rng.choice(vals), rng.choice(vals), rng.choice(vals)
        v = rng.choice(pos)
        a = fails(u, v, d1, d2)
        f1, f2 = F1(u, v, d1, d2), F2(u, v, d1, d2)
        n_fail += a
        n_f1 += f1
        n_f2 += f2
        if a != (f1 or f2):
            mismatches += 1
    # the two quadruples quoted in the paper
    q1 = (Fr(1, 5), Fr(1), Fr(-2), Fr(1, 2))   # Remark 13: holds, no item of Cor. 5 applies
    q2 = (Fr(2), Fr(1), Fr(2), Fr(1))          # Section 8: (F1) holds, (F2) fails, L = R = 4
    out = dict(samples=N, seed=SEED, mismatches=mismatches, failures=n_fail,
               F1_true=n_f1, F2_true=n_f2,
               remark_quadruple=dict(L=str(max(abs(q1[0] + q1[2]), abs(q1[1] + q1[3]))),
                                     R=str(max(abs(q1[0]), q1[1]) + max(abs(q1[2]), abs(q1[3]))),
                                     F1=F1(*q1), F2=F2(*q1)),
               conclusion_quadruple=dict(L=str(max(abs(q2[0] + q2[2]), abs(q2[1] + q2[3]))),
                                         R=str(max(abs(q2[0]), q2[1]) + max(abs(q2[2]), abs(q2[3]))),
                                         F1=F1(*q2), F2=F2(*q2)))
    save_json("characterisation_check.json", out)
    print(out)


if __name__ == "__main__":
    main()
