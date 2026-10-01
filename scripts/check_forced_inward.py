"""Random test of Proposition 2: (E1) and (E2) imply f_+ psi_2 < b (Appendix A).

Parameters are drawn at random; among the draws that satisfy (E1) (with a = 0,
the weakest choice) and (E2), we count how often f_+ psi_2 >= b.  The proposition
says this never happens.

Output: results/forced_inward_check.json
"""
import numpy as np

from _setup import save_json

N = 2_000_000
SEED = 20260929


def main():
    rng = np.random.default_rng(SEED)
    fm = rng.uniform(0.05, 5.0, N)
    fp = fm * rng.uniform(1.0, 4.0, N)
    psi1 = rng.uniform(0.01, 3.0, N)
    psi2 = rng.uniform(0.01, 3.0, N)
    A0 = rng.uniform(0.0, 6.0, N)
    b = rng.uniform(0.0, 10.0, N)
    E1 = A0 * psi2 < b
    E2 = np.maximum(np.abs(fp - A0), np.abs(fm - A0)) * np.maximum(psi1, psi2) < b - A0 * psi2
    ok = E1 & E2
    viol = ok & (fp * psi2 >= b)
    out = dict(draws=N, seed=SEED, satisfying_E1_E2=int(ok.sum()), violations=int(viol.sum()))
    save_json("forced_inward_check.json", out)
    print(out)


if __name__ == "__main__":
    main()
