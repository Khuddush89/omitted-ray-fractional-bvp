"""Run every computation of the paper and rebuild all tables and figures.

    python scripts/make_all.py            # everything (about 1-2 minutes)
    python scripts/make_all.py --quick    # skip the slow scans

Outputs go to results/ (JSON files and LaTeX table rows) and figures/ (PDF).
"""
import importlib
import sys
import time

import _setup  # noqa: F401  (sets up paths)

STEPS = [
    ("compute_constants", "closed-form constants, (E1)-(E4), Table 3"),
    ("solve_examples", "numerical solutions, Tables 7-9"),
    ("check_lemmas", "Tables 4-6 (Lemmas 5, 8 and Proposition 1)"),
    ("check_kernel", "Theorem 2 against exact solutions (Appendix A)"),
    ("check_characterisation", "Theorem 6 in exact rational arithmetic (Appendix A)"),
    ("check_forced_inward", "random test of Proposition 2 (Appendix A)"),
    ("box_scan", "relaxed search for Question 2 (Section 8)"),
    ("make_figures", "Figures 1-8"),
]
SLOW = {"box_scan", "check_characterisation"}


def main():
    quick = "--quick" in sys.argv
    t_all = time.time()
    for name, what in STEPS:
        if quick and name in SLOW:
            print(f"[skip] {name}: {what}")
            continue
        print(f"\n=== {name}: {what}")
        t0 = time.time()
        importlib.import_module(name).main()
        print(f"--- {name} done in {time.time() - t0:.1f} s")
    print(f"\nAll done in {time.time() - t_all:.1f} s")


if __name__ == "__main__":
    main()
