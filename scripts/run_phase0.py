#!/usr/bin/env python3
"""Phase-0 scorecard: run each model configuration through the fidelity benchmark and the
convergence suite, and print a table.

Usage:
    python scripts/run_phase0.py            # full: scorecards + convergence
    python scripts/run_phase0.py --quick    # scorecards only (skip convergence)
    python scripts/run_phase0.py --width8    # just the width-8 single number, fastest
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "harness"))

import petch
from benchmark import scorecard, print_scorecard, run_case, load_reference
import convergence as conv


def configs():
    """Model configurations to score. Extended as contributors are retired."""
    rows = []
    base = dict(petch.PAR); base['rate_scale'] = 0.29
    rows.append(("baseline langmuir (rate_scale=0.29)", base, petch.Flags(chemistry="langmuir")))
    # Belen row (contributor #1) — added once chemistry='belen' is implemented (Step 3).
    try:
        from petch.belen import surface_rate_belen  # noqa: F401
        belen_par = dict(petch.PAR)  # rate_scale stays 1.0
        rows.append(("belen, rate_scale=1.0 (contributor #1)", belen_par,
                     petch.Flags(chemistry="belen")))
    except Exception:
        pass
    return rows


def main():
    args = set(sys.argv[1:])
    ref = load_reference()
    print("=" * 74)
    print("PHASE-0 SCORECARD  —  2D feature-scale etch vs cached ViennaPS reference")
    print(f"  ViennaPS reference: width-8 depth {ref['depth_vps']:.3f} um, "
          f"ARDE {['%.2f'%d for d in ref['vps_depth']]} for widths {ref['widths']}")
    print("=" * 74)

    if "--width8" in args:
        for label, par, flags in configs():
            r = run_case(par, flags, 8.0)
            print(f"  [{label}]  width-8 depth = {r['_depth']:.3f} um  "
                  f"(ViennaPS {ref['depth_vps']:.3f})")
        return

    print("\n--- ACCURACY + SPEED ---")
    for label, par, flags in configs():
        print_scorecard(scorecard(par, flags, label=label))
        print()

    if "--quick" not in args:
        print("--- CONVERGENCE (compact geometry; behavior, not benchmark depth) ---")
        base = dict(petch.PAR); base['rate_scale'] = 0.29
        flags = petch.Flags(chemistry="langmuir")
        print("  grid refinement (dx -> 0):")
        conv.grid_refinement(base, flags)
        print("  ray refinement (N -> inf):")
        conv.ray_refinement(base, flags)

    print("\ndone.")


if __name__ == "__main__":
    main()
