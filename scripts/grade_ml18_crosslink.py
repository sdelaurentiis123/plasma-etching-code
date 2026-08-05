"""Grade the deposition-crosslinking confirmation run (ml18) against baseline.

The mouth defect is an *early transient* (RESULTS_EARLY_TRANSIENT_2026-08-04):
88 % of Krueger's full 60 s closure budget was spent by t = 8 s, with the
closure/etch ratio running 5.1x / 3.5x / 1.7x his 0.0310 over the first three
windows and matching him from t = 20 s onward.  The deposition-crosslinking
correction (RESULTS_LIP_CROSSLINK_2026-08-04, `3a931b1`) targets exactly that
transient, so the run is 12 s and the grade is closure/etch by window --- never
aperture at 60 s.

    python scripts/grade_ml18_crosslink.py
"""

from __future__ import annotations

import bisect
import json
import pathlib

CURATED = pathlib.Path("results/curated")
FEATURE = CURATED / "mixed_layer_feature_v1"

INITIAL_OPENING_NM = 90.0
# Krueger 2024 endpoint: 90 -> 38.8 nm neck, 825 nm depth, 60 s.
KRUEGER_RATIO = (0.5 * (INITIAL_OPENING_NM - 38.8) / 60.0) / (825.0 / 60.0)

WINDOWS = ((1, 4), (4, 8), (8, 12))
BASELINE = "ml16a-verbatim-lift"
RUN = "ml18-depxl-12s"


def load(tag):
    path = FEATURE / tag / "audit.json"
    if not path.exists():
        return None
    history = json.loads(path.read_text())["history"]
    times, openings, depths = [], [], []
    for record in history:
        metrics = record["metrics"]
        times.append(record["physical_time_s"])
        openings.append(metrics["mask_opening_nm"])
        depths.append(metrics["etch_depth_nm"])
    return times, openings, depths


def sample(times, values, t):
    index = min(bisect.bisect_left(times, t), len(values) - 1)
    return values[index]


def windows(times, openings, depths):
    rows = []
    for lo, hi in WINDOWS:
        if times[-1] < lo:
            continue
        hi_eff = min(hi, times[-1])
        if hi_eff - lo < 0.25:
            continue
        o_lo, o_hi = sample(times, openings, lo), sample(times, openings, hi_eff)
        d_lo, d_hi = sample(times, depths, lo), sample(times, depths, hi_eff)
        per_side = 0.5 * (o_lo - o_hi) / (hi_eff - lo)
        etch = (d_hi - d_lo) / (hi_eff - lo)
        rows.append((lo, hi_eff, per_side, etch,
                     per_side / etch if etch else float("nan")))
    return rows


def main():
    base = load(BASELINE)
    run = load(RUN)
    if run is None:
        raise SystemExit(f"missing {FEATURE / RUN / 'audit.json'}")

    print(f"Krueger reference closure/etch = {KRUEGER_RATIO:.4f}\n")
    print(f"{'window':>10} {'run ratio':>10} {'xKrueger':>9} "
          f"{'base ratio':>11} {'xKrueger':>9} {'improve':>8}")
    for row in windows(*run):
        lo, hi, per_side, etch, ratio = row
        b = [r for r in windows(*base) if r[0] == lo]
        label = f"{lo}-{hi:g}s"
        if b:
            b_ratio = b[0][4]
            print(f"{label:>10} {ratio:10.4f} {ratio / KRUEGER_RATIO:8.2f}x "
                  f"{b_ratio:11.4f} {b_ratio / KRUEGER_RATIO:8.2f}x "
                  f"{b_ratio / ratio:7.2f}x")
        else:
            print(f"{label:>10} {ratio:10.4f} {ratio / KRUEGER_RATIO:8.2f}x")

    print("\nmatched-time aperture / depth (run vs baseline):")
    times, openings, depths = run
    b_times, b_openings, b_depths = base
    for t in (0.2, 0.5, 1, 2, 4, 6, 8, 10, 12):
        if times[-1] < t:
            continue
        o, d = sample(times, openings, t), sample(times, depths, t)
        bo, bd = sample(b_times, b_openings, t), sample(b_times, b_depths, t)
        print(f"  t={t:5.1f}s  opening {o:7.2f} vs {bo:7.2f} "
              f"(+{o - bo:5.2f})   depth {d:7.1f} vs {bd:7.1f} "
              f"({100 * (d - bd) / bd:+5.1f}%)")

    lost = INITIAL_OPENING_NM - sample(times, openings, min(12, times[-1]))
    budget = INITIAL_OPENING_NM - 38.8
    print(f"\nclosure budget spent by t={min(12, times[-1]):.1f}s: "
          f"{lost:.1f} nm = {100 * lost / budget:.0f}% of Krueger's full-run "
          f"budget (baseline: 101% by t=12)")


if __name__ == "__main__":
    main()
