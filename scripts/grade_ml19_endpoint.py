"""Grade the ml19 60 s endpoint run against the full Krueger base gates.

ml18 (12 s) established that the deposition-crosslinking correction
(`3a931b1`) closes the early transient: closure/etch fell from 5.07x / 3.52x /
1.69x Krueger's 0.0310 to 1.83x / 1.59x / 1.18x, and the aperture at t = 12 s
went 38.25 -> 64.89 nm.  ml19 runs the same configuration to 60 s so the
endpoint metrics --- neck CD and its depth, etch depth, mask survival --- can
be graded against the published targets for the first time with every
mechanism either source-verbatim, measurement-gated, or declared-open.

Targets (Krueger 2024, JVST A 42, 043008 + thesis ch. 6):

  neck CD          38.8 nm @ 271 nm below mask top   (his MCFPM)
                   39.0 nm @ 200 nm                  (his SEM)
  etch depth       825 nm +/- 5 %
  mask remaining   850 nm (= initial thickness; armored a-C)
  closure/etch     0.0310 run average

    python scripts/grade_ml19_endpoint.py
"""

from __future__ import annotations

import bisect
import json
import pathlib

CURATED = pathlib.Path("results/curated")
FEATURE = CURATED / "mixed_layer_feature_v1"

INITIAL_OPENING_NM = 90.0
KRUEGER_RATIO = (0.5 * (INITIAL_OPENING_NM - 38.8) / 60.0) / (825.0 / 60.0)

RUN = "ml19-depxl-60s"
BASELINE = "ml16a-verbatim-lift"

WINDOWS = ((1, 4), (4, 8), (8, 12), (12, 20), (20, 40), (40, 60))

# (label, target, tolerance, metric key)
ENDPOINT_GATES = (
    ("etch depth (nm)", 825.0, 0.05, "etch_depth_nm"),
    ("mask remaining (nm)", 850.0, 0.02, "remaining_mask_thickness_nm"),
)


def load(tag):
    path = FEATURE / tag / "audit.json"
    if not path.exists():
        return None
    payload = json.loads(path.read_text())
    history = payload["history"]
    times, openings, depths = [], [], []
    for record in history:
        metrics = record["metrics"]
        times.append(record["physical_time_s"])
        openings.append(metrics["mask_opening_nm"])
        depths.append(metrics["etch_depth_nm"])
    # A run that stops on a topology refusal never writes ``final_metrics``;
    # the last history record carries the same fields.
    final = payload.get("final_metrics") or history[-1]["metrics"]
    return times, openings, depths, final


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


def verdict(value, target, tol):
    return "PASS" if abs(value - target) <= tol * target else "MISS"


def main():
    run = load(RUN)
    if run is None:
        raise SystemExit(f"missing {FEATURE / RUN / 'audit.json'}")
    times, openings, depths, final = run
    base = load(BASELINE)

    print(f"=== ml19 endpoint (t = {times[-1]:.1f} s) ===\n")

    print("ENDPOINT GATES")
    for label, target, tol, key in ENDPOINT_GATES:
        if key not in final:
            print(f"  {label:<26} MISSING")
            continue
        value = float(final[key])
        print(f"  {label:<26} {value:8.1f}  target {target:.0f} "
              f"+/-{100 * tol:.0f}%   {verdict(value, target, tol)}")

    # The MASK constriction is the object Krueger reports as ``w_m``.
    # ``neck_cd_nm`` is a global minimum over the whole feature, which at
    # partial etch depth is the advancing etch-front taper, not a mask neck
    # (same caveat recorded for ml18) --- so grade the mask throat explicitly.
    opening = final.get("mask_opening_nm")
    throat_z = final.get("mask_opening_throat_z_um")
    mask_top = final.get("mask_top_z_um")
    top = final.get("top_cd_nm")
    if opening is not None:
        print(f"\n  MASK constriction CD     {float(opening):8.2f}  "
              f"targets 45 (his w_m) / 38.8 (his sim neck)")
    if throat_z is not None and mask_top is not None:
        below = 1000.0 * (float(mask_top) - float(throat_z))
        print(f"  its depth below mask top {below:8.1f}  "
              f"targets 271 (sim) / 200 (SEM)")
    if top is not None:
        print(f"  top CD                   {float(top):8.2f}")
    neck = final.get("neck_cd_nm")
    neck_z = final.get("neck_depth_from_mask_top_nm")
    if neck is not None and neck_z is not None:
        print(f"  (global min {float(neck):.2f} @ {float(neck_z):.0f} nm "
              f"= etch-front taper, not a mask neck)")

    # Depth honesty: a run that stops early can sit inside the endpoint band
    # while running too fast.  Report the rate against Krueger's average.
    if times[-1] < 59.0:
        late = [i for i, t in enumerate(times) if t >= times[-1] - 6.0]
        rate = ((depths[-1] - depths[late[0]])
                / (times[-1] - times[late[0]]))
        krueger_rate = 825.0 / 60.0
        print(f"\n  RUN STOPPED EARLY at t={times[-1]:.1f}s --- endpoint band "
              f"is not a pass.")
        print(f"  late etch rate           {rate:8.2f} nm/s  vs Krueger "
              f"average {krueger_rate:.2f}  ({rate / krueger_rate:.2f}x)")
        print(f"  linear extrapolation to 60 s: "
              f"{depths[-1] + rate * (60.0 - times[-1]):.0f} nm "
              f"({100 * (depths[-1] + rate * (60.0 - times[-1]) - 825) / 825:+.0f}% "
              f"vs 825)")

    print(f"\nCLOSURE/ETCH BY WINDOW (Krueger reference {KRUEGER_RATIO:.4f})")
    base_rows = windows(*base[:3]) if base else []
    for lo, hi, per_side, etch, ratio in windows(times, openings, depths):
        b = [r for r in base_rows if r[0] == lo]
        label = f"{lo}-{hi:g}s"
        line = (f"  {label:>10} {ratio:8.4f}  {ratio / KRUEGER_RATIO:5.2f}x")
        if b:
            line += (f"   baseline {b[0][4]:7.4f}  "
                     f"{b[0][4] / KRUEGER_RATIO:5.2f}x")
        print(line)

    lost = INITIAL_OPENING_NM - sample(times, openings, times[-1])
    budget = INITIAL_OPENING_NM - 38.8
    print(f"\nclosure budget spent: {lost:.1f} nm = "
          f"{100 * lost / budget:.0f}% of Krueger's full-run budget")

    if base:
        print("\nmatched-time trajectory (ml19 vs ml16a baseline):")
        b_times, b_openings, b_depths = base[:3]
        for t in (4, 8, 12, 20, 30, 40, 50, 60):
            if times[-1] < t:
                continue
            o, d = sample(times, openings, t), sample(times, depths, t)
            bo, bd = sample(b_times, b_openings, t), sample(b_times, b_depths, t)
            print(f"  t={t:5.1f}s  opening {o:7.2f} vs {bo:7.2f}   "
                  f"depth {d:7.1f} vs {bd:7.1f} "
                  f"({100 * (d - bd) / bd:+5.1f}%)")


if __name__ == "__main__":
    main()
