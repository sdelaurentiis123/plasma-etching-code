"""Grade the ml21 final-mechanism run against the full Krueger base gates.

ml21 is the first run carrying every mechanism landed in the campaign:
deposition-driven crosslinking with the published per-material bond
multiplicity (`914bb8d`, `924cc04`), the Appendix-B angular classes on the
oxide/mask ion rows (`4c66df1`), the sub-resolution sliver dissolution that
unblocks 60 s completions (`a467fe3`), and the sqrt-2 axisymmetric lift
(`6e97ef3`).

Targets (Krueger 2024, JVST A 42, 043008 + thesis ch. 6):

  mask constriction  45 nm  (his w_m; band [40, 50])
  neck CD            38.8 nm @ 271 nm below mask top  (his MCFPM)
                     39.0 nm @ 200 nm                 (his SEM)
  etch depth         825 nm +/- 5 %
  mask remaining     850 nm (= initial thickness; armored a-C)
  closure/etch       0.0310 run average

Also emits a matched-simulated-time comparison against ml19 (the same
configuration WITHOUT multiplicity + angular classes), which isolates what the
final mechanism changed independently of how far each run got.

    python scripts/grade_ml21_final.py [--run TAG] [--baseline TAG]
"""

from __future__ import annotations

import argparse
import bisect
import json
import pathlib

CURATED = pathlib.Path("results/curated")
FEATURE = CURATED / "mixed_layer_feature_v1"

INITIAL_OPENING_NM = 90.0
# Krueger's run-average closure/etch: half the aperture loss over the etch depth.
KRUEGER_RATIO = (0.5 * (INITIAL_OPENING_NM - 38.8) / 60.0) / (825.0 / 60.0)

WINDOWS = ((1, 4), (4, 8), (8, 12), (12, 20), (20, 40), (40, 60))

ENDPOINT_GATES = (
    ("etch depth (nm)", 825.0, 0.05, "etch_depth_nm"),
    ("mask remaining (nm)", 850.0, 0.02, "remaining_mask_thickness_nm"),
)
MOUTH_BAND = (40.0, 50.0)
NECK_TARGETS = {"mcfpm": (38.8, 271.0), "sem": (39.0, 200.0)}


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
    # A run stopped by a refusal never writes ``final_metrics``; the last
    # history record carries the same fields.
    final = payload.get("final_metrics") or history[-1]["metrics"]
    return times, openings, depths, final, payload


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
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", default="ml21-final")
    ap.add_argument("--baseline", default="ml19-depxl-60s")
    args = ap.parse_args()

    run = load(args.run)
    if run is None:
        raise SystemExit(f"missing {FEATURE / args.run / 'audit.json'}")
    times, openings, depths, final, payload = run

    complete = "final_metrics" in payload
    print(f"=== {args.run} (reached t = {times[-1]:.2f} s; "
          f"{'COMPLETE' if complete else 'PARTIAL'}) ===\n")

    print("ENDPOINT GATES")
    for label, target, tol, key in ENDPOINT_GATES:
        if key not in final:
            print(f"  {label:<26} MISSING")
            continue
        value = float(final[key])
        note = "" if complete else "   (partial run: not an endpoint)"
        print(f"  {label:<26} {value:8.1f}  target {target:.0f} "
              f"+/-{100 * tol:.0f}%   {verdict(value, target, tol)}{note}")

    if "mask_opening_nm" in final:
        mouth = float(final["mask_opening_nm"])
        inband = MOUTH_BAND[0] <= mouth <= MOUTH_BAND[1]
        print(f"  {'mask constriction (nm)':<26} {mouth:8.1f}  band "
              f"[{MOUTH_BAND[0]:.0f}, {MOUTH_BAND[1]:.0f}]     "
              f"{'PASS' if inband else 'MISS'}")
    for key, (cd, depth) in NECK_TARGETS.items():
        have_cd = final.get("neck_cd_nm")
        have_z = final.get("neck_depth_from_mask_top_nm")
        if have_cd is None:
            continue
        print(f"  neck vs {key:<18} {float(have_cd):8.1f} @ "
              f"{float(have_z):.0f} nm   target {cd} @ {depth:.0f} nm")

    print(f"\nCLOSURE / ETCH  (Krueger run-average {KRUEGER_RATIO:.4f})")
    print(f"  {'window (s)':<14}{'per-side (nm/s)':>17}{'etch (nm/s)':>13}"
          f"{'ratio':>9}{'x Krueger':>11}")
    for lo, hi, per_side, etch, ratio in windows(times, openings, depths):
        print(f"  {f'{lo:.0f}-{hi:.1f}':<14}{per_side:17.4f}{etch:13.3f}"
              f"{ratio:9.4f}{ratio / KRUEGER_RATIO:11.2f}")

    base = load(args.baseline)
    if base is not None:
        b_times, b_open, b_depth, _, _ = base
        overlap = min(times[-1], b_times[-1])
        print(f"\nMATCHED-TIME vs {args.baseline} (overlap 0 - {overlap:.2f} s)")
        print(f"  {'t (s)':>8}{'opening ' + args.run[:6]:>16}{'opening base':>14}"
              f"{'d_open':>9}{'depth run':>11}{'depth base':>12}{'d_depth':>9}")
        for t in [x for x in (0.5, 1, 2, 5, 10, 20, 40, 60) if x <= overlap]:
            o_r, o_b = sample(times, openings, t), sample(b_times, b_open, t)
            d_r, d_b = sample(times, depths, t), sample(b_times, b_depth, t)
            print(f"  {t:8.1f}{o_r:16.2f}{o_b:14.2f}{o_r - o_b:9.2f}"
                  f"{d_r:11.2f}{d_b:12.2f}{d_r - d_b:9.2f}")

    ex = payload.get("extrusion_projection_max_deviation_mesh_units")
    if ex is not None:
        print(f"\nextrusion projection max deviation: {float(ex):.3e} mesh units")
    events = payload.get("topology_events") or []
    print(f"topology events: {len(events)}")


if __name__ == "__main__":
    main()
