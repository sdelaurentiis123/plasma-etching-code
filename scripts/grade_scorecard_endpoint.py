"""Grade the 8-condition scorecard at its published 60 s ENDPOINTS.

The 2026-08-06 pass could only grade a matched early time (t = 1.948 s) because
the box it drew exposed no CUDA device, so the endpoint-only criteria -- clog,
necking-absent, and absolute depth -- were explicitly refused.  This grader is
the endpoint counterpart: it reads the completed audits and grades every
criterion the scorecard actually publishes.

Criteria (scorecard-1, results/curated/mixed_layer_scorecard_1/):

  base      etch depth 825 nm +/- 5 %; mask opening 45 +/- 5 nm
  O2 0.5    clog preserved -- the aperture must SEAL (0.000 nm)
  O2 1.5    depth ranks maximum over the oxygen sweep
  O2 2.5    necking absent -- the aperture stays open
  4 kW      r(4/6) = depth(4 kW) / depth(base) in [0.84, 0.94]
  8 kW      r(8/6) = depth(8 kW) / depth(base) in [0.97, 1.06]

Endpoint states are only meaningful when the run actually reached 60 s, so
every row reports the simulated time it was graded at and a row whose run
stopped short is marked INCOMPLETE rather than silently graded early -- the
same refusal discipline as the matched-time grader.

    python scripts/grade_scorecard_endpoint.py [RESULTS_DIR]
"""

from __future__ import annotations

import json
import pathlib
import sys

TAGS = ("sc-base", "sc-o2-0.5", "sc-o2-1.5", "sc-o2-2.5", "sc-p4", "sc-p8")
LABEL = {
    "sc-base": "base (6 kW, O2 1.0)",
    "sc-o2-0.5": "O2 0.5",
    "sc-o2-1.5": "O2 1.5",
    "sc-o2-2.5": "O2 2.5",
    "sc-p4": "4 kW",
    "sc-p8": "8 kW",
}
DURATION_S = 60.0
COMPLETE_FRACTION = 0.98  # graded as an endpoint only within 2 % of 60 s


def load(directory, tag):
    path = pathlib.Path(directory) / tag / "audit.json"
    if not path.exists():
        return None
    audit = json.loads(path.read_text())
    history = audit.get("history") or []
    if not history:
        return None
    last = history[-1]
    return {
        "status": audit.get("status", ""),
        "t": float(last.get("physical_time_s", 0.0)),
        "metrics": last.get("metrics", {}),
    }


def verdict(ok):
    return "PASS" if ok else "MISS"


def main():
    directory = sys.argv[1] if len(sys.argv) > 1 else "results/curated/scorecard_endpoint"
    runs = {tag: load(directory, tag) for tag in TAGS}
    missing = [t for t, r in runs.items() if r is None]
    if missing:
        print(f"no audit for: {', '.join(missing)}")
    runs = {t: r for t, r in runs.items() if r is not None}
    if not runs:
        sys.exit("nothing to grade")

    print(f"{'condition':<22}{'t (s)':>8}{'depth (nm)':>13}{'opening (nm)':>14}"
          f"{'neck (nm)':>11}{'mask (nm)':>11}  state")
    for tag in TAGS:
        run = runs.get(tag)
        if run is None:
            continue
        m = run["metrics"]
        state = ("endpoint" if run["t"] >= COMPLETE_FRACTION * DURATION_S
                 else "INCOMPLETE")
        print(f"{LABEL[tag]:<22}{run['t']:>8.2f}{m.get('etch_depth_nm', 0):>13.2f}"
              f"{m.get('mask_opening_nm', 0):>14.3f}{m.get('neck_cd_nm', 0):>11.2f}"
              f"{m.get('remaining_mask_thickness_nm', 0):>11.2f}  {state}")

    def endpoint(tag):
        run = runs.get(tag)
        if run is None or run["t"] < COMPLETE_FRACTION * DURATION_S:
            return None
        return run

    print("\ncriterion                         target            measured      verdict")
    base = endpoint("sc-base")
    if base:
        d = base["metrics"]["etch_depth_nm"]
        ok = abs(d - 825.0) <= 0.05 * 825.0
        print(f"{'base etch depth':<34}{'825 +/- 5%':<18}{d:>10.1f}      {verdict(ok)}")
        o = base["metrics"]["mask_opening_nm"]
        ok = abs(o - 45.0) <= 5.0
        print(f"{'base mask opening':<34}{'45 +/- 5':<18}{o:>10.2f}      {verdict(ok)}")
        mk = base["metrics"]["remaining_mask_thickness_nm"]
        ok = abs(mk - 850.0) <= 0.02 * 850.0
        print(f"{'base mask remaining':<34}{'850 +/- 2%':<18}{mk:>10.2f}      {verdict(ok)}")
    else:
        print("base did not reach 60 s -- absolute gates refused")

    clog = endpoint("sc-o2-0.5")
    if clog:
        o = clog["metrics"]["mask_opening_nm"]
        ok = o <= 1e-9
        print(f"{'clog at O2 0.5':<34}{'sealed (0.000)':<18}{o:>10.3f}      {verdict(ok)}")

    neck = endpoint("sc-o2-2.5")
    if neck:
        o = neck["metrics"]["mask_opening_nm"]
        ok = o > 0.0
        print(f"{'necking absent at O2 2.5':<34}{'open':<18}{o:>10.3f}      {verdict(ok)}")

    o2 = {t: endpoint(t) for t in ("sc-o2-0.5", "sc-base", "sc-o2-1.5", "sc-o2-2.5")}
    if all(o2.values()):
        depths = {t: o2[t]["metrics"]["etch_depth_nm"] for t in o2}
        top = max(depths, key=depths.get)
        ok = top == "sc-o2-1.5"
        order = " | ".join(f"{LABEL[t].split()[-1]}:{depths[t]:.0f}" for t in o2)
        print(f"{'O2 depth rank max at 1.5':<34}{'1.5 highest':<18}"
              f"{LABEL[top].split()[-1]:>10}      {verdict(ok)}")
        print(f"    depths -> {order}")

    if base:
        for tag, band, name in (("sc-p4", (0.84, 0.94), "r(4/6)"),
                                ("sc-p8", (0.97, 1.06), "r(8/6)")):
            run = endpoint(tag)
            if not run:
                print(f"{name + ' (endpoint)':<34}{'run short':<18}"
                      f"{'--':>10}      REFUSED")
                continue
            r = run["metrics"]["etch_depth_nm"] / base["metrics"]["etch_depth_nm"]
            ok = band[0] <= r <= band[1]
            print(f"{name:<34}{f'[{band[0]}, {band[1]}]':<18}{r:>10.3f}      "
                  f"{verdict(ok)}")


if __name__ == "__main__":
    main()
