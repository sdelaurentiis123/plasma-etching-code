"""Grade the Gray-law run (ml23) against the pre-Gray baselines.

Two questions this answers, both from step-log trajectories so no endpoint is
required:

1. Matched-simulated-time comparison -- depth and aperture at the same physical
   time, against ml19 (deposition-crosslinking, pre-Gray) and ml21 (final
   mechanism, pre-Gray).  Rates are compared where the runs overlap, which is
   the only honest comparison when one run is far slower than another.

2. The ARDE sign.  The trench deepens as it etches, so the instantaneous floor
   rate as a function of accumulated depth *is* a rate-vs-aspect-ratio curve
   for a single run.  Pre-Gray configurations measured a FLAT or RISING rate
   (RESULTS_FLOOR_DELIVERY / the anti-ARDE finding); the regime claim is that
   Gray's laws make it fall.  Reported as the fitted slope of rate against
   depth over the run's own span, with the first two steps dropped (startup
   transient, dt still adapting).
"""

import re
import sys

_STEP = re.compile(
    r"step=(\d+)\s+t=([\d.]+)s\s+dt=([\d.eE+-]+)s\s+depth=([\d.]+)nm\s+"
    r"mask_opening=([\d.]+)nm")


def read_steps(path):
    out = []
    with open(path, errors="ignore") as handle:
        for line in handle:
            m = _STEP.search(line)
            if m:
                out.append({
                    "step": int(m.group(1)), "t": float(m.group(2)),
                    "dt": float(m.group(3)), "depth": float(m.group(4)),
                    "open": float(m.group(5))})
    return out


def rate_series(steps, skip=2):
    """Instantaneous floor rate and aperture closure rate per step."""
    series = []
    for prev, cur in zip(steps[skip:], steps[skip + 1:]):
        dt = cur["t"] - prev["t"]
        if dt <= 0:
            continue
        series.append({
            "t": cur["t"],
            "depth": cur["depth"],
            "rate": (cur["depth"] - prev["depth"]) / dt,
            "closure": (prev["open"] - cur["open"]) / dt / 2.0,  # per side
        })
    return series


def slope(xs, ys):
    n = len(xs)
    if n < 2:
        return float("nan")
    mx = sum(xs) / n
    my = sum(ys) / n
    den = sum((x - mx) ** 2 for x in xs)
    if den == 0:
        return float("nan")
    return sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / den


def at_time(steps, t):
    """Linear interpolation of the trajectory at a physical time."""
    if not steps or t < steps[0]["t"] or t > steps[-1]["t"]:
        return None
    for a, b in zip(steps, steps[1:]):
        if a["t"] <= t <= b["t"]:
            if b["t"] == a["t"]:
                return a
            f = (t - a["t"]) / (b["t"] - a["t"])
            return {"t": t,
                    "depth": a["depth"] + f * (b["depth"] - a["depth"]),
                    "open": a["open"] + f * (b["open"] - a["open"])}
    return steps[-1]


def main(paths):
    runs = {}
    for label, path in paths:
        steps = read_steps(path)
        if steps:
            runs[label] = steps
        else:
            print(f"  (no step records in {path})")

    if not runs:
        print("no runs parsed")
        return 1

    print("== span ==")
    for label, steps in runs.items():
        print(f"{label:18s} {len(steps):4d} steps, t=0..{steps[-1]['t']:.3f}s, "
              f"depth 0..{steps[-1]['depth']:.1f}nm, "
              f"aperture {steps[0]['open']:.2f}->{steps[-1]['open']:.2f}nm")

    # Matched-time table over the shared span.
    horizon = min(s[-1]["t"] for s in runs.values())
    marks = [t for t in (0.5, 1.0, 2.0, 4.0, 8.0, 12.0) if t <= horizon]
    if marks:
        print(f"\n== matched simulated time (shared span 0..{horizon:.2f}s) ==")
        head = "  t(s) | " + " | ".join(f"{l:>22s}" for l in runs)
        print(head)
        print("  " + "-" * (len(head) - 2))
        for t in marks:
            cells = []
            for steps in runs.values():
                p = at_time(steps, t)
                cells.append(f"{p['depth']:8.2f}nm {p['open']:7.2f}nm"
                             if p else "            --        ")
            print(f"  {t:5.2f} | " + " | ".join(f"{c:>22s}" for c in cells))

    print("\n== floor rate vs depth (the ARDE sign) ==")
    for label, steps in runs.items():
        ser = rate_series(steps)
        if len(ser) < 3:
            print(f"{label:18s} too few steps")
            continue
        xs = [p["depth"] for p in ser]
        ys = [p["rate"] for p in ser]
        m = slope(xs, ys)
        sign = ("DECELERATES (ARDE)" if m < -1e-4 else
                "ACCELERATES (anti-ARDE)" if m > 1e-4 else "FLAT (no ARDE)")
        print(f"{label:18s} rate {ys[0]:6.2f} -> {ys[-1]:6.2f} nm/s over "
              f"depth {xs[0]:6.1f} -> {xs[-1]:6.1f} nm | "
              f"slope {m:+.5f} (nm/s)/nm | {sign}")

    print("\n== aperture closure per side (nm/s) ==")
    for label, steps in runs.items():
        ser = rate_series(steps)
        if not ser:
            continue
        ys = [p["closure"] for p in ser]
        rates = [p["rate"] for p in ser]
        ratios = [c / r for c, r in zip(ys, rates) if r > 1e-9]
        mean_ratio = sum(ratios) / len(ratios) if ratios else float("nan")
        print(f"{label:18s} closure {ys[0]:6.3f} -> {ys[-1]:6.3f} | "
              f"mean closure/etch {mean_ratio:.4f} (Krueger 0.0310)")
    return 0


if __name__ == "__main__":
    args = sys.argv[1:]
    if not args:
        print("usage: grade_ml23_gray.py label=path [label=path ...]")
        raise SystemExit(2)
    raise SystemExit(main([tuple(a.split("=", 1)) for a in args]))
