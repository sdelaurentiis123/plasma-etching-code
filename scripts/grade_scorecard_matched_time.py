"""Grade the 8-condition scorecard at a matched simulated time.

The scorecard's 60 s endpoints are the published comparison, but a CPU-only
box runs ~90-160 s/step at dx = 10 nm, which puts a 60 s endpoint ~4 h away per
condition.  This harness grades what a matched EARLY time can legitimately
certify and refuses to imply more:

  gradeable at matched time
    - the oxygen ORDERING (film thickness is monotone in O, so the aperture
      ordering base < 1.5 < 2.5 and the 0.5 outlier must hold from the start)
    - the power DEPTH RATIOS r(4/6) and r(8/6), which are ratios of depths at
      the same simulated instant and so are dimensionally the published
      quantity even though the absolute depths are early
  NOT gradeable at matched time
    - clog/necking topology verdicts, which are endpoint states
    - absolute depth against 825 nm

Usage: python scripts/grade_scorecard_matched_time.py STEPLOG_DIR [T_MATCH]
"""
from __future__ import annotations

import pathlib
import re
import sys

STEP = re.compile(r"^step=(\d+)\s+t=([0-9.]+)s.*depth=([0-9.]+)nm\s+"
                  r"mask_opening=([0-9.]+)nm")


def trajectory(path):
    rows = []
    for line in pathlib.Path(path).read_text(errors="ignore").splitlines():
        m = STEP.match(line.strip())
        if m:
            rows.append((float(m.group(2)), float(m.group(3)), float(m.group(4))))
    return rows


def at_time(rows, t):
    """Linear interpolation of (depth, opening) at simulated time t."""
    if not rows or rows[-1][0] < t:
        return None
    prev = rows[0]
    for row in rows:
        if row[0] >= t:
            if row[0] == prev[0]:
                return row[1], row[2]
            f = (t - prev[0]) / (row[0] - prev[0])
            return (prev[1] + f * (row[1] - prev[1]),
                    prev[2] + f * (row[2] - prev[2]))
        prev = row
    return None


def main():
    d = pathlib.Path(sys.argv[1])
    t_match = float(sys.argv[2]) if len(sys.argv) > 2 else None
    logs = {p.stem: trajectory(p) for p in sorted(d.glob("sc-*.log"))}
    logs = {k: v for k, v in logs.items() if v}
    if not logs:
        print("no step records found")
        return
    reach = {k: v[-1][0] for k, v in logs.items()}
    if t_match is None:
        t_match = min(reach.values())
    print(f"matched simulated time: t = {t_match:.3f} s")
    print(f"{'condition':<12} {'reached':>8} {'depth':>9} {'opening':>9}")
    vals = {}
    for name, rows in sorted(logs.items()):
        got = at_time(rows, t_match)
        if got is None:
            print(f"{name:<12} {reach[name]:8.2f}   (short of matched time)")
            continue
        vals[name] = got
        print(f"{name:<12} {reach[name]:8.2f} {got[0]:9.3f} {got[1]:9.3f}")

    print("\n-- oxygen ordering (aperture must widen with O) --")
    order = [("sc-o2-0.5", 0.5), ("sc-base", 1.0), ("sc-o2-1.5", 1.5),
             ("sc-o2-2.5", 2.5)]
    seq = [(lbl, vals[k][1]) for k, lbl in order if k in vals]
    print("  " + "  ".join(f"O2 {l}: {v:.2f}" for l, v in seq))
    if len(seq) >= 2:
        mono = all(b >= a for (_, a), (_, b) in zip(seq, seq[1:]))
        print(f"  monotone in O2: {mono}")

    print("\n-- power depth ratios at matched time --")
    if "sc-base" in vals:
        base = vals["sc-base"][0]
        for k, lbl, lo, hi in (("sc-p4", "r(4/6)", 0.84, 0.94),
                               ("sc-p8", "r(8/6)", 0.97, 1.06)):
            if k in vals:
                r = vals[k][0] / base
                verdict = "in band" if lo <= r <= hi else "MISS"
                print(f"  {lbl} = {r:.3f}   band [{lo}, {hi}]   {verdict}")


if __name__ == "__main__":
    main()
