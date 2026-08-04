"""Early-transient dissection of the mask-opening closure.

Reconstructs the aperture / throat-position history from archived pilot audits
and grades where in *time* the closure budget is spent, against Krueger's
60 s endpoint (neck 38.8 nm at 271 nm depth, from the digitised Fig. 7a).

Also inverts the top-band balance: given the audited, geometry-free
O-removal-to-deposition ratio, what effective sticking would the near-vertical
lip need in order to hold open?  Run:

    python scripts/early_transient_analysis.py
"""

from __future__ import annotations

import bisect
import json
import pathlib

CURATED = pathlib.Path("results/curated")
FEATURE = CURATED / "mixed_layer_feature_v1"
LIP_AUDIT = CURATED / "lip_deposition_audit" / "audit_neck45_dx0.01.json"

# Krueger 2024, thesis section 6.4 (verbatim): "an ideal, straight walled,
# opening with an initial width of 90 nm.  The etch was performed for 60
# seconds."  Digitised Fig. 7a neck minimum: 38.8 nm at 271 nm below mask top.
INITIAL_OPENING_NM = 90.0
KRUEGER_NECK_NM = 38.8
KRUEGER_DEPTH_NM = 825.0
PROCESS_TIME_S = 60.0
KRUEGER_CLOSURE_BUDGET_NM = INITIAL_OPENING_NM - KRUEGER_NECK_NM

RUNS = (
    "ml9a-base-atoms",
    "ml13-base-cascade",
    "ml16a-verbatim-lift",
    "ml16b-ml13c-lift",
)


def load_trajectory(tag):
    """Return (times, openings, depths, throat_depth_below_mask_top_nm)."""
    path = FEATURE / tag / "audit.json"
    if not path.exists():
        return None
    history = json.loads(path.read_text())["history"]
    times, openings, depths, throats = [], [], [], []
    for record in history:
        metrics = record["metrics"]
        times.append(record["physical_time_s"])
        openings.append(metrics["mask_opening_nm"])
        depths.append(metrics["etch_depth_nm"])
        throats.append(
            (metrics["mask_top_z_um"] - metrics["mask_opening_throat_z_um"]) * 1e3)
    return times, openings, depths, throats


def sample(times, values, t):
    index = min(bisect.bisect_left(times, t), len(values) - 1)
    return values[index]


def closure_budget_table(times, openings):
    rows = []
    for t in (2, 4, 6, 8, 10, 12, 20, 30, 40, 50, 60):
        lost = INITIAL_OPENING_NM - sample(times, openings, t)
        rows.append((t, lost, 100.0 * lost / KRUEGER_CLOSURE_BUDGET_NM))
    return rows


def rate_table(times, openings, depths):
    """Per-side closure rate and closure/etch ratio over windows."""
    rows = []
    windows = ((1, 4), (4, 8), (8, 12), (12, 20), (20, 30),
               (30, 40), (40, 50), (50, 60))
    for lo, hi in windows:
        o_lo, o_hi = sample(times, openings, lo), sample(times, openings, hi)
        d_lo, d_hi = sample(times, depths, lo), sample(times, depths, hi)
        per_side = 0.5 * (o_lo - o_hi) / (hi - lo)
        etch = (d_hi - d_lo) / (hi - lo)
        rows.append((lo, hi, per_side, etch, per_side / etch if etch else float("nan")))
    return rows


def krueger_reference_ratio():
    """His run-average closure/etch ratio, both from the 60 s endpoint."""
    per_side = 0.5 * KRUEGER_CLOSURE_BUDGET_NM / PROCESS_TIME_S
    etch = KRUEGER_DEPTH_NM / PROCESS_TIME_S
    return per_side, etch, per_side / etch


def top_band_balance_inversion():
    """What the near-vertical top band would need to hold open.

    The O channel and the depositor channel are both thermal/isotropic, so
    delivery cancels in their ratio and the balance is set purely by
    p_ox * J_O / (s_eff * J_dep) -- a geometry-free statement.
    """
    audit = json.loads(LIP_AUDIT.read_text())
    band = audit["bands"][0]
    j_dep = audit["source_flux_m2_s"]["depositors"]
    j_oxy = audit["source_flux_m2_s"]["oxygen"]
    p_ox = 0.0423  # Krueger converged Table 6.5 O-based polymer etch
    s_eff = band["effective_sticking"]
    o_only = p_ox * j_oxy / (s_eff * j_dep)
    required_s_eff = p_ox * j_oxy / j_dep
    fully_crosslinked = p_ox * j_oxy / (0.02 * j_dep)
    return {
        "band_nm": band["band_nm"],
        "wall_tilt_deg": band["wall_tilt_deg_mean"],
        "measured_removal_over_deposition": band["removal_over_deposition"],
        "o_channel_only": o_only,
        "ion_contribution": band["removal_over_deposition"] - o_only,
        "effective_sticking": s_eff,
        "crosslinked_fraction": band["crosslinked_fraction"],
        "required_effective_sticking_for_balance": required_s_eff,
        "published_crosslinked_row": 0.02,
        "removal_over_deposition_if_fully_crosslinked": fully_crosslinked,
    }


def main():
    print("=" * 78)
    print("EARLY-TRANSIENT ANALYSIS -- where the closure budget is spent")
    print("=" * 78)
    per_side, etch, ratio = krueger_reference_ratio()
    print(f"Krueger 60 s reference: neck {KRUEGER_NECK_NM} nm, depth "
          f"{KRUEGER_DEPTH_NM} nm")
    print(f"  closure budget {KRUEGER_CLOSURE_BUDGET_NM:.1f} nm "
          f"({per_side:.3f} nm/s per side), etch {etch:.2f} nm/s, "
          f"closure/etch = {ratio:.4f}")

    for tag in RUNS:
        traj = load_trajectory(tag)
        if traj is None:
            print(f"\n[{tag}] audit not found -- skipped")
            continue
        times, openings, depths, throats = traj
        print(f"\n--- {tag} ---")
        print(f"  final: opening {openings[-1]:.1f} nm, depth {depths[-1]:.1f} nm, "
              f"throat {throats[-1]:.0f} nm below mask top")
        print("  closure budget spent (vs Krueger's full-run 51.2 nm):")
        for t, lost, pct in closure_budget_table(times, openings):
            marker = "  <== full budget" if pct >= 100.0 and t <= 12 else ""
            print(f"    t={t:5.1f} s   lost {lost:5.1f} nm   {pct:6.1f} %{marker}")
        print("  closure/etch ratio by window (Krueger run-average "
              f"{ratio:.4f}):")
        for lo, hi, side, etch_rate, r in rate_table(times, openings, depths):
            flag = "  ~ Krueger" if 0.8 * ratio <= r <= 1.3 * ratio else ""
            print(f"    [{lo:2d},{hi:2d}] s  per-side {side:6.3f} nm/s  "
                  f"etch {etch_rate:6.2f} nm/s  ratio {r:.4f}{flag}")

    print()
    print("=" * 78)
    print("TOP-BAND BALANCE INVERSION (geometry-free: both channels thermal)")
    print("=" * 78)
    inv = top_band_balance_inversion()
    for key, value in inv.items():
        if isinstance(value, float):
            print(f"  {key:44s} {value:.4f}")
        else:
            print(f"  {key:44s} {value}")
    print()
    print("  Reading: the near-vertical lip needs s_eff = "
          f"{inv['required_effective_sticking_for_balance']:.4f} to balance, which is "
          f"{inv['required_effective_sticking_for_balance'] / 0.02:.2f}x the published "
          "crosslinked row (0.02).")
    print("  Even a FULLY crosslinked lip film reaches only removal/deposition = "
          f"{inv['removal_over_deposition_if_fully_crosslinked']:.3f}.")


if __name__ == "__main__":
    main()
