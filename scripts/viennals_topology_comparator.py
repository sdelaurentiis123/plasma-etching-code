#!/usr/bin/env python3
"""Pinned black-box ViennaLS/Petch T1 topology conformance comparator.

This adapter compiles an original API-only ViennaLS probe against declared
local upstream revisions, runs both engines on the same analytic rounded
keyhole, and compares two independent initial-value problems:

* uniform-normal closure of an initially open 100 nm neck;
* uniform-normal reopening of a prescribed 100 nm solid cap.

The reverse branch deliberately does not start from either solver's first
discrete closure.  That state has a solver-dependent subcell cap thickness and
cannot support a meaningful reverse-time comparison.  Petch remains the
authority for conservative surface-state ledgers; ViennaLS is only an external
geometry-evolution comparator.
"""
from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
from time import perf_counter

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/petch-matplotlib")

import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[1]
PROBE = ROOT / "comparators" / "viennals_keyhole_probe.cpp"
PETCH_AUDIT = ROOT / "scripts" / "topology_common_refinement_audit.py"

PINNED_REVISIONS = {
    "ViennaLS": "c3a1a2ad5bb0b05f75cd4beb71b540f6c544a287",
    "ViennaHRLE": "8aa4aa7ac4a72d101bee8a140dbf64e20963e305",
    "ViennaCore": "640411618bef44717487f25ec7417b4f6daaabf7",
    "ViennaPS": "2956ed587984c6dc38be24c6e2390e10c9b2f0a7",
}
DEFAULT_ROOTS = {
    "ViennaLS": Path("/private/tmp/ViennaLS-v5.8.3"),
    "ViennaHRLE": Path("/private/tmp/ViennaHRLE-v1.1.2"),
    "ViennaCore": Path("/private/tmp/ViennaCore-v2.2.1"),
    "ViennaPS": Path("/private/tmp/ViennaPS-audit"),
}

COAT_SPEED_UM_S = 0.025
ETCH_SPEED_UM_S = 0.050
ANALYTIC_CLOSURE_TIME_S = 2.0
ANALYTIC_REOPENING_TIME_S = 1.0


def _sha256(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def _git_revision(path):
    return subprocess.run(
        ("git", "rev-parse", "HEAD"), cwd=path, check=True,
        capture_output=True, text=True, timeout=15,
    ).stdout.strip()


def _git_status(path):
    return subprocess.run(
        ("git", "status", "--porcelain"), cwd=path, check=True,
        capture_output=True, text=True, timeout=15,
    ).stdout


def verify_upstream_roots(roots):
    receipt = {}
    for name, expected in PINNED_REVISIONS.items():
        path = Path(roots[name]).resolve()
        if not path.is_dir():
            raise FileNotFoundError(f"missing pinned {name} checkout: {path}")
        actual = _git_revision(path)
        if actual != expected:
            raise RuntimeError(
                f"{name} revision mismatch: expected {expected}, found {actual}")
        if _git_status(path):
            raise RuntimeError(f"{name} checkout is dirty: {path}")
        receipt[name] = {
            "path": path.name,
            "revision": actual,
            "worktree_clean": True,
        }
    return receipt


def parse_viennals_csv(path, dx_um):
    rows = []
    with Path(path).open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            rows.append({
                "phase": row["phase"],
                "phase_time_s": float(row["phase_time_s"]),
                "physical_time_s": float(row["physical_time_s"]),
                "components": int(row["components"]),
                "void_points": int(row["void_points"]),
                "active_points": int(row["active_points"]),
            })
    closure = next((
        row["phase_time_s"] for row in rows
        if row["phase"] == "coat" and row["void_points"] > 0), None)
    sealed = next((
        row for row in rows if row["phase"] == "sealed_reference"), None)
    reopening = next((
        row["phase_time_s"] for row in rows
        if row["phase"] == "etch" and row["void_points"] == 0), None)
    if closure is None or sealed is None or sealed["void_points"] <= 0 or reopening is None:
        raise RuntimeError(f"incomplete ViennaLS topology receipt: {path}")
    return {
        "dx_um": float(dx_um),
        "dt_s": 5.0 * float(dx_um),
        "closure_time_s": float(closure),
        "reopening_time_s": float(reopening),
        "sealed_reference_void_points": int(sealed["void_points"]),
        "history": rows,
    }


def _adjacent_refinement(levels):
    ordered = sorted(levels, key=lambda item: item["dx_um"], reverse=True)
    pairs = []
    for coarse, fine in zip(ordered[:-1], ordered[1:]):
        closure_difference = abs(
            coarse["closure_time_s"] - fine["closure_time_s"])
        reopening_difference = abs(
            coarse["reopening_time_s"] - fine["reopening_time_s"])
        pairs.append({
            "coarser_dx_um": coarse["dx_um"],
            "finer_dx_um": fine["dx_um"],
            "closure_difference_s": closure_difference,
            "closure_bound_s": coarse["dx_um"] / COAT_SPEED_UM_S,
            "reopening_difference_s": reopening_difference,
            "reopening_bound_s": coarse["dx_um"] / ETCH_SPEED_UM_S,
            "passed": bool(
                closure_difference <= coarse["dx_um"] / COAT_SPEED_UM_S
                and reopening_difference <= coarse["dx_um"] / ETCH_SPEED_UM_S),
        })
    return {
        "criterion": "adjacent event times differ by no more than one coarse-cell crossing",
        "pairs": pairs,
        "passed": bool(pairs and pairs[-1]["passed"]),
    }


def compare_levels(petch_levels, vienna_levels):
    petch_by_dx = {float(item["dx_um"]): item for item in petch_levels}
    vienna_by_dx = {float(item["dx_um"]): item for item in vienna_levels}
    if set(petch_by_dx) != set(vienna_by_dx):
        raise ValueError("Petch and ViennaLS grids do not match")
    levels = []
    for dx in sorted(petch_by_dx, reverse=True):
        petch = petch_by_dx[dx]
        vienna = vienna_by_dx[dx]
        dt = max(float(petch["dt_s"]), float(vienna["dt_s"]))
        # One engine's event is localized only to one grid cell plus one
        # checkpoint interval.  The black-box paired bound adds the independent
        # localization budgets; this is a dimensional, declared bound rather
        # than an arbitrary percentage tolerance.
        petch_closure_error = abs(
            float(petch["closure_time_s"]) - ANALYTIC_CLOSURE_TIME_S)
        vienna_closure_error = abs(
            float(vienna["closure_time_s"]) - ANALYTIC_CLOSURE_TIME_S)
        petch_reopening_error = abs(
            float(petch["reopening_time_s"]) - ANALYTIC_REOPENING_TIME_S)
        vienna_reopening_error = abs(
            float(vienna["reopening_time_s"]) - ANALYTIC_REOPENING_TIME_S)
        closure_single_bound = dx / COAT_SPEED_UM_S + dt
        reopening_single_bound = dx / ETCH_SPEED_UM_S + dt
        closure_paired_difference = abs(
            float(petch["closure_time_s"]) - float(vienna["closure_time_s"]))
        reopening_paired_difference = abs(
            float(petch["reopening_time_s"]) - float(vienna["reopening_time_s"]))
        level_passed = bool(
            petch_closure_error <= closure_single_bound
            and vienna_closure_error <= closure_single_bound
            and petch_reopening_error <= reopening_single_bound
            and vienna_reopening_error <= reopening_single_bound
            and closure_paired_difference <= 2.0 * closure_single_bound
            and reopening_paired_difference <= 2.0 * reopening_single_bound)
        levels.append({
            "dx_um": dx,
            "dt_s": dt,
            "petch_closure_time_s": float(petch["closure_time_s"]),
            "viennals_closure_time_s": float(vienna["closure_time_s"]),
            "analytic_closure_time_s": ANALYTIC_CLOSURE_TIME_S,
            "petch_closure_error_s": petch_closure_error,
            "viennals_closure_error_s": vienna_closure_error,
            "single_engine_closure_bound_s": closure_single_bound,
            "paired_closure_difference_s": closure_paired_difference,
            "paired_closure_bound_s": 2.0 * closure_single_bound,
            "petch_reopening_time_s": float(petch["reopening_time_s"]),
            "viennals_reopening_time_s": float(vienna["reopening_time_s"]),
            "analytic_reopening_time_s": ANALYTIC_REOPENING_TIME_S,
            "petch_reopening_error_s": petch_reopening_error,
            "viennals_reopening_error_s": vienna_reopening_error,
            "single_engine_reopening_bound_s": reopening_single_bound,
            "paired_reopening_difference_s": reopening_paired_difference,
            "paired_reopening_bound_s": 2.0 * reopening_single_bound,
            "passed": level_passed,
        })
    petch_refinement = _adjacent_refinement(petch_levels)
    vienna_refinement = _adjacent_refinement(vienna_levels)
    return {
        "criterion": (
            "each event agrees with the analytic fixture within one cell crossing plus "
            "one checkpoint interval; independent Petch/ViennaLS budgets add"),
        "levels": levels,
        "petch_refinement": petch_refinement,
        "viennals_refinement": vienna_refinement,
        "authoritative_level": levels[-1],
        "passed": bool(
            levels[-1]["passed"]
            and petch_refinement["passed"]
            and vienna_refinement["passed"]),
    }


def _write_json(path, payload):
    target = Path(path)
    temporary = target.with_suffix(target.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    temporary.replace(target)


def _plot(comparison, path):
    levels = sorted(comparison["levels"], key=lambda item: item["dx_um"], reverse=True)
    dx = [item["dx_um"] * 1000.0 for item in levels]
    figure, axes = plt.subplots(1, 2, figsize=(10.8, 4.1), constrained_layout=True)
    for axis, event, analytic in (
            (axes[0], "closure", ANALYTIC_CLOSURE_TIME_S),
            (axes[1], "reopening", ANALYTIC_REOPENING_TIME_S)):
        axis.plot(dx, [item[f"petch_{event}_time_s"] for item in levels],
                  "o-", label="Petch")
        axis.plot(dx, [item[f"viennals_{event}_time_s"] for item in levels],
                  "s-", label="ViennaLS v5.8.3")
        axis.axhline(analytic, color="black", linestyle=":", label="analytic")
        axis.invert_xaxis()
        axis.set_xlabel("grid spacing (nm; finer →)")
        axis.set_ylabel("event time (s)")
        axis.set_title(event.capitalize())
        axis.grid(alpha=0.25)
        axis.legend(fontsize=8)
    figure.suptitle("Matched uniform-normal keyhole topology conformance")
    figure.savefig(path, dpi=180)
    plt.close(figure)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--levels-um", type=float, nargs="+", default=(0.05, 0.025, 0.0125))
    parser.add_argument(
        "--output-dir", type=Path,
        default=ROOT / "results" / "topology_petch_viennals_t1")
    parser.add_argument("--viennals-root", type=Path, default=DEFAULT_ROOTS["ViennaLS"])
    parser.add_argument("--viennahrle-root", type=Path, default=DEFAULT_ROOTS["ViennaHRLE"])
    parser.add_argument("--viennacore-root", type=Path, default=DEFAULT_ROOTS["ViennaCore"])
    parser.add_argument("--viennaps-root", type=Path, default=DEFAULT_ROOTS["ViennaPS"])
    parser.add_argument("--compile-timeout-s", type=float, default=120.0)
    parser.add_argument("--engine-timeout-s", type=float, default=600.0)
    args = parser.parse_args()
    if (not args.levels_um or any(value <= 0.0 for value in args.levels_um)
            or args.compile_timeout_s <= 0.0 or args.engine_timeout_s <= 0.0):
        raise ValueError("grid levels and timeouts must be positive")

    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    roots = {
        "ViennaLS": args.viennals_root,
        "ViennaHRLE": args.viennahrle_root,
        "ViennaCore": args.viennacore_root,
        "ViennaPS": args.viennaps_root,
    }
    upstream = verify_upstream_roots(roots)
    compiler = os.environ.get("CXX", "clang++")
    binary = output / "viennals_keyhole_probe"
    compile_command = (
        compiler, "-std=c++17", "-O2", "-DNDEBUG",
        f"-I{Path(roots['ViennaLS']).resolve() / 'include' / 'viennals'}",
        f"-I{Path(roots['ViennaHRLE']).resolve() / 'include' / 'viennahrle'}",
        f"-I{Path(roots['ViennaCore']).resolve() / 'include' / 'viennacore'}",
        str(PROBE), "-o", str(binary),
    )
    started = perf_counter()
    subprocess.run(
        compile_command, cwd=ROOT, check=True,
        capture_output=True, text=True, timeout=args.compile_timeout_s)

    petch_dir = output / "petch"
    petch_command = (
        sys.executable, str(PETCH_AUDIT), "--levels-um",
        *(str(value) for value in args.levels_um),
        "--fixture", "rounded_keyhole", "--etch-law", "uniform_normal",
        "--output-dir", str(petch_dir),
    )
    petch_process = subprocess.run(
        petch_command, cwd=ROOT, check=False,
        capture_output=True, text=True, timeout=args.engine_timeout_s)
    if petch_process.returncode not in (0, 2):
        raise RuntimeError(
            "Petch comparator process failed:\n" + petch_process.stderr[-4000:])
    petch_audit = json.loads((petch_dir / "audit.json").read_text())
    if not all(item.get("passed", False) for item in petch_audit["levels"]):
        raise RuntimeError("Petch did not complete every bounded topology branch")

    vienna_dir = output / "viennals"
    vienna_dir.mkdir(exist_ok=True)
    vienna_levels = []
    for dx in args.levels_um:
        csv_path = vienna_dir / f"dx_{float(dx):.8g}.csv"
        process = subprocess.run(
            (str(binary), str(dx), "8", "12", str(csv_path)),
            cwd=ROOT, check=False, capture_output=True, text=True,
            timeout=args.engine_timeout_s)
        if process.returncode != 0:
            raise RuntimeError(
                f"ViennaLS dx={dx} failed ({process.returncode}):\n"
                + process.stderr[-4000:])
        vienna_levels.append(parse_viennals_csv(csv_path, dx))

    comparison = compare_levels(petch_audit["levels"], vienna_levels)
    receipt = {
        "schema": "petch.viennals_topology_comparator.v1",
        "scope": "manufactured_geometry_conformance_not_chemistry_validation",
        "passed": comparison["passed"],
        "comparison": comparison,
        "petch": {
            "git_revision": subprocess.run(
                ("git", "rev-parse", "HEAD"), cwd=ROOT, check=True,
                capture_output=True, text=True).stdout.strip(),
            "audit_script_sha256": _sha256(PETCH_AUDIT),
            "audit_relative_path": "petch/audit.json",
            "conservation_authority": True,
        },
        "vienna": {
            "upstream": upstream,
            "probe_sha256": _sha256(PROBE),
            "compiler": compiler,
            "compiler_flags": ["-std=c++17", "-O2", "-DNDEBUG"],
            "compiler_version": subprocess.run(
                (compiler, "--version"), check=True, capture_output=True,
                text=True, timeout=15).stdout.splitlines()[0],
            "geometry_comparator_only": True,
            "surface_state_conservation_authority": False,
        },
        "operator": {
            "geometry": "analytic rounded chamber plus 100 nm neck",
            "closure_initial_condition": "open neck",
            "reopening_initial_condition": "independent analytic 100 nm sealed cap",
            "normal_growth_speed_um_s": COAT_SPEED_UM_S,
            "normal_etch_speed_um_s": ETCH_SPEED_UM_S,
            "timestep_rule": "dt_s = 5 * dx_um",
        },
        "runtime": {
            "wall_time_s": perf_counter() - started,
            "compile_timeout_s": args.compile_timeout_s,
            "engine_timeout_s": args.engine_timeout_s,
            "platform": platform.platform(),
            "python": sys.version.split()[0],
        },
    }
    _write_json(output / "audit.json", receipt)
    _plot(comparison, output / "petch_viennals_topology.png")
    print(json.dumps({
        "passed": receipt["passed"],
        "wall_time_s": receipt["runtime"]["wall_time_s"],
        "authoritative_level": comparison["authoritative_level"],
    }, indent=2))
    return 0 if receipt["passed"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
