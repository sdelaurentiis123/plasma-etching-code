#!/usr/bin/env python3
"""Grade Lan--Jeon C4F6 drift against an independent BOLSIG+ replay."""
from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import io
import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
MEASUREMENT_CSV = (
    ROOT / "data" / "experimental" / "lan_jeon_2014_c4f6"
    / "figure7_pure_c4f6_drift.csv"
)
LOCAL_AUDIT = ROOT / "results" / "curated" / "c4f6_electron_swarm_v1" / "audit.json"
DEFAULT_OUTPUT = ROOT / "results" / "curated" / "c4f6_bolsig_bulk_reference_v1"
BOLSIG_ZIP_SHA256 = "60a98cf2a4d67a3ff212fd30e0072eff95b2c1295eeeb2e8ae8d94f7ebde5f42"
BOLSIG_BINARY_SHA256 = "fe9d3c995e5b71033eb7db8c562b8d191e30aabaa8549ba0438de994607b8dc0"
BOLSIG_RAW_OUTPUT_SHA256 = "2df9619acc40910109d5dfe368f7a134936a0fc4659788c895aaf8238b35603a"
EXPORTED_COLLISION_SHA256 = "1403f7977be4c06379a097b4df60369c36539a3e59633f6b1a9c513348d208d0"
EXPORTED_INSTRUCTION_SHA256 = "9ce147ccc5d90994368bc1d1201c3b1a0da9ab1ea5fe75afdaec3e654899b3cd"
CSV_FIELDS = (
    "reduced_electric_field_Td",
    "measured_legacy_pt_Wv_m_s",
    "petch_flux_drift_m_s",
    "bolsig_flux_drift_m_s",
    "bolsig_bulk_drift_m_s",
    "bolsig_mean_energy_eV",
    "bolsig_flux_mobility_times_N_m_inv_V_inv_s_inv",
    "bolsig_bulk_mobility_times_N_m_inv_V_inv_s_inv",
)


def _digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def parse_bolsig_transport(path: Path) -> tuple[dict[str, float], ...]:
    if _digest(path) != BOLSIG_RAW_OUTPUT_SHA256:
        raise RuntimeError("BOLSIG+ raw comparator output checksum changed")
    text = path.read_text(encoding="utf-8", errors="strict")
    if "BOLSIG+ version: 03/2016" not in text or "Growth model" not in text:
        raise RuntimeError("unexpected BOLSIG+ output header")
    section = text.split("Transport coefficients", 1)[1].split("Rate coefficients", 1)[0]
    parsed = []
    for line in section.splitlines():
        fields = line.split()
        if len(fields) != 24 or not fields[0].isdigit():
            continue
        parsed.append({
            "printed_field_Td": float(fields[1]),
            "mean_energy_eV": float(fields[2]),
            "flux_mobility_times_N": float(fields[3]),
            "bulk_mobility_times_N": float(fields[4]),
        })
    if len(parsed) != 18:
        raise RuntimeError("BOLSIG+ comparator row topology changed")
    return tuple(parsed)


def _measurement_rows() -> tuple[dict[str, str], ...]:
    with MEASUREMENT_CSV.open(newline="", encoding="utf-8") as stream:
        rows = tuple(csv.DictReader(stream))
    if len(rows) != 18:
        raise RuntimeError("Lan--Jeon drift board topology changed")
    return rows


def build_rows(bolsig_rows) -> tuple[dict[str, float], ...]:
    measured = _measurement_rows()
    local = json.loads(LOCAL_AUDIT.read_text(encoding="utf-8"))
    local_flux = local["transport_definition_diagnostic"][
        "predicted_flux_drift_velocity_m_s"
    ]
    rows = []
    for source, reference, petch_flux in zip(measured, bolsig_rows, local_flux):
        field = float(source["reduced_electric_field_Td"])
        if abs(reference["printed_field_Td"] / field - 1.0) > 5.0e-4:
            raise RuntimeError("BOLSIG+ field differs from source board")
        flux_mobility = reference["flux_mobility_times_N"]
        bulk_mobility = reference["bulk_mobility_times_N"]
        rows.append({
            "reduced_electric_field_Td": field,
            "measured_legacy_pt_Wv_m_s": float(source["drift_velocity_m_s"]),
            "petch_flux_drift_m_s": float(petch_flux),
            "bolsig_flux_drift_m_s": flux_mobility * field * 1.0e-21,
            "bolsig_bulk_drift_m_s": bulk_mobility * field * 1.0e-21,
            "bolsig_mean_energy_eV": reference["mean_energy_eV"],
            "bolsig_flux_mobility_times_N_m_inv_V_inv_s_inv": flux_mobility,
            "bolsig_bulk_mobility_times_N_m_inv_V_inv_s_inv": bulk_mobility,
        })
    return tuple(rows)


def _summary(predicted, measured) -> dict:
    relative = np.asarray(predicted, dtype=float) / np.asarray(measured, dtype=float) - 1.0
    return {
        "mean_signed_relative_residual": float(np.mean(relative)),
        "mean_absolute_relative_residual": float(np.mean(np.abs(relative))),
        "maximum_absolute_relative_residual": float(np.max(np.abs(relative))),
        "signed_relative_residual": relative.tolist(),
    }


def audit(rows) -> dict:
    measured = [row["measured_legacy_pt_Wv_m_s"] for row in rows]
    petch_flux = [row["petch_flux_drift_m_s"] for row in rows]
    bolsig_flux = [row["bolsig_flux_drift_m_s"] for row in rows]
    bolsig_bulk = [row["bolsig_bulk_drift_m_s"] for row in rows]
    cross_solver = np.asarray(bolsig_flux) / np.asarray(petch_flux) - 1.0
    return {
        "schema": "petch.c4f6-bolsig-bulk-reference-audit.v1",
        "sources": {
            "lan_jeon_doi": "10.3938/jkps.64.1320",
            "pt_definition_doi": "10.1088/1361-6595/abe729",
            "pt_definition_bibkey": "casey-2021-pt-foundations",
            "collision_deck_bibkey": "lan-jeon-2014-c4f6",
            "bolsig_method_doi": "10.1088/0963-0252/14/4/011",
            "bolsig_version": "03/2016",
            "bolsig_download_url": "https://www.bolsig.laplace.univ-tlse.fr/wp-content/uploads/2016/03/bolsigplus032016-mac.zip",
            "bolsig_zip_sha256": BOLSIG_ZIP_SHA256,
            "bolsig_binary_sha256": BOLSIG_BINARY_SHA256,
            "bolsig_raw_output_sha256": BOLSIG_RAW_OUTPUT_SHA256,
            "exported_collision_sha256": EXPORTED_COLLISION_SHA256,
            "exported_instruction_sha256": EXPORTED_INSTRUCTION_SHA256,
        },
        "conditions": {
            "growth_model": "BOLSIG+ density-gradient expansion",
            "energy_grid": "800-point manual parabolic, 0--200 eV",
            "gas_temperature_K": 300.0,
            "mole_fraction": {"C4F6": 1.0},
            "field_count": len(rows),
            "collision_coefficients_retuned": False,
        },
        "cross_solver_flux_agreement": {
            "maximum_absolute_relative_difference": float(np.max(np.abs(cross_solver))),
            "mean_absolute_relative_difference": float(np.mean(np.abs(cross_solver))),
            "interpretation": "independent deterministic two-term flux replay",
        },
        "legacy_pt_Wv_comparison": {
            "petch_flux": _summary(petch_flux, measured),
            "bolsig_flux": _summary(bolsig_flux, measured),
            "bolsig_bulk": _summary(bolsig_bulk, measured),
        },
        "verdict": {
            "flux_replay_independently_corroborated": True,
            "density_gradient_bulk_resolves_legacy_Wv": False,
            "reason": (
                "Lan--Jeon reports legacy pulsed-Townsend average drift Wv. "
                "It is neither flux drift nor the universal bulk coefficient. "
                "The casey-2021-pt-foundations transformation requires the same-study PT "
                "Townsend property and longitudinal diffusion, which are not "
                "reported in the Lan--Jeon source."
            ),
            "supports_unique_c4f6_reactor_state": False,
            "supports_wafer_flux": False,
            "supports_krueger_depth": False,
        },
    }


def _csv_text(rows) -> str:
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=CSV_FIELDS, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return stream.getvalue()


def _read_committed_rows(path: Path) -> tuple[dict[str, float], ...]:
    with path.open(newline="", encoding="utf-8") as stream:
        records = tuple(csv.DictReader(stream))
    if len(records) != 18 or tuple(records[0]) != CSV_FIELDS:
        raise RuntimeError("committed BOLSIG+ comparator CSV topology changed")
    return tuple({key: float(value) for key, value in row.items()} for row in records)


def _readme(result: dict) -> str:
    cross = result["cross_solver_flux_agreement"]
    comparison = result["legacy_pt_Wv_comparison"]
    return f"""# C4F6 BOLSIG+ bulk-transport reference

The exact Lan--Jeon collision tables were exported without retuning and run in
the official deterministic BOLSIG+ `03/2016` density-gradient mode at all 18
Figure-7 fields. This is an independent comparator, not a production
dependency.

BOLSIG+ flux drift agrees with petch flux drift to at most
`{100.0 * cross['maximum_absolute_relative_difference']:.3f}%` across the board.
That independently corroborates the local two-term implementation.

Changing the comparison quantity from flux to BOLSIG+ bulk drift does **not**
reproduce the plotted legacy PT average drift `Wv`: the mean absolute residual
is `{100.0 * comparison['bolsig_bulk']['mean_absolute_relative_residual']:.2f}%`
and the maximum is
`{100.0 * comparison['bolsig_bulk']['maximum_absolute_relative_residual']:.2f}%`,
versus `{100.0 * comparison['bolsig_flux']['mean_absolute_relative_residual']:.2f}%`
for the independent flux replay. Bulk drift overcorrects the low-field,
attachment-dominated points and improves much of the high-field board.

This is a physical observable-definition result, not a numerical failure.
`casey-2021-pt-foundations` shows that the older PT transport property transforms as
`W_B = W_tilde + alpha_tilde D_L_tilde`. Lan--Jeon reports `Wv` but not the
same-study Townsend property and longitudinal diffusion needed for that
transformation. The source cross sections must not be retuned to force either
flux or modern bulk drift through a legacy quantity with missing co-observables.

The result validates the collision solver's local flux calculation. It does
not identify a C4F6 reactor state, wafer flux, or Krueger depth.
"""


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bolsig-output", type=Path)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    csv_path = args.output / "bolsig_transport.csv"
    if args.bolsig_output is None:
        if not csv_path.exists():
            raise SystemExit("--bolsig-output is required before the comparator CSV exists")
        rows = _read_committed_rows(csv_path)
    else:
        rows = build_rows(parse_bolsig_transport(args.bolsig_output))
    result = audit(rows)
    payloads = {
        "bolsig_transport.csv": _csv_text(rows),
        "audit.json": json.dumps(result, indent=2, sort_keys=True) + "\n",
        "README.md": _readme(result),
    }
    if args.check:
        for name, text in payloads.items():
            target = args.output / name
            if not target.exists() or target.read_text(encoding="utf-8") != text:
                raise SystemExit(f"committed C4F6 BOLSIG+ audit is stale: {target}")
        return
    args.output.mkdir(parents=True, exist_ok=True)
    for name, text in payloads.items():
        (args.output / name).write_text(text, encoding="utf-8")
    print(json.dumps(result["legacy_pt_Wv_comparison"], indent=2))


if __name__ == "__main__":
    main()
