#!/usr/bin/env python3
"""Freeze the source-scoped Lam chlorine neutral-transport diagnostic."""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
from pathlib import Path
from statistics import mean, median

from petch.reactor_global import (
    MalyshevMeasuredChlorineDissociationProvider,
    malyshev_1998_chlorine_in_chlorine_diffusivity,
    malyshev_1998_eq7_transport_diagnostic,
    malyshev_1998_eq7_wall_return_inversion,
    solve_cylindrical_neutral_wall_loss,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIRECTORY = ROOT / "results" / "curated" / "reactor_global_chlorine"
CSV_PATH = OUTPUT_DIRECTORY / "malyshev_1998_lam_transport_diagnostic.csv"
AUDIT_PATH = OUTPUT_DIRECTORY / "malyshev_1998_lam_transport_diagnostic.json"
REPORT_PATH = OUTPUT_DIRECTORY / "MALYSHEV_1998_TRANSPORT_DIAGNOSTIC.md"

GAS_TEMPERATURE_K = 333.0
GAS_TEMPERATURE_BASIS = (
    "source-reported initial 60 C wall/gas condition; source says gas heats "
    "with power but does not publish the powered gas temperature"
)
PUBLISHED_FITTED_GAMMA = 0.035

FIELDNAMES = (
    "source_figure",
    "window_to_wafer_gap_cm",
    "pressure_mTorr",
    "tcp_source_power_W",
    "measured_relative_cl2_density_percent",
    "cl_to_cl2_number_density_ratio",
    "declared_gas_temperature_K",
    "gas_temperature_role",
    "gauge_pressure_Pa",
    "initial_cl2_particle_fraction",
    "eq11_particle_pressure_multiplier",
    "bulk_particle_pressure_Pa",
    "bulk_neutral_density_m3",
    "neufeld_reduced_temperature",
    "neufeld_omega_11",
    "cl_in_cl2_diffusivity_m2_s",
    "cl_in_cl2_reduced_diffusivity_m_inv_s",
    "chlorine_mean_speed_m_s",
    "required_wall_return_frequency_s_inv",
    "absorbing_wall_limit_frequency_s_inv",
    "target_is_transport_attainable",
    "effective_wall_recombination_probability",
    "matched_wall_return_frequency_s_inv",
    "published_fitted_gamma",
    "published_gamma_transport_frequency_s_inv",
    "published_gamma_replayed_relative_cl2_percent",
    "published_gamma_replay_error_percentage_point",
    "supports_prediction",
    "supports_local_wall_probability_prediction",
    "supports_wafer_flux",
    "supports_feature_depth",
)


def build_outputs() -> tuple[str, dict[str, object], str]:
    provider = (
        MalyshevMeasuredChlorineDissociationProvider.from_package_data())
    rows: list[dict[str, str]] = []
    excluded = {
        "diagnostic_flow_check": 0,
        "nonphysical_or_zero_derived_dissociation": 0,
        "electron_temperature_support_missing": 0,
        "electron_density_support_missing": 0,
    }
    effective_probabilities = []
    replay_errors = []
    unattainable_conditions = []

    for marker in provider.markers:
        if marker.validation_role == "diagnostic_flow_check":
            excluded["diagnostic_flow_check"] += 1
            continue
        if not marker.supports_eq7_inversion:
            excluded["nonphysical_or_zero_derived_dissociation"] += 1
            continue
        try:
            inversion = malyshev_1998_eq7_wall_return_inversion(marker)
        except ValueError as error:
            message = str(error)
            if "Figure-3" in message:
                excluded["electron_temperature_support_missing"] += 1
            elif "Figure-11" in message:
                excluded["electron_density_support_missing"] += 1
            else:
                raise
            continue

        diagnostic = malyshev_1998_eq7_transport_diagnostic(
            inversion,
            gas_temperature_K=GAS_TEMPERATURE_K,
            gas_temperature_basis=GAS_TEMPERATURE_BASIS,
        )
        solver_inputs = {
            "geometry": diagnostic.geometry_state.active_geometry,
            "diffusivity_m2_s": (
                diagnostic.diffusivity_state.diffusivity_m2_s),
            "mean_thermal_speed_m_s": (
                diagnostic.incident_velocity_state.mean_speed_m_s),
        }
        published_gamma_state = solve_cylindrical_neutral_wall_loss(
            **solver_inputs,
            wall_reaction_probability=PUBLISHED_FITTED_GAMMA,
        )
        predicted_relative_cl2 = 100.0 / (
            1.0
            + inversion.electron_driven_cl2_destruction_frequency_s_inv
            / (
                2.0
                * published_gamma_state.exact_loss_frequency_s_inv
            )
        )
        replay_error = (
            predicted_relative_cl2
            - marker.relative_cl2_density_percent
        )
        replay_errors.append(replay_error)

        probability = (
            diagnostic.effective_wall_recombination_probability)
        matched_frequency = None
        if probability is not None:
            effective_probabilities.append(probability)
            matched_frequency = (
                diagnostic.matched_wall_state.exact_loss_frequency_s_inv)
        else:
            unattainable_conditions.append({
                "window_to_wafer_gap_cm": marker.window_to_wafer_gap_cm,
                "pressure_mTorr": marker.pressure_mTorr,
                "tcp_source_power_W": marker.tcp_source_power_W,
                "measured_relative_cl2_density_percent": (
                    marker.relative_cl2_density_percent),
                "required_wall_return_frequency_s_inv": (
                    inversion.required_wall_return_frequency_s_inv),
                "absorbing_wall_limit_frequency_s_inv": (
                    diagnostic.absorbing_wall_state
                    .exact_loss_frequency_s_inv),
                "interpretation": (
                    "near-zero inferred dissociation makes Eq. 7 singular; "
                    "the required return exceeds the exact absorbing-wall "
                    "limit at the declared temperature"
                ),
            })

        provenance = diagnostic.diffusivity_state.provenance
        rows.append({
            "source_figure": marker.source_figure,
            "window_to_wafer_gap_cm": (
                f"{marker.window_to_wafer_gap_cm:g}"),
            "pressure_mTorr": f"{marker.pressure_mTorr:g}",
            "tcp_source_power_W": f"{marker.tcp_source_power_W:.3f}",
            "measured_relative_cl2_density_percent": (
                f"{marker.relative_cl2_density_percent:.4f}"),
            "cl_to_cl2_number_density_ratio": (
                f"{inversion.cl_to_cl2_number_density_ratio:.8e}"),
            "declared_gas_temperature_K": f"{GAS_TEMPERATURE_K:.1f}",
            "gas_temperature_role": "initial_condition_sensitivity",
            "gauge_pressure_Pa": f"{diagnostic.gauge_pressure_Pa:.8e}",
            "initial_cl2_particle_fraction": (
                f"{diagnostic.initial_cl2_particle_fraction:.8e}"),
            "eq11_particle_pressure_multiplier": (
                f"{diagnostic.particle_pressure_multiplier:.8e}"),
            "bulk_particle_pressure_Pa": (
                f"{diagnostic.bulk_particle_pressure_Pa:.8e}"),
            "bulk_neutral_density_m3": (
                f"{diagnostic.bulk_neutral_density_m3:.8e}"),
            "neufeld_reduced_temperature": (
                f"{provenance['reduced_temperature']:.12g}"),
            "neufeld_omega_11": f"{provenance['omega_11']:.12g}",
            "cl_in_cl2_diffusivity_m2_s": (
                f"{diagnostic.diffusivity_state.diffusivity_m2_s:.8e}"),
            "cl_in_cl2_reduced_diffusivity_m_inv_s": (
                f"{diagnostic.diffusivity_state.diffusivity_m2_s * diagnostic.bulk_neutral_density_m3:.8e}"  # noqa: E501
            ),
            "chlorine_mean_speed_m_s": (
                f"{diagnostic.incident_velocity_state.mean_speed_m_s:.8e}"),
            "required_wall_return_frequency_s_inv": (
                f"{inversion.required_wall_return_frequency_s_inv:.8e}"),
            "absorbing_wall_limit_frequency_s_inv": (
                f"{diagnostic.absorbing_wall_state.exact_loss_frequency_s_inv:.8e}"  # noqa: E501
            ),
            "target_is_transport_attainable": (
                str(diagnostic.target_is_transport_attainable).lower()),
            "effective_wall_recombination_probability": (
                "" if probability is None else f"{probability:.12g}"),
            "matched_wall_return_frequency_s_inv": (
                "" if matched_frequency is None
                else f"{matched_frequency:.8e}"),
            "published_fitted_gamma": f"{PUBLISHED_FITTED_GAMMA:g}",
            "published_gamma_transport_frequency_s_inv": (
                f"{published_gamma_state.exact_loss_frequency_s_inv:.8e}"),
            "published_gamma_replayed_relative_cl2_percent": (
                f"{predicted_relative_cl2:.8f}"),
            "published_gamma_replay_error_percentage_point": (
                f"{replay_error:.8f}"),
            "supports_prediction": "false",
            "supports_local_wall_probability_prediction": "false",
            "supports_wafer_flux": "false",
            "supports_feature_depth": "false",
        })

    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=FIELDNAMES, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    csv_payload = stream.getvalue()

    diffusivity_model = malyshev_1998_chlorine_in_chlorine_diffusivity()
    room_temperature_diffusivity = (
        diffusivity_model.diffusivity_cm2_s_at_pressure(
            gas_temperature_K=298.15,
            pressure_Pa=101325.0,
        )
    )
    lam_temperature_diffusivity = (
        diffusivity_model.diffusivity_cm2_s_at_pressure(
            gas_temperature_K=GAS_TEMPERATURE_K,
            pressure_Pa=101325.0,
        )
    )
    lam_temperature_state = diffusivity_model.evaluate(
        total_neutral_density_m3=(
            101325.0 / (1.380649e-23 * GAS_TEMPERATURE_K)),
        gas_temperature_K=GAS_TEMPERATURE_K,
    )
    reduced_diffusivity = (
        lam_temperature_state.diffusivity_m2_s
        * lam_temperature_state.total_neutral_density_m3
    )
    old_reduced_diffusivity = 6.21e20
    source_measurement_anchor = 0.15
    anchor_relative_residual = (
        room_temperature_diffusivity / source_measurement_anchor - 1.0)

    audit = {
        "audit_id": "MALYSHEV-1998-LAM-NEUTRAL-TRANSPORT-DIAGNOSTIC-R1",
        "claim_class": (
            "temperature-declared retrospective transport diagnostic; not "
            "a predictive reactor boundary or feature-depth result"
        ),
        "coefficient_selection_target": None,
        "reactor_fit_target": None,
        "feature_depth_target": None,
        "sources": {
            "malyshev_1998": {
                "bibkey": "malyshev-1998-lam-cl2",
                "doi": "10.1063/1.368010",
                "uses": [
                    "Cl-in-Cl2 Chapman-Enskog expression and LJ parameters",
                    "source-declared 1.25 diffusion correction",
                    "333 K initial gas/wall state",
                    "Eq. 11 particle-density adjustment with measured flows",
                    "same-board fitted gamma=0.035 replay",
                ],
            },
            "neufeld_1972": {
                "bibkey": "neufeld-1972-collision-integrals",
                "doi": "10.1063/1.1678363",
                "uses": "original Eq. 2 and Table-I Omega(1,1)* row",
                "valid_reduced_temperature": [0.3, 100.0],
                "table_ii_maximum_relative_fit_error": 0.0011,
            },
            "chantry_1987": {
                "bibkey": "chantry-1987-wall-diffusion",
                "doi": "10.1063/1.339662",
                "uses": "partial-reflection Robin boundary and exact limits",
            },
        },
        "gas_temperature": {
            "declared_K": GAS_TEMPERATURE_K,
            "basis": GAS_TEMPERATURE_BASIS,
            "measured_at_each_power": False,
            "formal_uncertainty_available": False,
        },
        "diffusivity": {
            "room_temperature_298p15K_one_atm_cm2_s": (
                room_temperature_diffusivity),
            "source_measurement_anchor_cm2_s": source_measurement_anchor,
            "room_temperature_anchor_relative_residual": (
                anchor_relative_residual),
            "source_correction_factor": 1.25,
            "333K_one_atm_cm2_s": lam_temperature_diffusivity,
            "333K_reduced_diffusivity_m_inv_s": reduced_diffusivity,
            "economou_constant_m_inv_s": old_reduced_diffusivity,
            "new_to_economou_ratio": (
                reduced_diffusivity / old_reduced_diffusivity),
            "economou_excess_relative_to_new": (
                old_reduced_diffusivity / reduced_diffusivity - 1.0),
            "supports_prediction": False,
            "reason": (
                "the source measurement anchor and LJ parameters carry no "
                "complete physical uncertainty"
            ),
            "collision_integral_method": (
                "Neufeld 1972 evaluated correlation replaces the older "
                "Hirschfelder table cited by Malyshev"
            ),
        },
        "marker_accounting": {
            "audited_marker_total": len(provider.markers),
            "successful_measured_state_rows": len(rows),
            "transport_attainable_rows": len(effective_probabilities),
            "transport_unattainable_rows": len(unattainable_conditions),
            "excluded": excluded,
        },
        "effective_probability_at_declared_temperature": {
            "minimum": min(effective_probabilities),
            "median": median(effective_probabilities),
            "maximum": max(effective_probabilities),
            "interpretation": (
                "model-conditioned values required to reproduce each marker; "
                "not local measurements, a fitted law, or predictive inputs"
            ),
        },
        "unattainable_conditions": unattainable_conditions,
        "published_gamma_forward_replay": {
            "gamma": PUBLISHED_FITTED_GAMMA,
            "gamma_role": (
                "same Figures 7-8 board fit in the source; retrospective, "
                "not independent validation"
            ),
            "row_count": len(replay_errors),
            "mean_absolute_error_percentage_point": mean(
                abs(value) for value in replay_errors),
            "root_mean_square_error_percentage_point": math.sqrt(
                mean(value ** 2 for value in replay_errors)),
            "mean_error_percentage_point": mean(replay_errors),
            "maximum_absolute_error_percentage_point": max(
                abs(value) for value in replay_errors),
            "formal_gate": None,
        },
        "cross_source_wall_comparison": {
            "status": "not_scored",
            "reason": (
                "Stafford direct conditioned-wall data are at 300 K and do "
                "not certify Malyshev's unmeasured powered gas temperature or "
                "its distributed chamber surface state"
            ),
        },
        "boundaries": {
            "incident_velocity": (
                "thermalized isotropic Maxwellian sensitivity; Guha 2008 "
                "shows freshly dissociated chlorine may be nonthermal at low "
                "pressure, and Malyshev publishes no local velocity moment"
            ),
            "rare_gas_transport": (
                "the 5% rare-gas inventory is retained in the particle-"
                "pressure ledger from the reported flows, but the source's "
                "Cl-in-Cl2 approximation is replayed rather than inventing "
                "an unvalidated multicomponent transport correction"
            ),
            "no_depth": (
                "volume-average electron state and neutral wall return do not "
                "identify species-resolved sheath-edge ion flux, IEAD, wafer "
                "neutral flux, surface kinetics, etch rate, or feature depth"
            ),
        },
        "csv_sha256": hashlib.sha256(
            csv_payload.encode("utf-8")).hexdigest(),
    }

    report = f"""# Malyshev 1998 Lam neutral-transport diagnostic

**Verdict: the missing Cl/Cl2 diffusion law is now reconstructed, but the
powered gas temperature and distributed wall state still prevent prediction.**

The source-parameterized Chapman--Enskog reconstruction, with Neufeld's
evaluated collision integral replacing the older Hirschfelder table cited by
Malyshev, gives
`D(298.15 K, 1 atm) = {room_temperature_diffusivity:.6f} cm2/s`. This is
`{100.0 * anchor_relative_residual:.2f}%` above the source's rounded
`0.15 cm2/s` room-temperature anchor after retaining its declared `1.25`
factor exactly. At the source-reported initial `333 K`, it gives
`D(1 atm) = {lam_temperature_diffusivity:.6f} cm2/s` and
`N D = {reduced_diffusivity:.6e} m-1 s-1`. The latter is
`{100.0 * (1.0 - reduced_diffusivity / old_reduced_diffusivity):.2f}%` below
the old `6.21e20 m-1 s-1` Economou constant, which was not temperature-safe.

Each of the 23 measured-state Eq.-7 rows was then mapped through the exact
cylindrical Robin eigenmode at a declared `333 K` sensitivity. The powered
particle-pressure ledger retains the reported 95/5 Cl2/rare-gas inventory;
the rare gas is not incorrectly dissociated by Eq. 11. 22 rows admit a
physical effective wall probability, spanning
`{min(effective_probabilities):.5f}--{max(effective_probabilities):.5f}` with
median `{median(effective_probabilities):.5f}`. One 11 cm / 10 mTorr /
200.704 W row is unattainable: its near-zero inferred dissociation requires a
wall-return frequency above the perfectly absorbing-wall limit. It is reported
as a physics failure, not assigned `gamma > 1`.

The source's fitted `gamma = 0.035`, replayed without refitting through the new
transport and measured electron-state board, has MAE
`{mean(abs(value) for value in replay_errors):.3f}` percentage points and RMSE
`{math.sqrt(mean(value ** 2 for value in replay_errors)):.3f}` percentage
points in relative Cl2 density. This is not independent validation: Malyshev
fit that gamma to the same Figures 7--8 data using its older electron rate and
approximate transport closure.

The principal unresolved boundary is now explicit. The paper says the gas
starts at 333 K and heats with power, but publishes no powered gas-temperature
board. The incident-speed closure is also a thermalized sensitivity; Guha's
direct measurements warn that fresh Cl can remain nonthermal at low pressure.
Stafford's direct conditioned-wall data are at 300 K and cannot be silently
extrapolated to this distributed Lam wall. Consequently the effective
probabilities are sensitivity diagnostics only; no row supports a predictive
wall law, wafer flux, etch rate, or feature depth.
"""
    return csv_payload, audit, report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    csv_payload, audit, report = build_outputs()
    json_payload = json.dumps(audit, indent=2, sort_keys=True) + "\n"
    if args.check:
        if CSV_PATH.read_text() != csv_payload:
            raise SystemExit("transport diagnostic CSV is stale")
        if AUDIT_PATH.read_text() != json_payload:
            raise SystemExit("transport diagnostic JSON is stale")
        if REPORT_PATH.read_text() != report:
            raise SystemExit("transport diagnostic report is stale")
        return
    OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)
    CSV_PATH.write_text(csv_payload)
    AUDIT_PATH.write_text(json_payload)
    REPORT_PATH.write_text(report)


if __name__ == "__main__":
    main()
