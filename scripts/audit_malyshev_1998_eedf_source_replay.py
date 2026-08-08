#!/usr/bin/env python3
"""Audit a physical Cl2 EEPF/global-model replay against Lam state markers.

The collision bytes are user supplied and hash gated; they are never copied
to the repository. The committed outputs are derived sensitivities only.
Forward TCP power is not plasma-absorbed power in Malyshev et al., so the
absorbed fractions below are a preregistered bracket, not fitted parameters.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np

from petch.reactor_global import (
    AbsorbedPowerEstimate,
    DeterministicTwoTermBoltzmannSolver,
    EEDFChlorineAbsorbedPowerModel,
    EEDFChlorineCondition,
    ElectronEnergyGrid,
    FixedPositiveIonWallEnergyProvider,
    HAMILTON_2018_CL2_STATE_CROSS_SECTIONS_SHA256,
    LeeEconomouChlorineChargedTransportProvider,
    MalyshevMeasuredElectronDensityProvider,
    MalyshevMeasuredElectronTemperatureProvider,
    PASCAL_PER_MTORR,
    PositiveIonWallEnergyState,
    ReactionNetwork,
    ReactorScalarInput,
    StateDependentChlorineNeutralTransportProvider,
    build_lee_lieberman_chlorine_particle_network,
    derive_hamilton_2018_dissociation_replay,
    load_legacy_siglo_cl2_replay,
    load_legacy_siglo_comsol_chlorine_replay,
    load_legacy_siglo_hamilton_comsol_chlorine_replay,
    lymberopoulos_economou_1995_chlorine_reduced_ion_mobilities,
    malyshev_1998_chlorine_in_chlorine_diffusivity,
    malyshev_1998_lam_geometry,
    malyshev_1998_eq11_relative_cl2_density_percent,
    stafford_2010_bounded_hill_wall_recombination_provider,
    standard_volume_flow_molecules_s,
    thermalized_chlorine_incident_velocity_state,
)
from petch.reactor_global.chlorine_lam_dissociation import (
    MalyshevMeasuredChlorineDissociationProvider,
)


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "results" / "curated" / "reactor_global_chlorine"
JSON_NAME = "malyshev_1998_eedf_source_replay.json"
CSV_NAME = "malyshev_1998_eedf_source_replay.csv"
REPORT_NAME = "MALYSHEV_1998_EEDF_SOURCE_REPLAY.md"
ATOMIC_JSON_NAME = "malyshev_1998_eedf_atomic_cl_source_replay.json"
ATOMIC_CSV_NAME = "malyshev_1998_eedf_atomic_cl_source_replay.csv"
ATOMIC_REPORT_NAME = "MALYSHEV_1998_EEDF_ATOMIC_CL_SOURCE_REPLAY.md"
HAMILTON_ATOMIC_JSON_NAME = (
    "malyshev_1998_eedf_hamilton_atomic_cl_source_replay.json")
HAMILTON_ATOMIC_CSV_NAME = (
    "malyshev_1998_eedf_hamilton_atomic_cl_source_replay.csv")
HAMILTON_ATOMIC_REPORT_NAME = (
    "MALYSHEV_1998_EEDF_HAMILTON_ATOMIC_CL_SOURCE_REPLAY.md")
HAMILTON_JSON_NAME = "malyshev_1998_eedf_hamilton_source_replay.json"
HAMILTON_CSV_NAME = "malyshev_1998_eedf_hamilton_source_replay.csv"
HAMILTON_REPORT_NAME = "MALYSHEV_1998_EEDF_HAMILTON_SOURCE_REPLAY.md"
HAMILTON_ATOMIC_EE_JSON_NAME = (
    "malyshev_1998_eedf_hamilton_atomic_cl_ee_source_replay.json")
HAMILTON_ATOMIC_EE_CSV_NAME = (
    "malyshev_1998_eedf_hamilton_atomic_cl_ee_source_replay.csv")
HAMILTON_ATOMIC_EE_REPORT_NAME = (
    "MALYSHEV_1998_EEDF_HAMILTON_ATOMIC_CL_EE_SOURCE_REPLAY.md")
ABSORBED_FRACTIONS = (0.30, 0.50, 0.70)
SOURCE_POWERS_W = (300.0, 500.0)
GAP_CM = 6.5
PRESSURE_MTORR = 2.0
CHLORINE_PARTIAL_PRESSURE_FRACTION = 0.95
GAS_TEMPERATURE_K = 333.0
CL2_FLOW_SCCM = 114.0
STANDARD_TEMPERATURE_K = 273.15
STANDARD_PRESSURE_PA = 101325.0
ION_TEMPERATURE_EV = 0.12
SOURCE_FREQUENCY_HZ = 13.56e6
ION_WALL_ENERGY_EV = {"Cl2+": 12.0, "Cl+": 14.0}
DIRECT_WALL_RATIO_DOMAIN = (0.105610, 0.779646)
DECLARED_WALL_RATIO_DOMAIN = (1.0e-5, 30.0)


def _scalar(value: float, unit: str, source: str) -> ReactorScalarInput:
    return ReactorScalarInput(
        value=value,
        unit=unit,
        source=source,
        evidence_kind="sensitivity",
        relative_uncertainty=None,
    )


def _condition(source_power_W: float, absorbed_fraction: float):
    geometry_state = malyshev_1998_lam_geometry(GAP_CM)
    absorbed_power_W = source_power_W * absorbed_fraction
    return EEDFChlorineCondition(
        condition_id=(
            f"malyshev-gap{GAP_CM:g}-p{PRESSURE_MTORR:g}-"
            f"source{source_power_W:g}-absorbed{absorbed_fraction:.2f}"
        ),
        geometry=geometry_state.active_geometry,
        neutral_control_volume=geometry_state.neutral_control_volume,
        pressure=_scalar(
            PRESSURE_MTORR * PASCAL_PER_MTORR
            * CHLORINE_PARTIAL_PRESSURE_FRACTION,
            "Pa",
            "Malyshev 2 mTorr total pressure times printed 95% Cl2 feed",
        ),
        gas_temperature=_scalar(
            GAS_TEMPERATURE_K,
            "K",
            "Malyshev printed 333 K initial gas/wall temperature; powered "
            "gas temperature is unpublished",
        ),
        chlorine_molecule_feed=_scalar(
            standard_volume_flow_molecules_s(
                CL2_FLOW_SCCM,
                standard_temperature_K=STANDARD_TEMPERATURE_K,
                standard_pressure_Pa=STANDARD_PRESSURE_PA,
            ),
            "molecule s^-1",
            "Malyshev printed 114 sccm Cl2 at 2 mTorr; standard state declared",
        ),
        source_power=_scalar(
            source_power_W,
            "W",
            "Malyshev forward TCP power into the matching network",
        ),
        absorbed_power=AbsorbedPowerEstimate(
            lower_W=absorbed_power_W,
            upper_W=absorbed_power_W,
            point_W=absorbed_power_W,
            boundary_kind="preregistered_source_power_fraction_sensitivity",
            measurement_source=(
                "Malyshev forward TCP power into matching network"),
            loss_source=(
                f"assumed plasma-absorbed fraction {absorbed_fraction:.2f}; "
                "hardware/plasma losses not published"),
            measurement_evidence="measured",
            loss_evidence="sensitivity",
            provenance={
                "absorbed_fraction": absorbed_fraction,
                "coefficient_selection_target": None,
                "feature_depth_used": False,
            },
        ),
        # At 2 mTorr, omega/N is large enough that the 13.56 MHz field-heating
        # term is strongly collision limited.  The DC replay's 600 Td ceiling
        # is therefore not a physically meaningful RF ceiling; 10 kTd covers
        # the preregistered high-frequency local-field sensitivity without
        # constraining the power solution at its bound.
        reduced_field_bounds_Td=(20.0, 10_000.0),
        source_frequency=_scalar(
            SOURCE_FREQUENCY_HZ,
            "Hz",
            "Malyshev printed 13.56 MHz TCP source frequency",
        ),
        neutral_density_constraint=(
            "chlorine_nuclei_equivalent_molecules"),
    )


def _providers():
    wall = stafford_2010_bounded_hill_wall_recombination_provider(
        "anodized_aluminum",
        valid_cl_to_cl2_ratio=DECLARED_WALL_RATIO_DOMAIN,
        valid_gas_temperature_K=(300.0, GAS_TEMPERATURE_K),
        transfer_source=(
            "declared transfer to Malyshev anodized-Al ratio and 333 K "
            "initial wall state; no reactor, feature, or depth fit"),
    )
    neutral = StateDependentChlorineNeutralTransportProvider(
        wall_recombination_provider=wall,
        incident_velocity_state=thermalized_chlorine_incident_velocity_state(
            GAS_TEMPERATURE_K,
            source=(
                "333 K isotropic thermalized-Cl sensitivity; Malyshev "
                "publishes no powered incident velocity moment"),
            evidence_kind="sensitivity",
            relative_uncertainty=None,
            provenance={"coefficient_selection_target": None},
        ),
        diffusivity_model=malyshev_1998_chlorine_in_chlorine_diffusivity(),
    )
    all_mobilities = (
        lymberopoulos_economou_1995_chlorine_reduced_ion_mobilities())
    charged = LeeEconomouChlorineChargedTransportProvider(
        reduced_mobilities={
            species: all_mobilities[species] for species in ("Cl2+", "Cl+")
        },
        ion_temperature=_scalar(
            ION_TEMPERATURE_EV,
            "eV",
            "Lymberopoulos-Economou 1995 chlorine source-model value",
        ),
    )
    wall_energy = FixedPositiveIonWallEnergyProvider(
        PositiveIonWallEnergyState(
            energy_eV_per_lost_ion=ION_WALL_ENERGY_EV,
            source=(
                "fixed Cl2+/Cl+ sheath-wall energy sensitivity; no Lam "
                "species-resolved IED is published"),
            evidence_kind="sensitivity",
            relative_uncertainty=None,
        )
    )
    return charged, neutral, wall_energy


def _reference(source_power_W: float) -> dict[str, object]:
    temperature = (
        MalyshevMeasuredElectronTemperatureProvider.from_package_data()
        .evaluate(
            window_to_wafer_gap_cm=GAP_CM,
            pressure_mTorr=PRESSURE_MTORR,
            tcp_source_power_W=source_power_W,
        )
    )
    density = (
        MalyshevMeasuredElectronDensityProvider.from_package_data().evaluate(
            window_to_wafer_gap_cm=GAP_CM,
            pressure_mTorr=PRESSURE_MTORR,
            tcp_source_power_W=source_power_W,
        )
    )
    dissociation_markers = tuple(
        marker
        for marker in (
            MalyshevMeasuredChlorineDissociationProvider.from_package_data()
            .markers
        )
        if marker.window_to_wafer_gap_cm == GAP_CM
        and marker.pressure_mTorr == PRESSURE_MTORR
        and abs(marker.tcp_source_power_W - source_power_W)
        <= marker.digitization_power_uncertainty_W
        and marker.validation_role == "reactor_dissociation_validation_candidate"
    )
    if len(dissociation_markers) > 1:
        raise RuntimeError("ambiguous Malyshev dissociation marker")
    marker = dissociation_markers[0] if dissociation_markers else None
    return {
        "electron_temperature_eV": temperature.electron_temperature.value,
        "electron_temperature_method": temperature.method,
        "electron_density_m3": density.volume_average_electron_density.value,
        "electron_density_method": density.method,
        "relative_cl2_density_percent": (
            None if marker is None else marker.relative_cl2_density_percent),
        "relative_cl2_digitization_uncertainty_percentage_point": (
            None if marker is None else
            marker.digitization_relative_cl2_uncertainty_percentage_point),
        "reported_absolute_cl2_density_accuracy_percent": (
            None if marker is None else
            marker.reported_absolute_density_relative_uncertainty_percent),
    }


def run_replay(
    collision_deck: Path,
    *,
    energy_cells: int,
    atomic_cl_momentum: Path | None = None,
    hamilton_dissociation: bool = False,
    electron_electron_collisions: bool = False,
    absorbed_fractions: tuple[float, ...] = ABSORBED_FRACTIONS,
    source_powers_W: tuple[float, ...] = SOURCE_POWERS_W,
) -> dict[str, object]:
    if atomic_cl_momentum is None:
        replay = load_legacy_siglo_cl2_replay(
            collision_deck, maximum_energy_eV=200.0)
        if hamilton_dissociation:
            replay = derive_hamilton_2018_dissociation_replay(replay)
            model_variant = "legacy_siglo_hamilton_molecular_cl2_only"
        else:
            model_variant = "legacy_siglo_molecular_cl2_only"
        molecular_replay = replay
        atomic_momentum_hash = None
    else:
        if hamilton_dissociation:
            replay = load_legacy_siglo_hamilton_comsol_chlorine_replay(
                collision_deck,
                atomic_cl_momentum,
                maximum_energy_eV=200.0,
            )
            model_variant = (
                "legacy_siglo_hamilton_plus_comsol_nist_atomic_cl")
        else:
            replay = load_legacy_siglo_comsol_chlorine_replay(
                collision_deck,
                atomic_cl_momentum,
                maximum_energy_eV=200.0,
            )
            model_variant = "legacy_siglo_plus_comsol_nist_atomic_cl"
        molecular_replay = replay.molecular_replay
        atomic_momentum_hash = replay.atomic_momentum_payload_sha256
    if int(energy_cells) != energy_cells or energy_cells < 100:
        raise ValueError("energy_cells must be an integer of at least 100")
    segment_weights = np.asarray((0.20, 0.30, 0.25, 0.15, 0.10))
    segment_cells = np.floor(segment_weights * energy_cells).astype(int)
    segment_cells[-1] += int(energy_cells) - int(np.sum(segment_cells))
    thresholds = tuple(sorted({
        float(process.energy_loss_eV)
        for process in replay.derived_deck.processes
        if process.energy_loss_eV is not None
        and 0.0 < float(process.energy_loss_eV) < 200.0
    }))
    energy_grid = ElectronEnergyGrid.piecewise_linear(
        (0.0, 0.5, 5.0, 20.0, 80.0, 200.0),
        tuple(int(value) for value in segment_cells),
        inserted_boundaries_eV=thresholds,
    )
    lee = build_lee_lieberman_chlorine_particle_network()
    heavy = ReactionNetwork(species=lee.species, reactions=lee.reactions[6:8])
    model = EEDFChlorineAbsorbedPowerModel(
        DeterministicTwoTermBoltzmannSolver(
            energy_grid,
            replay.derived_deck,
        ),
        replay.collision_chemistry,
        heavy,
        electron_electron_coulomb_model=(
            "isotropic_classical_debye"
            if electron_electron_collisions else "none"),
    )
    if electron_electron_collisions:
        model_variant += "_plus_isotropic_ee"
    charged, neutral, wall_energy = _providers()
    absorbed_fractions = tuple(float(value) for value in absorbed_fractions)
    source_powers_W = tuple(float(value) for value in source_powers_W)
    if (
        not absorbed_fractions
        or any(value not in ABSORBED_FRACTIONS for value in absorbed_fractions)
        or not source_powers_W
        or any(value not in SOURCE_POWERS_W for value in source_powers_W)
    ):
        raise ValueError("requested audit subset is outside the preregistration")
    rows = []
    solver_failures = []
    reference_seeds = {}
    # The middle sensitivity is a numerical continuation anchor only. It was
    # declared before the board and is not selected from any observable.
    solve_order = tuple(
        value for value in (0.50, 0.30, 0.70)
        if value in absorbed_fractions
    )
    for fraction in solve_order:
        previous = None
        previous_exhaust = None
        previous_field = 300.0
        for source_power_W in source_powers_W:
            if previous is None and source_power_W in reference_seeds:
                seed = reference_seeds[source_power_W]
                previous = seed["densities"]
                previous_exhaust = seed["exhaust"]
                previous_field = seed["field"]
            condition = _condition(source_power_W, fraction)
            try:
                solution = model.solve(
                    condition,
                    charged_transport_provider=charged,
                    neutral_wall_transport_provider=neutral,
                    wall_energy_provider=wall_energy,
                    initial_densities_m3=previous,
                    initial_exhaust_loss_frequency_s_inv=previous_exhaust,
                    initial_reduced_electric_field_Td=previous_field,
                    residual_tolerance=2.0e-7,
                    maximum_evaluations=1200,
                    maximum_tail_population_fraction=1.0e-6,
                )
            except RuntimeError as error:
                if atomic_cl_momentum is not None:
                    raise
                solver_failures.append({
                    "absorbed_fraction_sensitivity": fraction,
                    "source_power_W": source_power_W,
                    "validation_role": (
                        "reactor_diagnostic_training_candidate"
                        if source_power_W == 300.0
                        else "held_out_reactor_diagnostic"),
                    "failure_class": type(error).__name__,
                    "failure_message": str(error),
                    "physical_interpretation": (
                        "molecular-only electron collision deck cannot close "
                        "a dissociated chlorine state; atomic-Cl electron "
                        "collisions are mandatory"),
                })
                previous = None
                previous_exhaust = None
                previous_field = 300.0
                continue
            previous = dict(solution.densities_m3)
            previous_exhaust = solution.exhaust_loss_frequency_s_inv
            previous_field = solution.reduced_electric_field_Td
            if fraction == 0.50:
                reference_seeds[source_power_W] = {
                    "densities": previous,
                    "exhaust": previous_exhaust,
                    "field": previous_field,
                }
            reference = _reference(source_power_W)
            densities = solution.densities_m3
            ratio = densities["Cl"] / densities["Cl2"]
            mean_energy_temperature_proxy = (
                2.0 / 3.0 * solution.mean_electron_energy_eV)
            modeled_cl2_particle_fraction = (
                100.0 * densities["Cl2"]
                / (densities["Cl2"] + densities["Cl"])
            )
            modeled_relative_cl2 = (
                malyshev_1998_eq11_relative_cl2_density_percent(
                    densities["Cl2"], densities["Cl"])
            )
            chlorine_nuclei_equivalent_density = (
                densities["Cl2"] + 0.5 * densities["Cl"])
            chlorine_particle_density = densities["Cl2"] + densities["Cl"]
            reference_temperature = float(
                reference["electron_temperature_eV"])
            reference_density = float(reference["electron_density_m3"])
            reference_relative_cl2 = reference["relative_cl2_density_percent"]
            reported_cl2_accuracy = reference[
                "reported_absolute_cl2_density_accuracy_percent"]
            relative_cl2_error_percent = (
                None if reference_relative_cl2 is None else 100.0 * (
                    modeled_relative_cl2 / float(reference_relative_cl2)
                    - 1.0
                )
            )
            rows.append({
                "absorbed_fraction_sensitivity": fraction,
                "source_power_W": source_power_W,
                "absorbed_power_W": source_power_W * fraction,
                "validation_role": (
                    "reactor_diagnostic_training_candidate"
                    if source_power_W == 300.0
                    else "held_out_reactor_diagnostic"
                ),
                "reduced_electric_field_Td": (
                    solution.reduced_electric_field_Td),
                "mean_electron_energy_eV": solution.mean_electron_energy_eV,
                "mean_energy_temperature_proxy_eV": (
                    mean_energy_temperature_proxy),
                "measured_oes_electron_temperature_eV": (
                    reference_temperature),
                "temperature_proxy_percent_error": 100.0 * (
                    mean_energy_temperature_proxy - reference_temperature
                ) / reference_temperature,
                "electron_density_m3": densities["e"],
                "measured_volume_average_electron_density_m3": (
                    reference_density),
                "electron_density_percent_error": 100.0 * (
                    densities["e"] - reference_density) / reference_density,
                "electronegativity": densities["Cl-"] / densities["e"],
                "cl_to_cl2_ratio": ratio,
                "chlorine_nuclei_equivalent_molecule_density_m3": (
                    chlorine_nuclei_equivalent_density),
                "target_chlorine_nuclei_equivalent_molecule_density_m3": (
                    condition.target_neutral_density_m3),
                "chlorine_particle_density_m3": chlorine_particle_density,
                "chlorine_particle_density_multiplier_vs_gauge_equivalent": (
                    chlorine_particle_density
                    / chlorine_nuclei_equivalent_density),
                "wall_ratio_inside_direct_marker_domain": (
                    DIRECT_WALL_RATIO_DOMAIN[0]
                    <= ratio <= DIRECT_WALL_RATIO_DOMAIN[1]),
                "modeled_cl2_particle_fraction_percent": (
                    modeled_cl2_particle_fraction),
                "modeled_relative_cl2_density_percent_eq11": (
                    modeled_relative_cl2),
                "measured_relative_cl2_density_percent": (
                    reference_relative_cl2),
                "relative_cl2_eq11_error_percentage_point": (
                    None if reference_relative_cl2 is None else
                    modeled_relative_cl2 - float(reference_relative_cl2)
                ),
                "relative_cl2_eq11_error_percent": relative_cl2_error_percent,
                "reported_absolute_cl2_density_accuracy_percent": (
                    reported_cl2_accuracy),
                "within_reported_cl2_density_accuracy": (
                    None
                    if relative_cl2_error_percent is None
                    or reported_cl2_accuracy is None
                    else abs(relative_cl2_error_percent)
                    <= float(reported_cl2_accuracy)
                ),
                "cl2plus_axial_flux_m2_s": (
                    solution.axial_positive_ion_flux_m2_s["Cl2+"]),
                "clplus_axial_flux_m2_s": (
                    solution.axial_positive_ion_flux_m2_s["Cl+"]),
                "total_positive_ion_axial_flux_m2_s": sum(
                    solution.axial_positive_ion_flux_m2_s.values()),
                "clplus_ion_fraction": (
                    solution.axial_positive_ion_flux_m2_s["Cl+"]
                    / sum(solution.axial_positive_ion_flux_m2_s.values())),
                "modeled_power_density_W_m3": (
                    solution.modeled_power_density_W_m3),
                "absorbed_power_density_W_m3": (
                    solution.absorbed_power_density_W_m3),
                "maximum_normalized_residual": (
                    solution.maximum_normalized_residual),
                "relative_electron_growth_closure": (
                    solution.collision_chemistry_state
                    .relative_electron_growth_closure),
                "electron_to_neutral_density_ratio": (
                    densities["e"]
                    / sum(
                        densities[name]
                        for name in replay.derived_deck.targets)
                ),
                "electron_electron_coulomb_model": (
                    solution.electron_solution
                    .electron_electron_coulomb_model),
                "coulomb_logarithm": (
                    solution.electron_solution.coulomb_logarithm),
                "coulomb_nonlinear_iterations": (
                    solution.electron_solution.nonlinear_iteration_count),
                "coulomb_nonlinear_weighted_residual": (
                    solution.electron_solution.nonlinear_weighted_residual),
                "electron_growth_root_evaluations": (
                    solution.electron_solution.growth_root_evaluations),
                "solver_evaluations": solution.solver_evaluations,
                "supports_reactor_state_prediction": False,
                "supports_wafer_flux": False,
                "supports_feature_depth": False,
            })
    rows.sort(key=lambda row: (
        row["absorbed_fraction_sensitivity"], row["source_power_W"]))
    return {
        "schema": "petch.malyshev_1998_eedf_source_replay.v2",
        "claim_class": (
            "physical source replay and preregistered sensitivity; not "
            "independent reactor validation"),
        "model_variant": model_variant,
        "raw_collision_payload_sha256": (
            molecular_replay.raw_payload_sha256),
        "atomic_momentum_payload_sha256": atomic_momentum_hash,
        "hamilton_state_cross_sections_sha256": (
            HAMILTON_2018_CL2_STATE_CROSS_SECTIONS_SHA256
            if hamilton_dissociation else None),
        "derived_collision_deck_sha256": replay.derived_deck.payload_sha256,
        "raw_collision_bytes_committed": False,
        "energy_grid": {
            "family": "threshold_aligned_piecewise_linear_v1",
            "breakpoints_eV": [0.0, 0.5, 5.0, 20.0, 80.0, 200.0],
            "nominal_cell_count": energy_cells,
            "segment_cell_counts": [int(value) for value in segment_cells],
            "inserted_collision_thresholds_eV": list(thresholds),
            "actual_cell_count": energy_grid.cell_count,
        },
        "preregistration": {
            "absorbed_fraction_sensitivities": list(ABSORBED_FRACTIONS),
            "executed_absorbed_fractions": list(absorbed_fractions),
            "executed_source_powers_W": list(source_powers_W),
            "reactor_diagnostic_training_candidate_power_W": 300.0,
            "held_out_reactor_diagnostic_power_W": 500.0,
            "fraction_selected_in_this_run": None,
            "feature_depth_used_for_selection": False,
            "krueger_825_nm_used_for_selection": False,
        },
        "condition": {
            "window_to_wafer_gap_cm": GAP_CM,
            "total_pressure_mTorr": PRESSURE_MTORR,
            "modeled_chlorine_partial_pressure_fraction": (
                CHLORINE_PARTIAL_PRESSURE_FRACTION),
            "cl2_flow_sccm": CL2_FLOW_SCCM,
            "gas_temperature_K_sensitivity": GAS_TEMPERATURE_K,
            "direct_wall_ratio_domain": list(DIRECT_WALL_RATIO_DOMAIN),
            "bounded_wall_ratio_sensitivity_domain": list(
                DECLARED_WALL_RATIO_DOMAIN),
            "ion_temperature_eV": ION_TEMPERATURE_EV,
            "source_frequency_Hz": SOURCE_FREQUENCY_HZ,
            "ion_wall_energy_eV_sensitivity": ION_WALL_ENERGY_EV,
        },
        "missing_physics": list(replay.missing_reactor_channels) + [
            "five_percent_rare_gas_collision_ion_and_transport_channels",
            "powered_condition_gas_temperature",
        ],
        "electron_electron_coulomb_model": (
            "isotropic_classical_debye"
            if electron_electron_collisions else "none"),
        "comparison_boundaries": {
            "rf_heating": (
                "Hagelaar-Pitchford high-frequency two-term operator at "
                "13.56 MHz with modern RMS E/N convention; this is a "
                "quasi-stationary local-field sensitivity, not a "
                "time-periodic or nonlocal ICP electron-heating solution"
            ),
            "coulomb": (
                "Hagelaar-Pitchford 2005 isotropic electron-electron "
                "Fokker-Planck term with classical Debye Coulomb logarithm; "
                "electron-ion and anisotropic momentum terms remain absent"
                if electron_electron_collisions else
                "electron-electron and electron-ion Coulomb collisions absent"
            ),
            "temperature": (
                "2/3 mean E is a diagnostic proxy, not the exact Malyshev "
                "OES observable"),
            "wall": (
                "bounded Langmuir-Hill coverage response fitted only to "
                "Stafford direct spinning-wall markers; the asymptote is the "
                "largest direct marker and transfer beyond ratio 0.78 and "
                "from 300 to 333 K remains sensitivity evidence"),
            "relative_cl2": (
                "Malyshev Eq. 11 observable 100*nCl2/(nCl2+nCl/2); "
                "the distinct particle fraction nCl2/(nCl2+nCl) is also "
                "reported but never compared to the measurement. The "
                "source's about +/-25% absolute-density accuracy is used as "
                "a bounded accuracy band, never relabeled as one sigma"),
            "wafer_flux": (
                "axial global-model flux is not a validated local wafer "
                "flux or IED/IAD"),
        },
        "supports_reactor_state_prediction": False,
        "supports_wafer_flux": False,
        "supports_feature_depth": False,
        "rows": rows,
        "solver_failures": solver_failures,
        "all_requested_conditions_closed": not solver_failures,
    }


def _write(result: dict[str, object], output_directory: Path) -> None:
    output_directory.mkdir(parents=True, exist_ok=True)
    atomic = result["atomic_momentum_payload_sha256"] is not None
    hamilton = "hamilton" in result["model_variant"]
    electron_electron = "isotropic_ee" in result["model_variant"]
    if hamilton and atomic and electron_electron:
        json_name = HAMILTON_ATOMIC_EE_JSON_NAME
        csv_name = HAMILTON_ATOMIC_EE_CSV_NAME
        report_name = HAMILTON_ATOMIC_EE_REPORT_NAME
    elif hamilton and atomic:
        json_name = HAMILTON_ATOMIC_JSON_NAME
        csv_name = HAMILTON_ATOMIC_CSV_NAME
        report_name = HAMILTON_ATOMIC_REPORT_NAME
    elif hamilton:
        json_name = HAMILTON_JSON_NAME
        csv_name = HAMILTON_CSV_NAME
        report_name = HAMILTON_REPORT_NAME
    else:
        json_name = ATOMIC_JSON_NAME if atomic else JSON_NAME
        csv_name = ATOMIC_CSV_NAME if atomic else CSV_NAME
        report_name = ATOMIC_REPORT_NAME if atomic else REPORT_NAME
    (output_directory / json_name).write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    rows = result["rows"]
    csv_rows = rows if rows else result["solver_failures"]
    with (output_directory / csv_name).open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=list(csv_rows[0]),
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(csv_rows)
    if not rows:
        failure_lines = [
            "| {absorbed_fraction_sensitivity:.2f} | {source_power_W:.0f} | "
            "{failure_class} | {physical_interpretation} |".format(**row)
            for row in result["solver_failures"]
        ]
        report = f"""# Malyshev 1998 molecular-only EEPF failure receipt

## Verdict

The molecular-only collision deck is **retired under the corrected Malyshev
Eq. 11 pressure closure**. It creates a dissociated chlorine population but
contains no atomic-Cl electron collision target, so the coupled particle and
charge system does not close. This failure is recorded rather than preserving
the stale source-replay numbers. The atomic-Cl replay is the minimum operative
deck.

| absorbed fraction | source W | failure | interpretation |
|---:|---:|---|---|
{chr(10).join(failure_lines)}

No reactor state, wafer flux, or feature-depth claim is supported by this
negative board. No feature observable selected any coefficient.
"""
        (output_directory / report_name).write_text(report, encoding="utf-8")
        return
    table_lines = []
    for row in rows:
        cl2_error = row["relative_cl2_eq11_error_percentage_point"]
        cl2_pass = row["within_reported_cl2_density_accuracy"]
        table_lines.append(
            "| {absorbed_fraction_sensitivity:.2f} | {source_power_W:.0f} | "
            "{reduced_electric_field_Td:.1f} | "
            "{electron_density_percent_error:+.1f}% | "
            "{temperature_proxy_percent_error:+.1f}% | "
            f"{'n/a' if cl2_error is None else f'{cl2_error:+.1f} pp'} | "
            f"{'n/a' if cl2_pass is None else ('PASS' if cl2_pass else 'MISS')} | "
            "{total_positive_ion_axial_flux_m2_s:.3e} | "
            "{maximum_normalized_residual:.1e} |".format(**row)
        )
    variant_text = str(result["model_variant"]).replace("_", " ")
    atomic_boundary = (
        "- Atomic-Cl ionization is included, but electron detachment from "
        "Cl- and tracked excited-state kinetics remain absent.\n"
        if atomic else
        "- Atomic-Cl ionization, electron detachment from Cl-, and tracked "
        "excited-state kinetics are absent.\n"
    )
    coulomb_boundary = (
        "- Isotropic electron-electron Coulomb drift/diffusion is included "
        "with the classical Debye logarithm. Electron-ion and anisotropic "
        "electron-electron momentum terms remain absent; this is a "
        "density-coupling sensitivity, not a complete Coulomb model.\n"
        if electron_electron else
        "- Electron-electron and electron-ion Coulomb collisions are "
        "absent.\n"
    )
    report = f"""# Malyshev 1998 deterministic Cl2 EEPF source replay

## Verdict

This board is a **physical source replay and sensitivity, not a validated
knobs-to-wafer model**. Its collision variant is: **{variant_text}**. It solves
the non-Maxwellian two-term EEPF together with particle and power balances and
never selects a coefficient from feature depth. Forward TCP power is not
measured absorbed plasma power, so 30%, 50%, and 70% are all reported rather
than optimized.

| absorbed fraction | source W | E/N Td | ne error | 2/3 mean-E proxy error vs OES | Eq.-11 Cl2 error | within reported Cl2 accuracy | axial positive-ion flux m-2 s-1 | max closure |
|---:|---:|---:|---:|---:|---:|:---:|---:|---:|
{chr(10).join(table_lines)}

## Use boundary

- The 300 W condition is only a future reactor-diagnostic training candidate;
  500 W is reserved as a held-out reactor diagnostic. No fraction is selected
  by this run.
- The 13.56 MHz electron-heating term uses the Hagelaar--Pitchford
  high-frequency two-term equation and the modern RMS-field convention. At
  this intermediate RF frequency it is a declared quasi-stationary
  local-field sensitivity, not a time-periodic or spatially nonlocal ICP
  heating solution.
- `2/3 <E>` is not the exact OES forward observable. Its error diagnoses EEPF
  shape/chemistry but is not an apples-to-apples temperature validation.
- Malyshev's about +/-25% Cl2 absolute-density accuracy is used as a reported
  accuracy band, not a statistical sigma; digitization error remains separate.
- The bounded Stafford coverage fit is trained only on the direct Cl/Cl2
  interval `{DIRECT_WALL_RATIO_DOMAIN}`. Its transfer to ratio domain
  `{DECLARED_WALL_RATIO_DOMAIN}` and from 300 to 333 K is sensitivity evidence;
  pressure, power, material, and the direct-marker provenance remain explicit.
{atomic_boundary}{coulomb_boundary}- The global axial flux is not yet a local wafer flux,
  and it carries no species-resolved IED/IAD.
- Raw collision bytes are not committed. Replay identity is the SHA-256 in the
  JSON receipt.

Consequently this board cannot support feature depth. It exists to locate the
reactor-state failure before any feature coupling is attempted.
"""
    (output_directory / report_name).write_text(report, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("collision_deck", type=Path)
    parser.add_argument("--atomic-cl-momentum", type=Path)
    parser.add_argument("--hamilton-dissociation", action="store_true")
    parser.add_argument("--electron-electron-collisions", action="store_true")
    parser.add_argument("--energy-cells", type=int, default=400)
    parser.add_argument(
        "--absorbed-fraction",
        type=float,
        choices=ABSORBED_FRACTIONS,
    )
    parser.add_argument(
        "--source-power-W",
        type=float,
        choices=SOURCE_POWERS_W,
    )
    parser.add_argument("--output-directory", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--no-write", action="store_true")
    parser.add_argument("--summary-only", action="store_true")
    arguments = parser.parse_args()
    result = run_replay(
        arguments.collision_deck,
        energy_cells=arguments.energy_cells,
        atomic_cl_momentum=arguments.atomic_cl_momentum,
        hamilton_dissociation=arguments.hamilton_dissociation,
        electron_electron_collisions=(
            arguments.electron_electron_collisions),
        absorbed_fractions=(
            ABSORBED_FRACTIONS
            if arguments.absorbed_fraction is None
            else (arguments.absorbed_fraction,)
        ),
        source_powers_W=(
            SOURCE_POWERS_W
            if arguments.source_power_W is None
            else (arguments.source_power_W,)
        ),
    )
    if not arguments.no_write:
        _write(result, arguments.output_directory)
    if arguments.summary_only:
        for row in result["rows"]:
            cl2_error = row["relative_cl2_eq11_error_percentage_point"]
            cl2_pass = row["within_reported_cl2_density_accuracy"]
            print(
                f"f={row['absorbed_fraction_sensitivity']:.2f} "
                f"P={row['source_power_W']:.0f} W "
                f"E/N={row['reduced_electric_field_Td']:.1f} Td "
                f"ne={row['electron_density_percent_error']:+.1f}% "
                f"energy={row['temperature_proxy_percent_error']:+.1f}% "
                f"Cl2={'n/a' if cl2_error is None else f'{cl2_error:+.1f} pp'} "
                f"gate={'n/a' if cl2_pass is None else ('PASS' if cl2_pass else 'MISS')}"
            )
        for failure in result["solver_failures"]:
            print(
                f"f={failure['absorbed_fraction_sensitivity']:.2f} "
                f"P={failure['source_power_W']:.0f} W "
                f"FAIL={failure['failure_class']}"
            )
    else:
        print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
