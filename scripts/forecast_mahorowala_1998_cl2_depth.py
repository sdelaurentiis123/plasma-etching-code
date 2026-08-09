#!/usr/bin/env python3
"""Project Mahorowala's blind Cl2 depths from reactor diagnostics only.

The reactor/equipment transfer is deliberately source-estimate conditioned:
the thesis's center ion-current and Cl/ion estimates set wafer normalization,
not any etch rate.  All eleven Table-2.2 depths remain untouched until the
projection is complete.  The resulting surface-plane calculation is the next
failure-localization rung, not yet a feature-depth prediction: Cl2+ surface
kinetics and deterministic in-feature transport remain separate gates.
"""
from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
from pathlib import Path

import numpy as np

from petch.chang_sawin_chlorine_si import (
    BaloochCl2IonSiMechanism,
    ChangSawinClIonSiMechanism,
)
from petch.chlorine_species_resolved_si import (
    SpeciesResolvedChlorineSiMechanism,
)
from petch.product_ion_chlorine_si import (
    LeeChangProductIonSiSurfaceSensitivity,
    PRODUCT_ION_MASS_AMU,
)
from petch.reactor_global import (
    AbsorbedPowerEstimate,
    ChlorineWaferTransferEvidence,
    DeterministicTwoTermBoltzmannSolver,
    DiagnosticConditionedChlorineWaferTransfer,
    DiagnosticConditionedEtchProductResidenceTransfer,
    EtchProductPlasmaCondition,
    DiagnosticConditionedRFSheathTransfer,
    EEDFChlorineAbsorbedPowerModel,
    EEDFChlorineCondition,
    EEDFChlorineFixedFeedback,
    ElectronEnergyGrid,
    FixedPositiveIonWallEnergyProvider,
    LeeEtchProductLinearReactor,
    LeeEconomouChlorineChargedTransportProvider,
    PASCAL_PER_MTORR,
    PositiveIonWallEnergyState,
    ReactionNetwork,
    ReactorScalarInput,
    StateDependentChlorineNeutralTransportProvider,
    build_lee_lieberman_chlorine_particle_network,
    load_legacy_siglo_hamilton_comsol_chlorine_replay,
    lee_1995_reactive_product_wall,
    lee_1995_reflective_product_wall,
    lymberopoulos_economou_1995_chlorine_reduced_ion_mobilities,
    malyshev_1998_chlorine_in_chlorine_diffusivity,
    malyshev_1998_lam_geometry,
    stafford_2010_bounded_hill_wall_recombination_provider,
    standard_volume_flow_molecules_s,
    thermalized_chlorine_incident_velocity_state,
)
from petch.reactor_global.geometry import ElectropositiveEdgeFactors
from petch.sheath import bohm_speed
from petch.surface_kinetics import EnergeticFlux, SurfaceFluxes


ROOT = Path(__file__).resolve().parents[1]
SOURCE_CSV = (
    ROOT / "data" / "experimental" / "mahorowala_1998_cl2"
    / "table2_2_oxide_mask_fixed_time.csv"
)
SOURCE_MANIFEST = (
    ROOT / "data" / "experimental" / "mahorowala_1998_cl2"
    / "audit_manifest.json"
)
CL2PLUS_SURFACE_MANIFEST = (
    ROOT / "data" / "experimental" / "chang_1998_figure5_7"
    / "digitization_manifest.json"
)
CL2PLUS_SURFACE_CSV = (
    ROOT / "data" / "experimental" / "chang_1998_figure5_7"
    / "cl2plus_poly_si_yield_fit.csv"
)
DEFAULT_OUTPUT = (
    ROOT / "results" / "curated" / "reactor_to_feature_chlorine"
)
JSON_NAME = "mahorowala_1998_diagnostic_conditioned_depth_projection.json"
CSV_NAME = "mahorowala_1998_diagnostic_conditioned_depth_projection.csv"
REPORT_NAME = "MAHOROWALA_1998_DIAGNOSTIC_CONDITIONED_DEPTH_PROJECTION.md"

ABSORBED_SOURCE_POWER_FRACTION = 0.36464679405196565
GAP_CM_EQUIPMENT_CLASS_PRIOR = 6.5
NEUTRAL_CONTROL_VOLUME_M3 = 0.043
PRESSURE_MTORR = 10.0
GAS_TEMPERATURE_K = 500.0
SOURCE_FREQUENCY_HZ = 13.56e6
ION_TEMPERATURE_EV = 0.12
STANDARD_TEMPERATURE_K = 273.15
STANDARD_PRESSURE_PA = 101325.0
ELECTRODE_AREA_M2 = 0.04
WAFER_DIAMETER_M = 0.200
EXPOSED_POLYSILICON_AREA_FRACTION = 0.60
REFERENCE_ION_CURRENT_A_CM2 = 2.0e-3
REFERENCE_NEUTRAL_TO_ION_RATIO = 100.0
REFERENCE_SICL2_TO_TOTAL_ION_RATIO = 0.30
PLASMA_POTENTIAL_EV = 20.0
ION_WALL_ENERGY_EV = {"Cl2+": 12.0, "Cl+": 14.0}
E_CHARGE_C = 1.602176634e-19


def _hash(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _scalar(value: float, unit: str, source: str) -> ReactorScalarInput:
    return ReactorScalarInput(
        value=value,
        unit=unit,
        source=source,
        evidence_kind="sensitivity",
        relative_uncertainty=None,
    )


def _condition(source_power_W: float, flow_sccm: float):
    geometry = malyshev_1998_lam_geometry(
        GAP_CM_EQUIPMENT_CLASS_PRIOR)
    absorbed_power_W = source_power_W * ABSORBED_SOURCE_POWER_FRACTION
    return EEDFChlorineCondition(
        condition_id=(
            f"mahorowala-source{source_power_W:g}-flow{flow_sccm:g}"
        ),
        geometry=geometry.active_geometry,
        neutral_control_volume=geometry.neutral_control_volume,
        pressure=_scalar(
            PRESSURE_MTORR * PASCAL_PER_MTORR,
            "Pa",
            "Mahorowala thesis Chapter 2: pure Cl2 at fixed 10 mTorr",
        ),
        gas_temperature=_scalar(
            GAS_TEMPERATURE_K,
            "K",
            "Mahorowala thesis Chapter 4 center-condition estimate",
        ),
        chlorine_molecule_feed=_scalar(
            standard_volume_flow_molecules_s(
                flow_sccm,
                standard_temperature_K=STANDARD_TEMPERATURE_K,
                standard_pressure_Pa=STANDARD_PRESSURE_PA,
            ),
            "molecule s^-1",
            "Mahorowala Table 2.2 Cl2 flow; standard state declared",
        ),
        source_power=_scalar(
            source_power_W,
            "W",
            "Mahorowala Table 2.2 inductive-power setpoint",
        ),
        absorbed_power=AbsorbedPowerEstimate(
            lower_W=absorbed_power_W,
            upper_W=absorbed_power_W,
            point_W=absorbed_power_W,
            boundary_kind="lam_alliance_diagnostic_conditioned_transfer",
            measurement_source="Mahorowala inductive-power setpoint",
            loss_source=(
                "constant fraction conditioned only on Malyshev Lam Alliance "
                "electron density; transferred across Lam platforms as sensitivity"
            ),
            measurement_evidence="measured",
            loss_evidence="sensitivity",
            provenance={
                "absorbed_fraction": ABSORBED_SOURCE_POWER_FRACTION,
                "source_equipment": "Lam Alliance",
                "target_equipment": "Lam TCP 9400SE",
                "etch_rate_or_depth_used": False,
            },
        ),
        reduced_field_bounds_Td=(20.0, 10_000.0),
        source_frequency=_scalar(
            SOURCE_FREQUENCY_HZ,
            "Hz",
            "13.56 MHz Lam ICP equipment-class assumption; thesis omits frequency",
        ),
        neutral_density_constraint="total_neutral_particles",
    )


def _providers():
    wall = stafford_2010_bounded_hill_wall_recombination_provider(
        "anodized_aluminum",
        valid_cl_to_cl2_ratio=(1.0e-5, 30.0),
        valid_gas_temperature_K=(300.0, GAS_TEMPERATURE_K),
        transfer_source=(
            "bounded Stafford anodized-Al response transferred to the "
            "Mahorowala 500 K source estimate; no reactor state or depth fit"
        ),
    )
    neutral = StateDependentChlorineNeutralTransportProvider(
        wall_recombination_provider=wall,
        incident_velocity_state=thermalized_chlorine_incident_velocity_state(
            GAS_TEMPERATURE_K,
            source="Mahorowala 500 K thermalized-Cl source estimate",
            evidence_kind="sensitivity",
            relative_uncertainty=None,
            provenance={"etch_rate_or_depth_used": False},
        ),
        diffusivity_model=malyshev_1998_chlorine_in_chlorine_diffusivity(),
    )
    mobilities = lymberopoulos_economou_1995_chlorine_reduced_ion_mobilities()
    charged = LeeEconomouChlorineChargedTransportProvider(
        reduced_mobilities={
            name: mobilities[name] for name in ("Cl2+", "Cl+")
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
            source="species-resolved plasma-wall energy sensitivity",
            evidence_kind="sensitivity",
            relative_uncertainty=None,
        )
    )
    return charged, neutral, wall_energy


def _model(collision_deck: Path, atomic_cl_momentum: Path, energy_cells: int):
    replay = load_legacy_siglo_hamilton_comsol_chlorine_replay(
        collision_deck,
        atomic_cl_momentum,
        maximum_energy_eV=200.0,
    )
    weights = np.asarray((0.20, 0.30, 0.25, 0.15, 0.10))
    cells = np.floor(weights * int(energy_cells)).astype(int)
    cells[-1] += int(energy_cells) - int(np.sum(cells))
    thresholds = tuple(sorted({
        float(process.energy_loss_eV)
        for process in replay.derived_deck.processes
        if process.energy_loss_eV is not None
        and 0.0 < float(process.energy_loss_eV) < 200.0
    }))
    grid = ElectronEnergyGrid.piecewise_linear(
        (0.0, 0.5, 5.0, 20.0, 80.0, 200.0),
        tuple(int(value) for value in cells),
        inserted_boundaries_eV=thresholds,
    )
    lee = build_lee_lieberman_chlorine_particle_network()
    heavy = ReactionNetwork(species=lee.species, reactions=lee.reactions[6:8])
    model = EEDFChlorineAbsorbedPowerModel(
        DeterministicTwoTermBoltzmannSolver(grid, replay.derived_deck),
        replay.collision_chemistry,
        heavy,
    )
    return model, replay, grid


def _source_rows():
    with SOURCE_CSV.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def _solve_states(model, conditions):
    charged, neutral, wall_energy = _providers()
    center = (400.0, 100.0)
    order = sorted(
        conditions,
        key=lambda item: (
            item != center,
            abs(item[0] - center[0]) / 150.0
            + abs(item[1] - center[1]) / 75.0,
        ),
    )
    states = {}
    previous = None
    previous_exhaust = None
    previous_field = 500.0
    for source_power_W, flow_sccm in order:
        solution = model.solve(
            _condition(source_power_W, flow_sccm),
            charged_transport_provider=charged,
            neutral_wall_transport_provider=neutral,
            wall_energy_provider=wall_energy,
            initial_densities_m3=(
                None if previous is None else previous.densities_m3),
            initial_exhaust_loss_frequency_s_inv=(
                None if previous is None else previous_exhaust),
            initial_reduced_electric_field_Td=previous_field,
            residual_tolerance=2.0e-7,
            maximum_evaluations=1200,
            maximum_tail_population_fraction=1.0e-6,
        )
        states[(source_power_W, flow_sccm)] = solution
        previous = solution
        previous_exhaust = solution.exhaust_loss_frequency_s_inv
        previous_field = solution.reduced_electric_field_Td
    return states


def _surface_projection(
        mechanism, boundary, ion_flux_m2_s, *, energy_eV=None, weight=None):
    if energy_eV is None:
        energy_eV = np.asarray([boundary.mean_impact_energy_eV])
    if weight is None:
        weight = np.asarray([1.0])
    fluxes = SurfaceFluxes(
        {"Cl": boundary.atomic_chlorine_flux_m2_s},
        (
            EnergeticFlux(
                "Cl+",
                ion_flux_m2_s,
                np.asarray(energy_eV, dtype=float),
                np.ones(np.asarray(energy_eV).shape),
                np.asarray(weight, dtype=float),
            ),
        ),
    )
    result = mechanism.advance(
        mechanism.initial_state(), fluxes, 75.0, strict=False)
    return result


def _cl2plus_surface_projection(
        mechanism, boundary, surface_chlorination_fraction, *,
        energy_eV=None, weight=None, strict=True):
    if energy_eV is None:
        energy_eV = np.asarray([boundary.mean_impact_energy_eV])
    if weight is None:
        weight = np.asarray([1.0])
    fluxes = SurfaceFluxes(
        {"Cl": boundary.atomic_chlorine_flux_m2_s},
        (
            EnergeticFlux(
                "Cl2+",
                boundary.positive_ion_flux_m2_s["Cl2+"],
                np.asarray(energy_eV, dtype=float),
                np.ones(np.asarray(energy_eV).shape),
                np.asarray(weight, dtype=float),
            ),
        ),
    )
    return mechanism.advance(
        mechanism.initial_state(),
        fluxes,
        75.0,
        surface_chlorination_fraction=surface_chlorination_fraction,
        strict=strict,
    )


def _combined_surface_projection(
        mechanism, boundary, cl_iead, cl2_iead, *, sicl2_flux_m2_s=0.0):
    fluxes = SurfaceFluxes(
        {
            "Cl": boundary.atomic_chlorine_flux_m2_s,
            "SiCl2": float(sicl2_flux_m2_s),
        },
        (
            EnergeticFlux(
                "Cl+",
                boundary.positive_ion_flux_m2_s["Cl+"],
                cl_iead.energy_eV,
                np.ones(cl_iead.energy_eV.shape),
                cl_iead.weight,
            ),
            EnergeticFlux(
                "Cl2+",
                boundary.positive_ion_flux_m2_s["Cl2+"],
                cl2_iead.energy_eV,
                np.ones(cl2_iead.energy_eV.shape),
                cl2_iead.weight,
            ),
        ),
    )
    return mechanism.advance(
        mechanism.initial_state(), fluxes, 75.0, strict=False)


def _solve_product_surface_fixed_point(
        mechanism, boundary, cl_iead, cl2_iead, transfer, exhaust_frequency):
    """Solve Gamma_SiCl2 = residence_transfer(predicted Si source)."""

    def evaluate(product_flux):
        surface = _combined_surface_projection(
            mechanism,
            boundary,
            cl_iead,
            cl2_iead,
            sicl2_flux_m2_s=product_flux,
        )
        gross_rate = float(
            surface.clplus_removal_rate_si_m2_s
            + surface.cl2plus_removal_rate_si_m2_s
        )
        residence = transfer.predict(
            gross_si_source_rate_m2_s=gross_rate,
            exhaust_loss_frequency_s_inv=exhaust_frequency,
        )
        return surface, residence, product_flux - residence.sicl2_flux_m2_s

    lower = 0.0
    lower_surface, lower_residence, lower_residual = evaluate(lower)
    upper = max(
        transfer.reference_total_ion_flux_m2_s,
        lower_residence.sicl2_flux_m2_s,
        1.0,
    )
    upper_surface, upper_residence, upper_residual = evaluate(upper)
    for _ in range(64):
        if upper_residual >= 0.0:
            break
        upper *= 2.0
        upper_surface, upper_residence, upper_residual = evaluate(upper)
    else:  # pragma: no cover - monotone source closure should always bracket
        raise RuntimeError("failed to bracket SiCl2 residence fixed point")
    if lower_residual > 0.0 or upper_residual < 0.0:
        raise RuntimeError("invalid SiCl2 residence fixed-point bracket")
    surface = upper_surface
    residence = upper_residence
    for _ in range(96):
        midpoint = 0.5 * (lower + upper)
        surface, residence, residual = evaluate(midpoint)
        if residual > 0.0:
            upper = midpoint
        else:
            lower = midpoint
        if upper - lower <= 2.0e-12 * max(midpoint, 1.0):
            break
    product_flux = 0.5 * (lower + upper)
    surface, residence, residual = evaluate(product_flux)
    return surface, residence, residual


def _solve_table4_product_surface_fixed_point(
        mechanism, boundary, cl_iead, cl2_iead, product_reactor,
        product_condition):
    """Close predicted Si removal through the full Table-IV product network."""
    exposed_area = (
        np.pi * (0.5 * WAFER_DIAMETER_M) ** 2
        * EXPOSED_POLYSILICON_AREA_FRACTION
    )

    def evaluate(product_flux):
        surface = _combined_surface_projection(
            mechanism,
            boundary,
            cl_iead,
            cl2_iead,
            sicl2_flux_m2_s=product_flux,
        )
        gross_rate = float(
            surface.clplus_removal_rate_si_m2_s
            + surface.cl2plus_removal_rate_si_m2_s
        )
        plasma = product_reactor.solve(
            product_condition,
            gross_si_removal_flux_m2_s=gross_rate,
            exposed_silicon_area_m2=exposed_area,
        )
        return (
            surface,
            plasma,
            product_flux - plasma.wafer_neutral_flux_m2_s["SiCl2"],
        )

    lower = 0.0
    lower_surface, lower_plasma, lower_residual = evaluate(lower)
    upper = max(
        boundary.total_positive_ion_flux_m2_s,
        lower_plasma.wafer_neutral_flux_m2_s["SiCl2"],
        1.0,
    )
    upper_surface, upper_plasma, upper_residual = evaluate(upper)
    for _ in range(64):
        if upper_residual >= 0.0:
            break
        upper *= 2.0
        upper_surface, upper_plasma, upper_residual = evaluate(upper)
    else:  # pragma: no cover
        raise RuntimeError("failed to bracket Table-IV product fixed point")
    if lower_residual > 0.0 or upper_residual < 0.0:
        raise RuntimeError("invalid Table-IV product fixed-point bracket")
    surface = upper_surface
    plasma = upper_plasma
    for _ in range(96):
        midpoint = 0.5 * (lower + upper)
        surface, plasma, residual = evaluate(midpoint)
        if residual > 0.0:
            upper = midpoint
        else:
            lower = midpoint
        if upper - lower <= 2.0e-12 * max(midpoint, 1.0):
            break
    product_flux = 0.5 * (lower + upper)
    surface, plasma, residual = evaluate(product_flux)
    return surface, plasma, residual


def _blend_product_feedback(previous, predicted, fraction):
    fraction = float(fraction)
    if not 0.0 < fraction <= 1.0:
        raise ValueError("feedback relaxation fraction must be in (0, 1]")
    return EEDFChlorineFixedFeedback(
        chlorine_species_source_m3_s={
            name: (
                (1.0 - fraction)
                * previous.chlorine_species_source_m3_s[name]
                + fraction
                * predicted.chlorine_species_source_m3_s[name]
            )
            for name in ("Cl2", "Cl")
        },
        extra_neutral_density_m3=(
            (1.0 - fraction) * previous.extra_neutral_density_m3
            + fraction * predicted.extra_neutral_density_m3
        ),
        extra_positive_charge_density_m3=(
            (1.0 - fraction) * previous.extra_positive_charge_density_m3
            + fraction * predicted.extra_positive_charge_density_m3
        ),
        extra_collisional_power_density_W_m3=(
            (1.0 - fraction)
            * previous.extra_collisional_power_density_W_m3
            + fraction
            * predicted.extra_collisional_power_density_W_m3
        ),
        extra_charged_wall_power_density_W_m3=(
            (1.0 - fraction)
            * previous.extra_charged_wall_power_density_W_m3
            + fraction
            * predicted.extra_charged_wall_power_density_W_m3
        ),
        source=(
            "under-relaxed deterministic Lee Table-IV block coupling; "
            "relaxation changes convergence path, not the fixed point"
        ),
        supports_prediction=False,
    )


def _feedback_relative_distance(left, right):
    left_values = np.asarray((
        left.chlorine_species_source_m3_s["Cl2"],
        left.chlorine_species_source_m3_s["Cl"],
        left.extra_neutral_density_m3,
        left.extra_positive_charge_density_m3,
        left.extra_collisional_power_density_W_m3,
        left.extra_charged_wall_power_density_W_m3,
    ))
    right_values = np.asarray((
        right.chlorine_species_source_m3_s["Cl2"],
        right.chlorine_species_source_m3_s["Cl"],
        right.extra_neutral_density_m3,
        right.extra_positive_charge_density_m3,
        right.extra_collisional_power_density_W_m3,
        right.extra_charged_wall_power_density_W_m3,
    ))
    scale = np.maximum(np.maximum(np.abs(left_values), np.abs(right_values)), 1.0)
    return float(np.max(np.abs(left_values - right_values) / scale))


def _table4_product_condition(
        *, state, condition, charged_provider, wall_boundary):
    equivalent_temperature = (2.0 / 3.0) * state.mean_electron_energy_eV
    transport_condition = condition.transport_condition(equivalent_temperature)
    charged_state = charged_provider.predict(
        transport_condition, state.densities_m3)
    reference_transport = charged_state.positive_ion_transport["Cl+"]
    reference_bohm_speed = bohm_speed(equivalent_temperature, 35.453)
    common_edge = ElectropositiveEdgeFactors(
        axial=(
            reference_transport.axial_flux_velocity_m_s
            / reference_bohm_speed
        ),
        radial=(
            reference_transport.radial_flux_velocity_m_s
            / reference_bohm_speed
        ),
    )
    return EtchProductPlasmaCondition(
        geometry=transport_condition.geometry,
        neutral_control_volume_m3=(
            transport_condition.neutral_control_volume.value),
        electron_density_m3=state.densities_m3["e"],
        chlorine_atom_density_m3=state.densities_m3["Cl"],
        chlorine_negative_ion_density_m3=state.densities_m3["Cl-"],
        electron_temperature_eV=equivalent_temperature,
        gas_temperature_K=GAS_TEMPERATURE_K,
        exhaust_loss_frequency_s_inv=state.exhaust_loss_frequency_s_inv,
        common_edge_factors=common_edge,
        wall_boundary=wall_boundary,
        source=(
            "Mahorowala 200 mm wafer and 60% exposed poly-Si; "
            "Lee-Graves-Lieberman Table-IV product network and "
            "common-positive-ion edge-factor closure"
        ),
    )


def _solve_coupled_table4_product_state(
        *, model, condition, initial_state, transfer, sheath_transfer,
        combined_mechanism, product_reactor, product_wall, charged_provider,
        neutral_provider, wall_energy_provider, applied_bias_power_W,
        maximum_iterations=16, feedback_relaxation=0.4,
        coupling_tolerance=5.0e-4, initial_feedback=None,
        initial_product_ion_fluxes=None):
    """Close product pressure, charge, chlorine, and minimum-power feedback."""
    state = initial_state
    feedback = (
        EEDFChlorineFixedFeedback.zero()
        if initial_feedback is None else initial_feedback
    )
    if not isinstance(feedback, EEDFChlorineFixedFeedback):
        raise TypeError("initial_feedback must be EEDFChlorineFixedFeedback")
    product_ion_fluxes = (
        {name: 0.0 for name in PRODUCT_ION_MASS_AMU}
        if initial_product_ion_fluxes is None
        else {
            str(name): float(value)
            for name, value in initial_product_ion_fluxes.items()
        }
    )
    if (
        set(product_ion_fluxes) != set(PRODUCT_ION_MASS_AMU)
        or any(not np.isfinite(value) or value < 0.0
               for value in product_ion_fluxes.values())
    ):
        raise ValueError("invalid initial product-ion flux mapping")
    final = None
    coupling_residual = np.inf
    for iteration in range(1, int(maximum_iterations) + 1):
        full_boundary = transfer.predict(
            model_positive_ion_flux_m2_s=state.axial_positive_ion_flux_m2_s,
            model_atomic_chlorine_density_m3=state.densities_m3["Cl"],
            gas_temperature_K=GAS_TEMPERATURE_K,
            applied_bias_power_W=applied_bias_power_W,
        )
        product_ion_flux = sum(product_ion_fluxes.values())
        chlorine_power_fraction = (
            full_boundary.total_positive_ion_flux_m2_s
            / (
                full_boundary.total_positive_ion_flux_m2_s
                + product_ion_flux
            )
        )
        boundary = transfer.predict(
            model_positive_ion_flux_m2_s=state.axial_positive_ion_flux_m2_s,
            model_atomic_chlorine_density_m3=state.densities_m3["Cl"],
            gas_temperature_K=GAS_TEMPERATURE_K,
            applied_bias_power_W=(
                applied_bias_power_W * chlorine_power_fraction),
        )
        equivalent_temperature = (2.0 / 3.0) * state.mean_electron_energy_eV
        sheath = sheath_transfer.predict(
            positive_ion_flux_m2_s=boundary.positive_ion_flux_m2_s,
            electron_temperature_eV=equivalent_temperature,
            electron_density_m3=state.densities_m3["e"],
            delivered_bias_power_W=(
                boundary.bias_power_delivered_to_ions_W),
        )
        product_condition = _table4_product_condition(
            state=state,
            condition=condition,
            charged_provider=charged_provider,
            wall_boundary=product_wall,
        )
        surface, product_state, surface_residual = (
            _solve_table4_product_surface_fixed_point(
                combined_mechanism,
                boundary,
                sheath.distributions["Cl+"],
                sheath.distributions["Cl2+"],
                product_reactor,
                product_condition,
            )
        )
        product_ion_fluxes = dict(
            product_state.wafer_positive_ion_flux_m2_s)
        predicted_feedback = product_state.chlorine_feedback_lower_bound()
        coupling_residual = _feedback_relative_distance(
            feedback, predicted_feedback)
        feedback = _blend_product_feedback(
            feedback, predicted_feedback, feedback_relaxation)
        state = model.solve(
            condition,
            charged_transport_provider=charged_provider,
            neutral_wall_transport_provider=neutral_provider,
            wall_energy_provider=wall_energy_provider,
            initial_densities_m3=state.densities_m3,
            initial_exhaust_loss_frequency_s_inv=(
                state.exhaust_loss_frequency_s_inv),
            initial_reduced_electric_field_Td=state.reduced_electric_field_Td,
            residual_tolerance=2.0e-7,
            maximum_evaluations=1200,
            maximum_tail_population_fraction=1.0e-6,
            fixed_feedback=feedback,
        )
        final = (state, boundary, sheath, surface, product_state, feedback)
        if iteration >= 4 and coupling_residual <= coupling_tolerance:
            break
    if final is None:  # pragma: no cover
        raise RuntimeError("coupled product solve did not execute")
    if coupling_residual > coupling_tolerance:
        raise RuntimeError(
            "coupled product solve failed fixed-point gate: "
            f"residual={coupling_residual:.6g}, "
            f"iterations={maximum_iterations}"
        )
    return (
        *final,
        iteration,
        coupling_residual,
        product_ion_fluxes,
    )


def _renormalize_transfer_with_center_products(
        transfer, *, center_state, center_product_ion_flux_m2_s,
        target_total_ion_flux_m2_s, target_atomic_chlorine_flux_m2_s):
    """Keep the two published center anchors exact after adding products."""
    product_flux = float(center_product_ion_flux_m2_s)
    total_reference_ion_flux = float(target_total_ion_flux_m2_s)
    chlorine_reference_ion_flux = total_reference_ion_flux - product_flux
    if chlorine_reference_ion_flux <= 0.0:
        raise RuntimeError(
            "predicted center product-ion flux exceeds the published total "
            "ion-current anchor"
        )
    total_reference_chlorine_flux = float(target_atomic_chlorine_flux_m2_s)
    return DiagnosticConditionedChlorineWaferTransfer(
        reference_model_positive_ion_flux_m2_s=(
            center_state.axial_positive_ion_flux_m2_s),
        reference_model_atomic_chlorine_density_m3=(
            center_state.densities_m3["Cl"]),
        reference_gas_temperature_K=transfer.reference_gas_temperature_K,
        reference_wafer_total_ion_flux_m2_s=chlorine_reference_ion_flux,
        reference_wafer_neutral_to_ion_flux_ratio=(
            total_reference_chlorine_flux / chlorine_reference_ion_flux),
        electrode_area_m2=transfer.electrode_area_m2,
        bias_power_to_ion_fraction=transfer.bias_power_to_ion_fraction,
        plasma_potential_eV=transfer.plasma_potential_eV,
        evidence=transfer.evidence,
    )


def run(
        collision_deck: Path, atomic_cl_momentum: Path, *, energy_cells=415,
        bias_power_to_ion_fraction=1.0,
        cl2plus_coverage_mode="saturated_card",
        coupled_product_feedback=False,
        coupled_product_maximum_iterations=16):
    source_rows = _source_rows()
    conditions = {
        (float(row["inductive_power_W"]), float(row["cl2_flow_sccm"]))
        for row in source_rows
    }
    model, replay, grid = _model(
        collision_deck, atomic_cl_momentum, int(energy_cells))
    states = _solve_states(model, conditions)
    charged_provider, neutral_provider, wall_energy_provider = _providers()
    center = states[(400.0, 100.0)]
    reference_wafer_ion_flux = (
        REFERENCE_ION_CURRENT_A_CM2 * 1.0e4 / E_CHARGE_C)
    transfer = DiagnosticConditionedChlorineWaferTransfer(
        reference_model_positive_ion_flux_m2_s=(
            center.axial_positive_ion_flux_m2_s),
        reference_model_atomic_chlorine_density_m3=(
            center.densities_m3["Cl"]),
        reference_gas_temperature_K=GAS_TEMPERATURE_K,
        reference_wafer_total_ion_flux_m2_s=reference_wafer_ion_flux,
        reference_wafer_neutral_to_ion_flux_ratio=(
            REFERENCE_NEUTRAL_TO_ION_RATIO),
        electrode_area_m2=ELECTRODE_AREA_M2,
        bias_power_to_ion_fraction=float(bias_power_to_ion_fraction),
        plasma_potential_eV=PLASMA_POTENTIAL_EV,
        evidence=ChlorineWaferTransferEvidence(
            reference_total_ion_flux=(
                "Mahorowala thesis Chapter 4 estimate: about 2 mA/cm2"),
            reference_neutral_to_ion_ratio=(
                "Mahorowala thesis Chapter 4 estimate: about 100"),
            electrode_area=(
                "Mahorowala thesis Chapter 4 estimate: 400 cm2"),
            bias_power_coupling=(
                "Mahorowala center estimate assigns 80 W over 400 cm2 and "
                "2 mA/cm2 to about 100 eV sheath gain"),
            plasma_potential=(
                "Mahorowala thesis Chapter 4 assumes about 20 V"),
            equipment_transfer=(
                "Lam Alliance geometry/power trend transferred to Lam TCP "
                "9400SE and normalized only on non-etch center estimates"),
            reference_facts_measured=False,
            bias_coupling_measured=False,
            plasma_potential_measured=False,
        ),
    )
    mechanism = ChangSawinClIonSiMechanism()
    cl2plus_mechanism = BaloochCl2IonSiMechanism()
    combined_mechanism = SpeciesResolvedChlorineSiMechanism(
        cl2plus_coverage_mode=str(cl2plus_coverage_mode))
    sheath_transfer = DiagnosticConditionedRFSheathTransfer(
        ion_mass_amu={"Cl+": 35.45, "Cl2+": 70.90},
        electrode_area_m2=ELECTRODE_AREA_M2,
        plasma_potential_eV=PLASMA_POTENTIAL_EV,
        frequency_hz=13.56e6,
        collapse_fraction=1.0,
        phase_count=96,
        steps_per_period=128,
        steps_per_transit=128,
        source=(
            "13.56 MHz Lam Alliance bias frequency transferred as a Lam "
            "equipment-class sensitivity; fully modulated Child sheath"
        ),
    )
    rows = []
    internals = {}
    for source in source_rows:
        source_power = float(source["inductive_power_W"])
        flow = float(source["cl2_flow_sccm"])
        bias = float(source["rf_bias_power_W"])
        state = states[(source_power, flow)]
        boundary = transfer.predict(
            model_positive_ion_flux_m2_s=state.axial_positive_ion_flux_m2_s,
            model_atomic_chlorine_density_m3=state.densities_m3["Cl"],
            gas_temperature_K=GAS_TEMPERATURE_K,
            applied_bias_power_W=bias,
        )
        cl_only = _surface_projection(
            mechanism, boundary, boundary.positive_ion_flux_m2_s["Cl+"])
        cl2plus = _cl2plus_surface_projection(
            cl2plus_mechanism, boundary, cl_only.chlorination_fraction)
        sheath = sheath_transfer.predict(
            positive_ion_flux_m2_s=boundary.positive_ion_flux_m2_s,
            electron_temperature_eV=(
                (2.0 / 3.0) * state.mean_electron_energy_eV),
            electron_density_m3=state.densities_m3["e"],
            delivered_bias_power_W=boundary.bias_power_delivered_to_ions_W,
        )
        cl_iead = sheath.distributions["Cl+"]
        cl2_iead = sheath.distributions["Cl2+"]
        cl_rf = _surface_projection(
            mechanism,
            boundary,
            boundary.positive_ion_flux_m2_s["Cl+"],
            energy_eV=cl_iead.energy_eV,
            weight=cl_iead.weight,
        )
        cl2plus_rf = _cl2plus_surface_projection(
            cl2plus_mechanism,
            boundary,
            cl_rf.chlorination_fraction,
            energy_eV=cl2_iead.energy_eV,
            weight=cl2_iead.weight,
            strict=False,
        )
        combined_rf = _combined_surface_projection(
            combined_mechanism, boundary, cl_iead, cl2_iead)
        all_as_cl = _surface_projection(
            mechanism, boundary, boundary.total_positive_ion_flux_m2_s)
        observed = (
            None if not source["derived_poly_si_removed_nm"]
            else float(source["derived_poly_si_removed_nm"])
        )
        cl_depth = float(cl_only.etch_velocity_m_s * 75.0 * 1.0e9)
        cl2plus_depth = float(cl2plus.etch_velocity_m_s * 75.0 * 1.0e9)
        species_depth = cl_depth + cl2plus_depth
        rf_cl_depth = float(cl_rf.etch_velocity_m_s * 75.0 * 1.0e9)
        rf_cl2plus_depth = float(
            cl2plus_rf.etch_velocity_m_s * 75.0 * 1.0e9)
        rf_species_depth = rf_cl_depth + rf_cl2plus_depth
        combined_rf_depth = float(
            combined_rf.etch_velocity_m_s * 75.0 * 1.0e9)
        all_depth = float(all_as_cl.etch_velocity_m_s * 75.0 * 1.0e9)
        run_id = int(source["run"])
        internals[run_id] = {
            "boundary": boundary,
            "cl_iead": cl_iead,
            "cl2_iead": cl2_iead,
            "state": state,
            "clean_surface": combined_rf,
            "source_power_W": source_power,
            "flow_sccm": flow,
        }
        rows.append({
            "run": run_id,
            "quantitative_status": source["quantitative_status"],
            "inductive_power_W": source_power,
            "rf_bias_power_W": bias,
            "cl2_flow_sccm": flow,
            "observed_feature_depth_nm": observed,
            "reactor_reduced_field_Td": state.reduced_electric_field_Td,
            "reactor_mean_electron_energy_eV": state.mean_electron_energy_eV,
            "reactor_electron_density_m3": state.densities_m3["e"],
            "reactor_atomic_chlorine_density_m3": state.densities_m3["Cl"],
            "reactor_total_axial_ion_flux_m2_s": sum(
                state.axial_positive_ion_flux_m2_s.values()),
            "wafer_total_positive_ion_flux_m2_s": (
                boundary.total_positive_ion_flux_m2_s),
            "wafer_clplus_flux_m2_s": (
                boundary.positive_ion_flux_m2_s["Cl+"]),
            "wafer_cl2plus_flux_m2_s": (
                boundary.positive_ion_flux_m2_s["Cl2+"]),
            "wafer_clplus_fraction": (
                boundary.positive_ion_flux_m2_s["Cl+"]
                / boundary.total_positive_ion_flux_m2_s),
            "wafer_atomic_chlorine_flux_m2_s": (
                boundary.atomic_chlorine_flux_m2_s),
            "wafer_neutral_to_total_ion_ratio": (
                boundary.neutral_to_total_ion_flux_ratio),
            "mean_sheath_energy_gain_eV": (
                boundary.mean_sheath_energy_gain_eV),
            "mean_impact_energy_eV": boundary.mean_impact_energy_eV,
            "ion_power_closure_relative_residual": (
                boundary.power_closure_relative_residual),
            "surface_chlorination_fraction": float(
                cl_only.chlorination_fraction),
            "rf_sheath_surface_chlorination_fraction": float(
                cl_rf.chlorination_fraction),
            "clplus_only_surface_plane_depth_nm_75s": cl_depth,
            "cl2plus_surface_plane_depth_nm_75s": cl2plus_depth,
            "species_resolved_surface_plane_depth_nm_75s": species_depth,
            "rf_sheath_clplus_surface_plane_depth_nm_75s": rf_cl_depth,
            "rf_sheath_cl2plus_surface_plane_depth_nm_75s": rf_cl2plus_depth,
            "rf_sheath_species_resolved_surface_plane_depth_nm_75s": (
                rf_species_depth),
            "joined_species_resolved_surface_plane_depth_nm_75s": (
                combined_rf_depth),
            "all_positive_ions_as_clplus_sensitivity_depth_nm_75s": all_depth,
            "clplus_only_signed_error_percent": (
                None if observed is None else 100.0 * (cl_depth / observed - 1.0)),
            "all_ions_as_clplus_signed_error_percent": (
                None if observed is None else 100.0 * (all_depth / observed - 1.0)),
            "species_resolved_signed_error_percent": (
                None if observed is None
                else 100.0 * (species_depth / observed - 1.0)),
            "rf_sheath_species_resolved_signed_error_percent": (
                None if observed is None
                else 100.0 * (rf_species_depth / observed - 1.0)),
            "joined_species_resolved_signed_error_percent": (
                None if observed is None
                else 100.0 * (combined_rf_depth / observed - 1.0)),
            "rf_sheath_frequency_hz": sheath.frequency_hz,
            "rf_sheath_bias_dc_component_v": sheath.bias_dc_component_v,
            "rf_sheath_dc_v": sheath.sheath_dc_v,
            "rf_sheath_amplitude_v": sheath.sheath_rf_amplitude_v,
            "rf_sheath_clplus_mean_energy_eV": cl_iead.mean_energy_eV,
            "rf_sheath_clplus_standard_deviation_eV": (
                cl_iead.standard_deviation_eV),
            "rf_sheath_cl2plus_mean_energy_eV": cl2_iead.mean_energy_eV,
            "rf_sheath_cl2plus_standard_deviation_eV": (
                cl2_iead.standard_deviation_eV),
            "rf_sheath_clplus_probability_inside_35_100eV": (
                cl_iead.probability_inside(35.0, 100.0)),
            "rf_sheath_cl2plus_probability_inside_26_625eV": (
                cl2_iead.probability_inside(25.998846756576185, 625.0)),
            "rf_sheath_power_closure_relative_residual": (
                sheath.power_closure_relative_residual),
            "surface_energy_inside_measured_35_100eV": (
                35.0 <= boundary.mean_impact_energy_eV <= 100.0),
            "surface_neutral_ion_ratio_inside_measured_domain": (
                boundary.neutral_to_total_ion_flux_ratio >= 5.0),
            "surface_cl2plus_energy_inside_measured_26_625eV": (
                25.998846756576185
                <= boundary.mean_impact_energy_eV <= 625.0),
            "surface_high_chlorination_scope": (
                float(cl_only.chlorination_fraction) >= 0.85),
            "excluded_cl2plus_surface_kinetics": False,
            "feature_transport_applied": False,
            "formal_feature_depth_pass": False,
        })

    center_internal = internals[1]
    center_boundary = center_internal["boundary"]
    reference_product_flux = (
        REFERENCE_SICL2_TO_TOTAL_ION_RATIO
        * center_boundary.total_positive_ion_flux_m2_s
    )
    center_product_surface = _combined_surface_projection(
        combined_mechanism,
        center_boundary,
        center_internal["cl_iead"],
        center_internal["cl2_iead"],
        sicl2_flux_m2_s=reference_product_flux,
    )
    reference_product_gross_rate = float(
        center_product_surface.clplus_removal_rate_si_m2_s
        + center_product_surface.cl2plus_removal_rate_si_m2_s
    )
    geometry = malyshev_1998_lam_geometry(GAP_CM_EQUIPMENT_CLASS_PRIOR)
    product_transfers = {
        "reflective": DiagnosticConditionedEtchProductResidenceTransfer(
            reference_total_ion_flux_m2_s=(
                center_boundary.total_positive_ion_flux_m2_s),
            reference_sicl2_to_total_ion_ratio=(
                REFERENCE_SICL2_TO_TOTAL_ION_RATIO),
            reference_gross_si_source_rate_m2_s=(
                reference_product_gross_rate),
            reference_exhaust_loss_frequency_s_inv=(
                center.exhaust_loss_frequency_s_inv),
            reactor_volume_m3=geometry.neutral_control_volume.value,
            reactor_physical_area_m2=(
                geometry.active_geometry.physical_area_m2),
            gas_temperature_K=GAS_TEMPERATURE_K,
            wall_reactivity=0.0,
        ),
        "reactive": DiagnosticConditionedEtchProductResidenceTransfer(
            reference_total_ion_flux_m2_s=(
                center_boundary.total_positive_ion_flux_m2_s),
            reference_sicl2_to_total_ion_ratio=(
                REFERENCE_SICL2_TO_TOTAL_ION_RATIO),
            reference_gross_si_source_rate_m2_s=(
                reference_product_gross_rate),
            reference_exhaust_loss_frequency_s_inv=(
                center.exhaust_loss_frequency_s_inv),
            reactor_volume_m3=geometry.neutral_control_volume.value,
            reactor_physical_area_m2=(
                geometry.active_geometry.physical_area_m2),
            gas_temperature_K=GAS_TEMPERATURE_K,
            wall_reactivity=1.0,
        ),
    }
    for row in rows:
        internal = internals[row["run"]]
        for wall_name, product_transfer in product_transfers.items():
            surface, residence, fixed_point_residual = (
                _solve_product_surface_fixed_point(
                    combined_mechanism,
                    internal["boundary"],
                    internal["cl_iead"],
                    internal["cl2_iead"],
                    product_transfer,
                    internal["state"].exhaust_loss_frequency_s_inv,
                )
            )
            prefix = f"sicl2_{wall_name}_wall"
            depth = float(surface.etch_velocity_m_s * 75.0 * 1.0e9)
            row[f"{prefix}_flux_m2_s"] = residence.sicl2_flux_m2_s
            row[f"{prefix}_to_total_ion_ratio"] = (
                residence.sicl2_flux_m2_s
                / internal["boundary"].total_positive_ion_flux_m2_s
            )
            row[f"{prefix}_residence_time_s"] = residence.residence_time_s
            row[f"{prefix}_wall_loss_frequency_s_inv"] = (
                residence.wall_loss_frequency_s_inv
            )
            row[f"{prefix}_surface_chlorination_fraction"] = float(
                surface.chlorination_fraction
            )
            row[f"{prefix}_surface_plane_depth_nm_75s"] = depth
            observed = row["observed_feature_depth_nm"]
            row[f"{prefix}_signed_error_percent"] = (
                None if observed is None
                else 100.0 * (depth / observed - 1.0)
            )
            row[f"{prefix}_fixed_point_relative_residual"] = (
                fixed_point_residual
                / max(residence.sicl2_flux_m2_s, 1.0)
            )

    product_reactor = LeeEtchProductLinearReactor()
    product_walls = {
        "reflective": lee_1995_reflective_product_wall(),
        "reactive": lee_1995_reactive_product_wall(),
    }
    for row in rows:
        internal = internals[row["run"]]
        state = internal["state"]
        equivalent_temperature = (2.0 / 3.0) * state.mean_electron_energy_eV
        transport_condition = _condition(
            internal["source_power_W"], internal["flow_sccm"]
        ).transport_condition(equivalent_temperature)
        charged_state = charged_provider.predict(
            transport_condition, state.densities_m3)
        reference_transport = charged_state.positive_ion_transport["Cl+"]
        reference_bohm_speed = bohm_speed(equivalent_temperature, 35.453)
        common_edge = ElectropositiveEdgeFactors(
            axial=(
                reference_transport.axial_flux_velocity_m_s
                / reference_bohm_speed
            ),
            radial=(
                reference_transport.radial_flux_velocity_m_s
                / reference_bohm_speed
            ),
        )
        for wall_name, wall in product_walls.items():
            product_condition = EtchProductPlasmaCondition(
                geometry=transport_condition.geometry,
                neutral_control_volume_m3=(
                    transport_condition.neutral_control_volume.value
                ),
                electron_density_m3=state.densities_m3["e"],
                chlorine_atom_density_m3=state.densities_m3["Cl"],
                chlorine_negative_ion_density_m3=state.densities_m3["Cl-"],
                electron_temperature_eV=equivalent_temperature,
                gas_temperature_K=GAS_TEMPERATURE_K,
                exhaust_loss_frequency_s_inv=(
                    state.exhaust_loss_frequency_s_inv
                ),
                common_edge_factors=common_edge,
                wall_boundary=wall,
                source=(
                    "Mahorowala 200 mm wafer and 60% exposed poly-Si; "
                    "Lee-Graves-Lieberman Table-IV product network and "
                    "common-positive-ion edge-factor closure"
                ),
            )
            surface, product_state, fixed_point_residual = (
                _solve_table4_product_surface_fixed_point(
                    combined_mechanism,
                    internal["boundary"],
                    internal["cl_iead"],
                    internal["cl2_iead"],
                    product_reactor,
                    product_condition,
                )
            )
            prefix = f"table4_product_{wall_name}_wall"
            depth = float(surface.etch_velocity_m_s * 75.0 * 1.0e9)
            product_flux = product_state.wafer_neutral_flux_m2_s["SiCl2"]
            row[f"{prefix}_sicl2_flux_m2_s"] = product_flux
            row[f"{prefix}_sicl2_to_total_ion_ratio"] = (
                product_flux / internal["boundary"].total_positive_ion_flux_m2_s
            )
            row[f"{prefix}_total_product_ion_flux_m2_s"] = sum(
                product_state.wafer_positive_ion_flux_m2_s.values()
            )
            row[f"{prefix}_total_neutral_density_m3"] = (
                product_state.total_neutral_density_m3)
            row[f"{prefix}_total_positive_ion_density_m3"] = (
                product_state.total_positive_ion_density_m3)
            row[f"{prefix}_chlorine_atom_source_m3_s"] = (
                product_state.chlorine_atom_source_m3_s)
            row[f"{prefix}_table4_threshold_power_lower_bound_W_m3"] = (
                product_state.table4_threshold_power_lower_bound_W_m3)
            row[f"{prefix}_elemental_si_density_m3"] = (
                product_state.densities_m3["Si"]
            )
            row[f"{prefix}_surface_chlorination_fraction"] = float(
                surface.chlorination_fraction
            )
            row[f"{prefix}_surface_plane_depth_nm_75s"] = depth
            observed = row["observed_feature_depth_nm"]
            row[f"{prefix}_signed_error_percent"] = (
                None if observed is None
                else 100.0 * (depth / observed - 1.0)
            )
            row[f"{prefix}_fixed_point_relative_residual"] = (
                fixed_point_residual / max(product_flux, 1.0)
            )
            row[f"{prefix}_silicon_inventory_relative_residual"] = (
                product_state.silicon_inventory_relative_residual
            )
            row[f"{prefix}_linear_balance_maximum_relative_residual"] = (
                product_state.linear_balance_maximum_relative_residual
            )
            row[f"{prefix}_matrix_condition_number"] = (
                product_state.matrix_condition_number
            )

    if coupled_product_feedback:
        product_ion_sheath_transfer = DiagnosticConditionedRFSheathTransfer(
            ion_mass_amu=dict(PRODUCT_ION_MASS_AMU),
            electrode_area_m2=ELECTRODE_AREA_M2,
            plasma_potential_eV=PLASMA_POTENTIAL_EV,
            frequency_hz=13.56e6,
            collapse_fraction=1.0,
            phase_count=96,
            steps_per_period=128,
            steps_per_transit=128,
            source=(
                "13.56 MHz Lam equipment-class product-ion sheath "
                "projection at the common solved bias voltage"
            ),
        )
        product_ion_surface = LeeChangProductIonSiSurfaceSensitivity()
        coupled_transfer = transfer
        center_target_total_ion_flux = (
            transfer.reference_wafer_total_ion_flux_m2_s)
        center_target_atomic_chlorine_flux = (
            transfer.reference_wafer_neutral_to_ion_flux_ratio
            * center_target_total_ion_flux)
        center_internal = internals[1]
        center_condition = _condition(
            center_internal["source_power_W"], center_internal["flow_sccm"])
        center_state = center_internal["state"]
        center_feedback = None
        center_product_fluxes = None
        center_transfer_residual = np.inf
        center_result = None
        for center_transfer_iteration in range(1, 9):
            center_result = _solve_coupled_table4_product_state(
                model=model,
                condition=center_condition,
                initial_state=center_state,
                transfer=coupled_transfer,
                sheath_transfer=sheath_transfer,
                combined_mechanism=combined_mechanism,
                product_reactor=product_reactor,
                product_wall=product_walls["reflective"],
                charged_provider=charged_provider,
                neutral_provider=neutral_provider,
                wall_energy_provider=wall_energy_provider,
                applied_bias_power_W=80.0,
                maximum_iterations=coupled_product_maximum_iterations,
                initial_feedback=center_feedback,
                initial_product_ion_fluxes=center_product_fluxes,
            )
            center_state = center_result[0]
            center_feedback = center_result[5]
            center_product_fluxes = center_result[8]
            updated_transfer = _renormalize_transfer_with_center_products(
                coupled_transfer,
                center_state=center_state,
                center_product_ion_flux_m2_s=sum(
                    center_product_fluxes.values()),
                target_total_ion_flux_m2_s=(
                    center_target_total_ion_flux),
                target_atomic_chlorine_flux_m2_s=(
                    center_target_atomic_chlorine_flux),
            )
            center_transfer_residual = max(
                abs(
                    updated_transfer.ion_flux_scale
                    / coupled_transfer.ion_flux_scale - 1.0
                ),
                abs(
                    updated_transfer.radical_flux_scale
                    / coupled_transfer.radical_flux_scale - 1.0
                ),
            )
            coupled_transfer = updated_transfer
            print(
                "coupled center normalization: "
                f"iteration={center_transfer_iteration} "
                f"residual={center_transfer_residual:.6g}",
                flush=True,
            )
            if center_transfer_residual <= 5.0e-4:
                break
        if center_transfer_residual > 5.0e-4:
            raise RuntimeError(
                "coupled center normalization failed fixed-point gate: "
                f"residual={center_transfer_residual:.6g}"
            )
        # One final solve uses the just-updated transfer, so the cached center
        # state and both diagnostic anchors belong to the same fixed point.
        center_result = _solve_coupled_table4_product_state(
            model=model,
            condition=center_condition,
            initial_state=center_state,
            transfer=coupled_transfer,
            sheath_transfer=sheath_transfer,
            combined_mechanism=combined_mechanism,
            product_reactor=product_reactor,
            product_wall=product_walls["reflective"],
            charged_provider=charged_provider,
            neutral_provider=neutral_provider,
            wall_energy_provider=wall_energy_provider,
            applied_bias_power_W=80.0,
            maximum_iterations=coupled_product_maximum_iterations,
            initial_feedback=center_feedback,
            initial_product_ion_fluxes=center_product_fluxes,
        )
        coupled_center_total_ion_flux = (
            center_result[1].total_positive_ion_flux_m2_s
            + sum(center_result[8].values())
        )
        coupled_center_ion_anchor_residual = (
            coupled_center_total_ion_flux
            / center_target_total_ion_flux - 1.0
        )
        coupled_center_chlorine_anchor_residual = (
            center_result[1].atomic_chlorine_flux_m2_s
            / center_target_atomic_chlorine_flux - 1.0
        )
        if max(
            abs(coupled_center_ion_anchor_residual),
            abs(coupled_center_chlorine_anchor_residual),
        ) > 2.0e-3:
            raise RuntimeError(
                "coupled center diagnostic anchors failed closure gate: "
                f"ion={coupled_center_ion_anchor_residual:.6g}, "
                f"chlorine={coupled_center_chlorine_anchor_residual:.6g}"
            )
        for row in rows:
            internal = internals[row["run"]]
            print(
                "coupled Table-IV block: "
                f"run={row['run']} source={internal['source_power_W']:g} "
                f"flow={internal['flow_sccm']:g}",
                flush=True,
            )
            condition = _condition(
                internal["source_power_W"], internal["flow_sccm"])
            result = (
                center_result
                if row["run"] == 1
                else _solve_coupled_table4_product_state(
                    model=model,
                    condition=condition,
                    initial_state=internal["state"],
                    transfer=coupled_transfer,
                    sheath_transfer=sheath_transfer,
                    combined_mechanism=combined_mechanism,
                    product_reactor=product_reactor,
                    product_wall=product_walls["reflective"],
                    charged_provider=charged_provider,
                    neutral_provider=neutral_provider,
                    wall_energy_provider=wall_energy_provider,
                    applied_bias_power_W=float(row["rf_bias_power_W"]),
                    maximum_iterations=coupled_product_maximum_iterations,
                )
            )
            (
                coupled_state,
                coupled_boundary,
                coupled_sheath,
                coupled_surface,
                coupled_product_state,
                coupled_feedback,
                coupled_iterations,
                coupled_residual,
                coupled_product_ion_fluxes,
            ) = result
            prefix = "coupled_table4_product_reflective_wall"
            depth = float(coupled_surface.etch_velocity_m_s * 75.0 * 1.0e9)
            observed = row["observed_feature_depth_nm"]
            chlorine_ion_flux = (
                coupled_boundary.total_positive_ion_flux_m2_s)
            coupled_product_ion_flux = sum(
                coupled_product_ion_fluxes.values())
            product_ion_sheath = (
                product_ion_sheath_transfer.project_from_bias_dc_component(
                    positive_ion_flux_m2_s=coupled_product_ion_fluxes,
                    electron_temperature_eV=(
                        (2.0 / 3.0)
                        * coupled_state.mean_electron_energy_eV),
                    electron_density_m3=coupled_state.densities_m3["e"],
                    bias_dc_component_v=coupled_sheath.bias_dc_component_v,
                )
            )
            product_distributions = dict(product_ion_sheath.distributions)
            product_surface_limits = {
                wall_name: product_ion_surface.evaluate(
                    product_distributions,
                    chlorination_fraction=float(
                        coupled_surface.chlorination_fraction),
                    wall_limit=wall_name,
                    allow_energy_extrapolation=True,
                )
                for wall_name in ("reflective", "reactive")
            }
            row[f"{prefix}_surface_plane_depth_nm_75s"] = depth
            row[f"{prefix}_signed_error_percent"] = (
                None if observed is None
                else 100.0 * (depth / observed - 1.0)
            )
            row[f"{prefix}_iterations"] = coupled_iterations
            row[f"{prefix}_fixed_point_relative_residual"] = coupled_residual
            row[f"{prefix}_chlorine_ion_flux_m2_s"] = chlorine_ion_flux
            row[f"{prefix}_product_ion_flux_m2_s"] = coupled_product_ion_flux
            row[f"{prefix}_total_ion_flux_m2_s"] = (
                chlorine_ion_flux + coupled_product_ion_flux)
            row[f"{prefix}_chlorine_bias_power_fraction"] = (
                chlorine_ion_flux
                / (chlorine_ion_flux + coupled_product_ion_flux)
            )
            row[f"{prefix}_neutral_product_pressure_fraction"] = (
                coupled_product_state.total_neutral_density_m3
                / condition.target_neutral_density_m3
            )
            row[f"{prefix}_product_positive_charge_fraction"] = (
                coupled_product_state.total_positive_ion_density_m3
                / (
                    coupled_state.densities_m3["e"]
                    + coupled_state.densities_m3["Cl-"]
                )
            )
            row[f"{prefix}_minimum_product_power_fraction"] = (
                coupled_feedback.extra_collisional_power_density_W_m3
                / condition.absorbed_power_density_W_m3
            )
            row[f"{prefix}_reactor_electron_density_m3"] = (
                coupled_state.densities_m3["e"])
            row[f"{prefix}_reactor_atomic_chlorine_density_m3"] = (
                coupled_state.densities_m3["Cl"])
            row[f"{prefix}_reactor_mean_electron_energy_eV"] = (
                coupled_state.mean_electron_energy_eV)
            row[f"{prefix}_product_ion_surface_channel_included"] = False
            total_reconstructed_bias_power = (
                coupled_sheath.reconstructed_bias_power_W
                + product_ion_sheath.reconstructed_bias_power_W
            )
            target_bias_power = (
                float(row["rf_bias_power_W"])
                * float(bias_power_to_ion_fraction)
            )
            row[f"{prefix}_all_ion_bias_power_relative_residual"] = (
                (total_reconstructed_bias_power - target_bias_power)
                / max(total_reconstructed_bias_power, target_bias_power, 1.0)
            )
            for species, flux in coupled_product_ion_fluxes.items():
                row[f"{prefix}_{species.lower()}_flux_m2_s"] = flux
                row[f"{prefix}_{species.lower()}_mean_energy_eV"] = (
                    product_ion_sheath.distributions[species].mean_energy_eV)
            base_removal_rate = float(
                coupled_surface.clplus_removal_rate_si_m2_s
                + coupled_surface.cl2plus_removal_rate_si_m2_s)
            for wall_name, product_limit in product_surface_limits.items():
                sensitivity_rate = (
                    base_removal_rate
                    + product_limit.net_removal_rate_si_m2_s)
                sensitivity_depth = float(
                    sensitivity_rate
                    / combined_mechanism.clplus.parameters.bulk_si_atom_density_m3
                    * 75.0 * 1.0e9
                )
                sensitivity_prefix = (
                    f"{prefix}_product_ion_{wall_name}_surface_sensitivity"
                )
                row[f"{sensitivity_prefix}_depth_nm_75s"] = sensitivity_depth
                row[f"{sensitivity_prefix}_signed_error_percent"] = (
                    None if observed is None
                    else 100.0 * (sensitivity_depth / observed - 1.0)
                )
                row[f"{sensitivity_prefix}_gross_rate_si_m2_s"] = (
                    product_limit.total_gross_removal_rate_si_m2_s)
                row[f"{sensitivity_prefix}_deposition_rate_si_m2_s"] = (
                    product_limit.total_deposition_rate_si_m2_s)
    usable = [
        row for row in rows if row["observed_feature_depth_nm"] is not None]
    for key in (
        "clplus_only_signed_error_percent",
        "species_resolved_signed_error_percent",
        "rf_sheath_species_resolved_signed_error_percent",
        "joined_species_resolved_signed_error_percent",
        "sicl2_reflective_wall_signed_error_percent",
        "sicl2_reactive_wall_signed_error_percent",
        "table4_product_reflective_wall_signed_error_percent",
        "table4_product_reactive_wall_signed_error_percent",
        "all_ions_as_clplus_signed_error_percent",
    ):
        values = np.asarray([row[key] for row in usable])
        for row in rows:
            row.setdefault("board_" + key.replace("signed_error_percent", "mape_percent"),
                           float(np.mean(np.abs(values))))
    return {
        "schema": "petch.mahorowala_1998_diagnostic_depth_projection.v4",
        "claim_class": (
            "diagnostic-conditioned reactor-to-surface-plane projection; "
            "not a formal feature-depth prediction"
        ),
        "source": {
            "table_csv": str(SOURCE_CSV.relative_to(ROOT)),
            "table_csv_sha256": _hash(SOURCE_CSV),
            "audit_manifest": str(SOURCE_MANIFEST.relative_to(ROOT)),
            "audit_manifest_sha256": _hash(SOURCE_MANIFEST),
            "cl2plus_surface_csv": str(CL2PLUS_SURFACE_CSV.relative_to(ROOT)),
            "cl2plus_surface_csv_sha256": _hash(CL2PLUS_SURFACE_CSV),
            "cl2plus_surface_manifest": str(
                CL2PLUS_SURFACE_MANIFEST.relative_to(ROOT)),
            "cl2plus_surface_manifest_sha256": _hash(CL2PLUS_SURFACE_MANIFEST),
        },
        "raw_collision_payload_sha256": replay.molecular_replay.raw_payload_sha256,
        "atomic_momentum_payload_sha256": replay.atomic_momentum_payload_sha256,
        "derived_collision_deck_sha256": replay.derived_deck.payload_sha256,
        "raw_collision_bytes_committed": False,
        "energy_grid_cell_count": grid.cell_count,
        "conditioning": {
            "source_power_absorbed_fraction": ABSORBED_SOURCE_POWER_FRACTION,
            "reference_wafer_total_ion_flux_m2_s": reference_wafer_ion_flux,
            "reference_neutral_to_ion_flux_ratio": REFERENCE_NEUTRAL_TO_ION_RATIO,
            "reference_sicl2_to_total_ion_flux_ratio": (
                REFERENCE_SICL2_TO_TOTAL_ION_RATIO),
            "electrode_area_m2": ELECTRODE_AREA_M2,
            "bias_power_to_ion_fraction": float(bias_power_to_ion_fraction),
            "cl2plus_coverage_mode": str(cl2plus_coverage_mode),
            "coupled_product_feedback": bool(coupled_product_feedback),
            "plasma_potential_eV": PLASMA_POTENTIAL_EV,
            "ion_flux_scale": transfer.ion_flux_scale,
            "radical_flux_scale": transfer.radical_flux_scale,
            "etch_rate_or_depth_used_for_conditioning": False,
            "reference_facts_are_source_estimates_not_measurements": True,
            "wafer_diameter_m": WAFER_DIAMETER_M,
            "exposed_polysilicon_area_fraction": (
                EXPOSED_POLYSILICON_AREA_FRACTION
            ),
        },
        "rows": rows,
        "summary": {
            "usable_depth_count": len(usable),
            "clplus_only_surface_plane_mape_percent": float(np.mean(np.abs([
                row["clplus_only_signed_error_percent"] for row in usable
            ]))),
            "species_resolved_surface_plane_mape_percent": float(np.mean(np.abs([
                row["species_resolved_signed_error_percent"] for row in usable
            ]))),
            "rf_sheath_species_resolved_surface_plane_mape_percent": float(
                np.mean(np.abs([
                    row["rf_sheath_species_resolved_signed_error_percent"]
                    for row in usable
                ]))),
            "joined_species_resolved_surface_plane_mape_percent": float(
                np.mean(np.abs([
                    row["joined_species_resolved_signed_error_percent"]
                    for row in usable
                ]))),
            "sicl2_reflective_wall_surface_plane_mape_percent": float(
                np.mean(np.abs([
                    row["sicl2_reflective_wall_signed_error_percent"]
                    for row in usable
                ]))),
            "sicl2_reactive_wall_surface_plane_mape_percent": float(
                np.mean(np.abs([
                    row["sicl2_reactive_wall_signed_error_percent"]
                    for row in usable
                ]))),
            "table4_product_reflective_wall_surface_plane_mape_percent": float(
                np.mean(np.abs([
                    row["table4_product_reflective_wall_signed_error_percent"]
                    for row in usable
                ]))),
            "table4_product_reactive_wall_surface_plane_mape_percent": float(
                np.mean(np.abs([
                    row["table4_product_reactive_wall_signed_error_percent"]
                    for row in usable
                ]))),
            "all_ions_as_clplus_sensitivity_mape_percent": float(np.mean(np.abs([
                row["all_ions_as_clplus_signed_error_percent"] for row in usable
            ]))),
            "measured_energy_domain_count": sum(
                row["surface_energy_inside_measured_35_100eV"] for row in usable),
            "rf_sheath_frequency_evidence": (
                "Lam equipment-class sensitivity; target thesis omits frequency"),
            "rf_sheath_waveform_evidence": (
                "fully modulated collisionless Child-sheath sensitivity; not measured"),
            "etch_product_residence_evidence": (
                "Mahorowala 0.3 source-plane SiCl2/ion estimate; "
                "Lee-Graves-Lieberman reflective/reactive wall limits"),
            "table4_product_network_evidence": (
                "Lee-Graves-Lieberman 1995 Table IV exact particle rates; "
                "Mahorowala reported wafer diameter/exposed-poly fraction; "
                "no depth or etch-rate conditioning"
            ),
            "formal_feature_depth_pass_count": 0,
            **({
                "coupled_table4_product_reflective_wall_surface_plane_mape_percent": float(
                    np.mean(np.abs([
                        row[
                            "coupled_table4_product_reflective_wall_"
                            "signed_error_percent"
                        ]
                        for row in usable
                    ]))
                ),
                "coupled_product_feedback_evidence": (
                    "self-consistent neutral-pressure, positive-charge, "
                    "released-Cl, and Table-IV Arrhenius-threshold power "
                    "lower-bound feedback; product-ion surface response "
                    "remains excluded"
                ),
                "coupled_center_transfer_iterations": (
                    center_transfer_iteration),
                "coupled_center_transfer_fixed_point_relative_residual": (
                    center_transfer_residual),
                "coupled_center_total_ion_anchor_relative_residual": (
                    coupled_center_ion_anchor_residual),
                "coupled_center_atomic_chlorine_anchor_relative_residual": (
                    coupled_center_chlorine_anchor_residual),
                **{
                    (
                        "coupled_table4_product_reflective_wall_product_ion_"
                        f"{wall_name}_surface_sensitivity_mape_percent"
                    ): float(np.mean(np.abs([
                        row[
                            "coupled_table4_product_reflective_wall_"
                            f"product_ion_{wall_name}_surface_sensitivity_"
                            "signed_error_percent"
                        ]
                        for row in usable
                    ])))
                    for wall_name in ("reflective", "reactive")
                },
                "product_ion_surface_sensitivity_evidence": (
                    "Lee Eq. 3 includes all SiClx+ in etching; Chang Cl+ "
                    "yield transferred across product-ion species; not "
                    "supported by a species-resolved product-ion beam"
                ),
            } if coupled_product_feedback else {}),
        },
        "missing_for_formal_depth": [
            "measured per-run or validated bias-power-to-IEAD transfer",
            "measured Cl2+ incidence-angle response for feature sidewalls",
            "deterministic in-feature ion/radical transport and evolving geometry",
            "uncertainties on the center ion current, radical flux, and plasma potential",
            "same-tool reactor geometry and absorbed-power diagnostics",
            "self-consistent SiClx collision-power and chlorine-source feedback",
            "measured chamber-wall Si/Cl coverage state between Lee limits",
            "surface reaction probabilities for SiClx+ product-ion redeposition",
            "measured net Si removal/deposition yield of SiClx+ product ions",
        ],
    }


def _write(result, output: Path):
    output.mkdir(parents=True, exist_ok=True)
    (output / JSON_NAME).write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8")
    rows = result["rows"]
    with (output / CSV_NAME).open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    summary = result["summary"]
    lines = [
        "# Mahorowala 1998 diagnostic-conditioned depth projection",
        "",
        "This is a reactor-to-surface-plane failure-localization run. No etch rate or depth conditioned the reactor, wafer normalization, or surface law. It is not yet a formal feature-depth prediction.",
        "",
        f"- Cl+ only surface-plane MAPE: `{summary['clplus_only_surface_plane_mape_percent']:.2f}%`",
        f"- measured species-resolved Cl+/Cl2+ surface-plane MAPE: `{summary['species_resolved_surface_plane_mape_percent']:.2f}%`",
        f"- power-closed RF-sheath species-resolved surface-plane MAPE: `{summary['rf_sheath_species_resolved_surface_plane_mape_percent']:.2f}%`",
        f"- one-step joined Cl+/Cl2+ surface-plane MAPE: `{summary['joined_species_resolved_surface_plane_mape_percent']:.2f}%`",
        f"- SiCl2 feedback, reflective-wall limit MAPE: `{summary['sicl2_reflective_wall_surface_plane_mape_percent']:.2f}%`",
        f"- SiCl2 feedback, reactive-wall limit MAPE: `{summary['sicl2_reactive_wall_surface_plane_mape_percent']:.2f}%`",
        f"- full Table-IV SiClx network, reflective-wall limit MAPE: `{summary['table4_product_reflective_wall_surface_plane_mape_percent']:.2f}%`",
        f"- full Table-IV SiClx network, reactive-wall limit MAPE: `{summary['table4_product_reactive_wall_surface_plane_mape_percent']:.2f}%`",
        *(
            [
                "- block-coupled Table-IV pressure/charge/Cl/minimum-power "
                "feedback MAPE: "
                f"`{summary['coupled_table4_product_reflective_wall_surface_plane_mape_percent']:.2f}%`"
            ]
            if "coupled_table4_product_reflective_wall_surface_plane_mape_percent"
            in summary else []
        ),
        *(
            [
                "- coupled + Lee all-product-ion reflective sensitivity MAPE: "
                f"`{summary['coupled_table4_product_reflective_wall_product_ion_reflective_surface_sensitivity_mape_percent']:.2f}%`",
                "- coupled + Lee all-product-ion reactive sensitivity MAPE: "
                f"`{summary['coupled_table4_product_reflective_wall_product_ion_reactive_surface_sensitivity_mape_percent']:.2f}%`",
            ]
            if "coupled_table4_product_reflective_wall_product_ion_reflective_surface_sensitivity_mape_percent"
            in summary else []
        ),
        f"- all ions treated as Cl+ sensitivity MAPE: `{summary['all_ions_as_clplus_sensitivity_mape_percent']:.2f}%`",
        f"- points inside the measured 35--100 eV surface domain: `{summary['measured_energy_domain_count']}/{summary['usable_depth_count']}`",
        "- formal held-out feature-depth passes: `0`",
        "",
        "| run | source W | bias W | flow sccm | observed nm | RF species nm | reflective-product nm | reactive-product nm | Cl+ IEAD mean+/-sd eV | Cl2+ IEAD mean+/-sd eV | Cl/ion |",
        "|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        observed = row["observed_feature_depth_nm"]
        lines.append(
            "| {run} | {inductive_power_W:.0f} | {rf_bias_power_W:.0f} | "
            "{cl2_flow_sccm:.0f} | {observed} | "
            "{rf_sheath_species_resolved_surface_plane_depth_nm_75s:.1f} | "
            "{sicl2_reflective_wall_surface_plane_depth_nm_75s:.1f} | "
            "{sicl2_reactive_wall_surface_plane_depth_nm_75s:.1f} | "
            "{rf_sheath_clplus_mean_energy_eV:.1f}+/-{rf_sheath_clplus_standard_deviation_eV:.1f} | "
            "{rf_sheath_cl2plus_mean_energy_eV:.1f}+/-{rf_sheath_cl2plus_standard_deviation_eV:.1f} | "
            "{wafer_neutral_to_total_ion_ratio:.1f} |".format(
                observed="n/a" if observed is None else f"{observed:.1f}",
                **row,
            )
        )
    lines.extend((
        "",
        "## Exact blockers to a formal depth grade",
        "",
        *(f"- {item}" for item in result["missing_for_formal_depth"]),
        "",
    ))
    (output / REPORT_NAME).write_text("\n".join(lines), encoding="utf-8")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--collision-deck", type=Path, required=True)
    parser.add_argument("--atomic-cl-momentum", type=Path, required=True)
    parser.add_argument("--energy-cells", type=int, default=415)
    parser.add_argument(
        "--bias-power-to-ion-fraction", type=float, default=1.0,
        help=(
            "declared bias-power fraction entering the wafer ion-power "
            "ledger; no depth fit"
        ),
    )
    parser.add_argument(
        "--cl2plus-coverage-mode",
        choices=("saturated_card", "source_bounded_linear"),
        default="saturated_card",
        help=(
            "Cl2+ yield treatment: the measured saturated-surface card or "
            "an explicitly nonpredictive interpolation between Chang's "
            "printed bare and saturated slopes"
        ),
    )
    parser.add_argument(
        "--coupled-product-feedback",
        action="store_true",
        help=(
            "iterate the exact Table-IV product block back into reactor "
            "pressure, charge, released chlorine, and minimum collision power"
        ),
    )
    parser.add_argument(
        "--coupled-product-maximum-iterations", type=int, default=16,
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    arguments = parser.parse_args()
    result = run(
        arguments.collision_deck,
        arguments.atomic_cl_momentum,
        energy_cells=arguments.energy_cells,
        bias_power_to_ion_fraction=arguments.bias_power_to_ion_fraction,
        cl2plus_coverage_mode=arguments.cl2plus_coverage_mode,
        coupled_product_feedback=arguments.coupled_product_feedback,
        coupled_product_maximum_iterations=(
            arguments.coupled_product_maximum_iterations),
    )
    _write(result, arguments.output)
    print(json.dumps(result["summary"], sort_keys=True))


if __name__ == "__main__":
    main()
