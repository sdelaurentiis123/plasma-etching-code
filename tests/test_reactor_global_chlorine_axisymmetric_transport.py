from dataclasses import replace

import numpy as np
import pytest

from petch.reactor_global import (
    AxisymmetricFiniteVolumeGrid,
    ChlorineAxisymmetricTransportInput,
    CylindricalReactor,
    DeterministicChlorineAxisymmetricTransport,
    normalized_annular_skin_source,
)


def _provider():
    geometry = CylindricalReactor(radius_m=0.15, length_m=0.08)
    grid = AxisymmetricFiniteVolumeGrid.uniform(
        geometry, radial_cell_count=6, axial_cell_count=5)
    provider = DeterministicChlorineAxisymmetricTransport(
        grid=grid,
        ion_temperature_eV=0.12,
        positive_negative_recombination_m3_s=5.0e-14,
    )
    return provider, grid


def _state(geometry, density=None):
    return ChlorineAxisymmetricTransportInput(
        condition_id="manufactured-cl2-icp",
        geometry=geometry,
        charged_density_m3=(
            {"Cl+": 3.0e16, "Cl2+": 8.0e16, "Cl-": 7.0e16}
            if density is None else density
        ),
        total_neutral_density_m3=3.2e20,
        mean_electron_energy_eV=4.5,
        source="manufactured global chlorine state",
    )


def test_annular_source_is_normalized_and_peaks_away_from_axis():
    _, grid = _provider()
    source = normalized_annular_skin_source(
        grid,
        axial_skin_depth_m=0.02,
        ring_radius_m=0.10,
        radial_width_m=0.025,
    )
    average = np.sum(source * grid.cell_volume_m3) / grid.geometry.volume_m3
    assert average == pytest.approx(1.0, abs=2.0e-15)
    radial_profile = source[:, -1]
    assert np.argmax(radial_profile) > 0
    assert radial_profile[-1] < np.max(radial_profile)


def test_chlorine_axisymmetric_provider_closes_particle_ledgers_and_jvp():
    provider, grid = _provider()
    state = _state(grid.geometry)
    source = normalized_annular_skin_source(
        grid,
        axial_skin_depth_m=0.025,
        ring_radius_m=0.09,
        radial_width_m=0.035,
    )
    result = provider.predict(
        state,
        source,
        source_moment_name="declared_annular",
        source_moment_provenance="manufactured fixed annular source",
        source_moment_measured_or_validated=False,
        wafer_radius_m=0.10,
        relative_tolerance=2.0e-9,
        maximum_iterations=800,
    )
    assert result.total_wafer_positive_ion_flux_m2_s > 0.0
    assert result.lift_result.solution.maximum_species_ledger_relative_residual < 1.0e-12
    assert result.supports_absolute_wafer_flux_prediction is False
    assert result.supports_implicit_differentiation is True

    target = np.asarray([
        state.charged_density_m3[name]
        for name in ("Cl+", "Cl2+", "Cl-")
    ])
    direction = target * np.array([0.08, -0.05, 0.04])
    tangent, flux_tangent = provider.target_inventory_jvp(result, direction)
    step = 2.0e-4

    def perturbed(sign):
        density = dict(zip(
            ("Cl+", "Cl2+", "Cl-"), target + sign * step * direction))
        return provider.predict(
            replace(state, charged_density_m3=density),
            source,
            source_moment_name="declared_annular",
            source_moment_provenance="manufactured fixed annular source",
            source_moment_measured_or_validated=False,
            wafer_radius_m=0.10,
            initial_electrostatic_potential_V=(
                result.lift_result.electrostatic_potential_V),
            relative_tolerance=2.0e-9,
            maximum_iterations=800,
        )

    plus = perturbed(1.0)
    minus = perturbed(-1.0)
    for name in ("Cl+", "Cl2+"):
        finite = (
            plus.wafer_positive_ion_flux_m2_s[name]
            - minus.wafer_positive_ion_flux_m2_s[name]
        ) / (2.0 * step)
        assert flux_tangent[name] == pytest.approx(
            finite, rel=4.0e-5, abs=2.0e8)
    assert tangent.maximum_linearized_fixed_point_relative_residual < 1.0e-8
