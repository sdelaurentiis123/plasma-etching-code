from pathlib import Path
from dataclasses import replace

import numpy as np
import pytest

from petch.bosch_process_data import load_bosch_process_traces
from petch.reactor_global.bosch_spts_cylindrical import (
    BoschSPTSCylindricalParameters,
    BoschSPTSWaferIonTransmissionLaw,
    DeterministicBoschSPTSCylindricalReactorToWafer,
    bosch_real_zernike_design,
    bosch_real_zernike_modes,
)
from petch.reactor_global.bosch_spts_reduced import (
    BoschSPTSReducedParameters, DeterministicBoschSPTSReactorToWafer,
)


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "experimental" / "zenodo_17122442"


@pytest.fixture(scope="module")
def trace():
    return load_bosch_process_traces(
        DATA / "Process_data.nc", DATA / "Dictionary_process.nc")[0]


def test_axisymmetric_cylindrical_tier_matches_axisymmetric_area_average(trace):
    reduced = BoschSPTSReducedParameters(
        radial_cell_count=12, axial_cell_count=10)
    axisymmetric = DeterministicBoschSPTSReactorToWafer(reduced).solve(trace)
    cylindrical = DeterministicBoschSPTSCylindricalReactorToWafer(
        BoschSPTSCylindricalParameters(
            reduced=reduced, azimuthal_cell_count=12)).solve(
                trace, x_m=np.array([0.0, 0.03, -0.04, 0.06]),
                y_m=np.array([0.0, 0.04, 0.02, -0.03]))

    assert cylindrical.maximum_cylindrical_species_ledger_relative_residual < 1e-12
    assert cylindrical.maximum_cylindrical_linear_system_relative_residual < 1e-12
    assert cylindrical.wafer_area_average_flux_m2_s == pytest.approx(
        axisymmetric.wafer_area_average_flux_m2_s, rel=3.0e-12, abs=1.0)


def test_source_harmonics_produce_positive_azimuthal_wafer_response(trace):
    reduced = BoschSPTSReducedParameters(
        radial_cell_count=10, axial_cell_count=8,
        source_ring_radius_m=0.12, source_radial_width_m=0.025)
    model = DeterministicBoschSPTSCylindricalReactorToWafer(
        BoschSPTSCylindricalParameters(
            reduced=reduced, azimuthal_cell_count=16,
            source_cosine_coefficients=((0.3, -0.1), (), (0.2, 0.15)),
            source_sine_coefficients=((-0.2, 0.1), (), (0.1, -0.1))))
    radius = 0.08
    phi = np.linspace(0.0, 2.0 * np.pi, 8, endpoint=False)
    solution = model.solve(
        trace, x_m=radius * np.cos(phi), y_m=radius * np.sin(phi))

    assert np.all(solution.point_flux_m2_s >= 0.0)
    assert np.any(np.ptp(solution.point_flux_m2_s[:, 0], axis=1) > 0.0)
    assert np.any(np.ptp(solution.point_flux_m2_s[:, 2], axis=1) > 0.0)
    assert solution.provenance["target_depth_used"] is False


def test_cylindrical_point_flux_jvp_is_exact(trace):
    reduced = BoschSPTSReducedParameters(radial_cell_count=8, axial_cell_count=7)
    model = DeterministicBoschSPTSCylindricalReactorToWafer(
        BoschSPTSCylindricalParameters(
            reduced=reduced, azimuthal_cell_count=8,
            source_cosine_coefficients=((0.1,), (0.05,), (-0.1,)),
            source_sine_coefficients=((0.0,), (-0.05,), (0.08,))))
    x = np.array([0.0, 0.02, -0.04, 0.06])
    y = np.array([0.0, -0.03, 0.05, 0.01])
    solution = model.solve(trace, x_m=x, y_m=y)
    tangent = 0.03 * solution.reactor.volume_average_density_m3
    epsilon = 1.0e-5
    unit = solution.point_flux_m2_s / np.maximum(
        solution.reactor.volume_average_density_m3[:, :, None], 1.0)
    finite_difference = (
        (solution.reactor.volume_average_density_m3 + epsilon * tangent)[:, :, None]
        * unit
        - (solution.reactor.volume_average_density_m3 - epsilon * tangent)[:, :, None]
        * unit) / (2.0 * epsilon)

    assert model.density_to_point_flux_jvp(
        tangent, x_m=x, y_m=y) == pytest.approx(
            finite_difference, rel=2.0e-9, abs=1.0e5)


def test_source_harmonics_reuse_factorization_and_match_fresh_model(trace):
    reduced = BoschSPTSReducedParameters(radial_cell_count=8, axial_cell_count=7)
    cosine = ((0.2, -0.05), (-0.1, 0.08), (0.15, 0.04))
    sine = ((-0.1, 0.07), (0.05, -0.04), (0.12, -0.06))
    x = np.array([0.01, 0.03, -0.04, 0.06, -0.07])
    y = np.array([0.02, -0.04, 0.05, 0.01, -0.02])
    reusable = DeterministicBoschSPTSCylindricalReactorToWafer(
        BoschSPTSCylindricalParameters(
            reduced=reduced, azimuthal_cell_count=8))
    response = reusable.source_response(
        x_m=x, y_m=y, source_cosine_coefficients=cosine,
        source_sine_coefficients=sine)
    reused = reusable.solve(
        trace, x_m=x, y_m=y, source_response=response)
    fresh = DeterministicBoschSPTSCylindricalReactorToWafer(
        BoschSPTSCylindricalParameters(
            reduced=reduced, azimuthal_cell_count=8,
            source_cosine_coefficients=cosine,
            source_sine_coefficients=sine)).solve(trace, x_m=x, y_m=y)

    assert response.provenance["transport_factorization_reused"] is True
    assert reused.point_flux_m2_s == pytest.approx(
        fresh.point_flux_m2_s, rel=2.0e-14, abs=1.0)
    assert reused.wafer_area_average_flux_m2_s == pytest.approx(
        fresh.wafer_area_average_flux_m2_s, rel=2.0e-14, abs=1.0)


def test_radial_source_geometry_reuses_unchanged_transport_operator(trace):
    reduced = BoschSPTSReducedParameters(radial_cell_count=8, axial_cell_count=7)
    x = np.array([0.01, 0.03, -0.04, 0.06, -0.07])
    y = np.array([0.02, -0.04, 0.05, 0.01, -0.02])
    reusable = DeterministicBoschSPTSCylindricalReactorToWafer(
        BoschSPTSCylindricalParameters(
            reduced=reduced, azimuthal_cell_count=8))
    response = reusable.source_response(
        x_m=x, y_m=y, source_ring_radius_m=0.12,
        source_radial_width_m=0.025,
        source_central_fraction=(0.2, 0.7, 0.1))
    reused = reusable.solve(
        trace, x_m=x, y_m=y, source_response=response)
    changed_reduced = replace(
        reduced, source_ring_radius_m=0.12, source_radial_width_m=0.025,
        source_central_fraction=(0.2, 0.7, 0.1))
    fresh = DeterministicBoschSPTSCylindricalReactorToWafer(
        BoschSPTSCylindricalParameters(
            reduced=changed_reduced, azimuthal_cell_count=8)).solve(
                trace, x_m=x, y_m=y)

    assert reused.point_flux_m2_s == pytest.approx(
        fresh.point_flux_m2_s, rel=2.0e-14, abs=1.0)
    assert response.provenance["transport_factorization_reused"] is True


def test_species_radial_source_moments_match_fresh_model(trace):
    reduced = BoschSPTSReducedParameters(radial_cell_count=8, axial_cell_count=7)
    rings = (0.06, 0.11, 0.14)
    widths = (0.05, 0.025, 0.015)
    x = np.array([0.01, 0.03, -0.04, 0.06, -0.07])
    y = np.array([0.02, -0.04, 0.05, 0.01, -0.02])
    reusable = DeterministicBoschSPTSCylindricalReactorToWafer(
        BoschSPTSCylindricalParameters(
            reduced=reduced, azimuthal_cell_count=8))
    response = reusable.source_response(
        x_m=x, y_m=y, species_source_ring_radius_m=rings,
        species_source_radial_width_m=widths)
    reused = reusable.solve(
        trace, x_m=x, y_m=y, source_response=response)
    fresh = DeterministicBoschSPTSCylindricalReactorToWafer(
        BoschSPTSCylindricalParameters(
            reduced=reduced, azimuthal_cell_count=8,
            species_source_ring_radius_m=rings,
            species_source_radial_width_m=widths)).solve(
                trace, x_m=x, y_m=y)

    assert reused.point_flux_m2_s == pytest.approx(
        fresh.point_flux_m2_s, rel=2.0e-14, abs=1.0)
    assert response.provenance["parameters"]["species_source_ring_radius_m"] == list(rings)
    assert response.provenance["parameters"]["species_source_radial_width_m"] == list(widths)


def test_edge_sheath_focus_is_positive_and_conserves_total_ion_current(trace):
    reduced = BoschSPTSReducedParameters(radial_cell_count=16, axial_cell_count=10)
    x = np.array([0.0, 0.04, 0.08, 0.095, -0.095])
    y = np.zeros_like(x)
    model = DeterministicBoschSPTSCylindricalReactorToWafer(
        BoschSPTSCylindricalParameters(
            reduced=reduced, azimuthal_cell_count=12))
    base = model.source_response(x_m=x, y_m=y)
    focused = model.source_response(
        x_m=x, y_m=y, ion_edge_focus_amplitude=0.3,
        ion_edge_focus_onset_radius_m=0.09,
        ion_edge_focus_width_m=0.002)

    assert focused.unit_wafer_average_flux_per_density_m_s == pytest.approx(
        base.unit_wafer_average_flux_per_density_m_s, rel=2.0e-14, abs=1e-16)
    assert focused.unit_point_flux_per_density_m_s[:2] == pytest.approx(
        base.unit_point_flux_per_density_m_s[:2], rel=2.0e-14, abs=1e-16)
    base_ratio = (
        base.unit_point_flux_per_density_m_s[2, 3]
        / base.unit_point_flux_per_density_m_s[2, 0])
    focused_ratio = (
        focused.unit_point_flux_per_density_m_s[2, 3]
        / focused.unit_point_flux_per_density_m_s[2, 0])
    assert focused_ratio > base_ratio
    assert np.all(focused.unit_point_flux_per_density_m_s >= 0.0)
    assert focused.provenance["total_wafer_ion_current_conserved"] is True


def test_wall_conditioning_changes_only_neutral_nonwafer_transfer(trace):
    base_reduced = BoschSPTSReducedParameters(
        radial_cell_count=10, axial_cell_count=8)
    conditioned_reduced = replace(
        base_reduced, neutral_wall_loss_multiplier=2.0)
    x = np.array([0.0, 0.04, 0.08, -0.08])
    y = np.zeros_like(x)
    base = DeterministicBoschSPTSCylindricalReactorToWafer(
        BoschSPTSCylindricalParameters(
            reduced=base_reduced, azimuthal_cell_count=8)).source_response(
                x_m=x, y_m=y)
    conditioned = DeterministicBoschSPTSCylindricalReactorToWafer(
        BoschSPTSCylindricalParameters(
            reduced=conditioned_reduced, azimuthal_cell_count=8)).source_response(
                x_m=x, y_m=y)

    assert not np.array_equal(
        conditioned.unit_point_flux_per_density_m_s[:2],
        base.unit_point_flux_per_density_m_s[:2],
    )
    assert conditioned.unit_point_flux_per_density_m_s[2] == pytest.approx(
        base.unit_point_flux_per_density_m_s[2], rel=2e-14, abs=1e-16)


def test_real_zernike_basis_has_frozen_complete_ordering():
    assert [len(bosch_real_zernike_modes(order)) for order in range(1, 11)] == [
        2, 5, 9, 14, 20, 27, 35, 44, 54, 65,
    ]
    design = bosch_real_zernike_design(
        2, np.array([0.0]), np.array([1.234]))

    assert bosch_real_zernike_modes(2) == (
        (1, 1, "cos"),
        (1, 1, "sin"),
        (2, 0, "cos"),
        (2, 2, "cos"),
        (2, 2, "sin"),
    )
    assert design[0] == pytest.approx([0.0, 0.0, -1.0, 0.0, 0.0])


def test_wafer_ion_transmission_is_positive_and_conserves_current(trace):
    reduced = BoschSPTSReducedParameters(
        radial_cell_count=18, axial_cell_count=10)
    x = np.array([0.0, 0.04, -0.06, 0.08, -0.095])
    y = np.array([0.0, 0.03, 0.02, -0.04, 0.0])
    model = DeterministicBoschSPTSCylindricalReactorToWafer(
        BoschSPTSCylindricalParameters(
            reduced=reduced, azimuthal_cell_count=16,
            ion_edge_focus_amplitude=0.25,
            ion_edge_focus_onset_radius_m=0.09,
            ion_edge_focus_width_m=0.003))
    base = model.source_response(x_m=x, y_m=y)
    law = BoschSPTSWaferIonTransmissionLaw(
        static_maximum_order=2,
        static_coefficients=(0.015, -0.01, 0.025, 0.008, -0.006),
        dynamic_maximum_order=2,
        dynamic_coefficients=(0.001, -0.002, 0.004, 0.001, -0.001),
    )
    mapped = model.source_response(
        x_m=x,
        y_m=y,
        ion_transmission_law=law,
        c4f8_platen_vpp_rms_V=637.4409584828442,
    )

    assert mapped.unit_wafer_average_flux_per_density_m_s == pytest.approx(
        base.unit_wafer_average_flux_per_density_m_s,
        rel=2.0e-14,
        abs=1.0e-16,
    )
    assert mapped.unit_point_flux_per_density_m_s[:2] == pytest.approx(
        base.unit_point_flux_per_density_m_s[:2],
        rel=2.0e-14,
        abs=1.0e-16,
    )
    assert np.all(mapped.unit_point_flux_per_density_m_s > 0.0)
    assert not np.allclose(
        mapped.unit_point_flux_per_density_m_s[2],
        base.unit_point_flux_per_density_m_s[2],
    )
    provenance = mapped.provenance["wafer_ion_transmission"]
    assert provenance["enabled"] is True
    assert provenance["relative_total_ion_current_residual"] <= 2.0e-14
    assert mapped.provenance["total_wafer_ion_current_conserved"] is True


def test_wafer_ion_transmission_voltage_is_smooth_and_domain_bounded(trace):
    reduced = BoschSPTSReducedParameters(
        radial_cell_count=14, axial_cell_count=9)
    x = np.array([0.0, 0.04, 0.08, -0.095])
    y = np.zeros_like(x)
    model = DeterministicBoschSPTSCylindricalReactorToWafer(
        BoschSPTSCylindricalParameters(
            reduced=reduced, azimuthal_cell_count=12))
    law = BoschSPTSWaferIonTransmissionLaw(
        static_maximum_order=1,
        static_coefficients=(0.0, 0.0),
        dynamic_maximum_order=2,
        dynamic_coefficients=(0.0, 0.0, 0.01, 0.0, 0.0),
    )
    low = model.source_response(
        x_m=x, y_m=y, ion_transmission_law=law,
        c4f8_platen_vpp_rms_V=626.9533265149638)
    high = model.source_response(
        x_m=x, y_m=y, ion_transmission_law=law,
        c4f8_platen_vpp_rms_V=643.534317555529)
    low_edge_center = (
        low.unit_point_flux_per_density_m_s[2, -1]
        / low.unit_point_flux_per_density_m_s[2, 0])
    high_edge_center = (
        high.unit_point_flux_per_density_m_s[2, -1]
        / high.unit_point_flux_per_density_m_s[2, 0])

    assert high_edge_center > low_edge_center
    with pytest.raises(ValueError, match="outside the frozen domain"):
        model.source_response(
            x_m=x, y_m=y, ion_transmission_law=law,
            c4f8_platen_vpp_rms_V=650.0)
    with pytest.raises(ValueError, match="must be paired"):
        model.source_response(x_m=x, y_m=y, ion_transmission_law=law)


def test_wafer_ion_transmission_rejects_excessive_log_field():
    with pytest.raises(ValueError, match="exceeds its frozen bound"):
        BoschSPTSWaferIonTransmissionLaw(
            static_maximum_order=10,
            static_coefficients=(0.1,) * 65,
        )
