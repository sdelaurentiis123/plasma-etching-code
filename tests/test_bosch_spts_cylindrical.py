from pathlib import Path

import numpy as np
import pytest

from petch.bosch_process_data import load_bosch_process_traces
from petch.reactor_global.bosch_spts_cylindrical import (
    BoschSPTSCylindricalParameters,
    DeterministicBoschSPTSCylindricalReactorToWafer,
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
