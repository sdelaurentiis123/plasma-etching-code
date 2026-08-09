import math

import pytest

from petch.reactor_global.chlorine_etch_product_reactor import (
    EtchProductPlasmaCondition,
    LeeEtchProductLinearReactor,
    lee_1995_reactive_product_wall,
    lee_1995_reflective_product_wall,
)
from petch.reactor_global.geometry import (
    CylindricalReactor,
    ElectropositiveEdgeFactors,
)


def _condition(wall):
    geometry = CylindricalReactor(radius_m=0.215, length_m=0.065)
    return EtchProductPlasmaCondition(
        geometry=geometry,
        neutral_control_volume_m3=0.043,
        electron_density_m3=2.0e16,
        chlorine_atom_density_m3=8.0e19,
        chlorine_negative_ion_density_m3=1.0e16,
        electron_temperature_eV=3.6,
        gas_temperature_K=500.0,
        exhaust_loss_frequency_s_inv=4.0,
        common_edge_factors=ElectropositiveEdgeFactors(
            axial=0.12,
            radial=0.08,
        ),
        wall_boundary=wall,
        source="manufactured fixed base plasma",
    )


@pytest.mark.parametrize("wall", [
    lee_1995_reflective_product_wall(),
    lee_1995_reactive_product_wall(),
])
def test_linear_product_solve_is_positive_and_silicon_closed(wall):
    solution = LeeEtchProductLinearReactor().solve(
        _condition(wall),
        gross_si_removal_flux_m2_s=1.0e19,
        exposed_silicon_area_m2=math.pi * 0.1 ** 2 * 0.6,
    )

    assert all(value >= 0.0 for value in solution.densities_m3.values())
    assert abs(solution.silicon_inventory_relative_residual) < 2.0e-12
    assert solution.linear_balance_maximum_relative_residual < 2.0e-12
    assert solution.wafer_neutral_flux_m2_s["SiCl2"] > 0.0
    assert solution.wafer_positive_ion_flux_m2_s["Si+"] > 0.0
    assert solution.chlorine_atom_source_m3_s > 0.0
    assert solution.table4_threshold_power_lower_bound_W_m3 > 0.0
    assert sum(solution.positive_ion_wall_loss_molecule_s.values()) > 0.0
    feedback = solution.chlorine_feedback_lower_bound()
    assert feedback.extra_neutral_density_m3 == pytest.approx(
        solution.total_neutral_density_m3)
    assert feedback.extra_positive_charge_density_m3 == pytest.approx(
        solution.total_positive_ion_density_m3)
    assert not feedback.supports_prediction


def test_reactive_wall_reduces_elemental_silicon_residence():
    model = LeeEtchProductLinearReactor()
    reflective = model.solve(
        _condition(lee_1995_reflective_product_wall()),
        gross_si_removal_flux_m2_s=1.0e19,
        exposed_silicon_area_m2=0.02,
    )
    reactive = model.solve(
        _condition(lee_1995_reactive_product_wall()),
        gross_si_removal_flux_m2_s=1.0e19,
        exposed_silicon_area_m2=0.02,
    )

    assert reactive.densities_m3["Si"] < reflective.densities_m3["Si"]
    assert reactive.wall_deposition_loss_molecule_s["Si"] > 0.0
