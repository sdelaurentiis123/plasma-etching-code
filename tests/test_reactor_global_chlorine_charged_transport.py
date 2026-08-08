import numpy as np
import pytest

from petch.reactor_global import (
    LYMBEROPOULOS_1995_REDUCED_CL2PLUS_MOBILITY_M_INV_V_INV_S_INV,
    LYMBEROPOULOS_1995_REDUCED_CLMINUS_MOBILITY_M_INV_V_INV_S_INV,
    LYMBEROPOULOS_1995_REDUCED_CLPLUS_MOBILITY_M_INV_V_INV_S_INV,
    ChlorineFixedPressureCondition,
    CylindricalReactor,
    LeeEconomouChlorineChargedTransportProvider,
    ReactorScalarInput,
    lymberopoulos_economou_1995_chlorine_reduced_ion_mobilities,
    ramamurthi_economou_2002_chlorine_reduced_ion_mobilities,
    standard_volume_flow_molecules_s,
)
from petch.reactor_global.network import E_CHARGE_C
from petch.reactor_global.transport import ATOMIC_MASS_UNIT_KG


def _scalar(value, unit, *, evidence_kind="published_model"):
    return ReactorScalarInput(
        value=value,
        unit=unit,
        source="independent charged-transport regression input",
        evidence_kind=evidence_kind,
        relative_uncertainty=None,
    )


def _condition():
    geometry = CylindricalReactor(radius_m=0.1525, length_m=0.075)
    return ChlorineFixedPressureCondition(
        condition_id="charged-transport-source-reproduction",
        geometry=geometry,
        neutral_control_volume=_scalar(geometry.volume_m3, "m3"),
        pressure=_scalar(1.333223684, "Pa"),
        gas_temperature=_scalar(500.0, "K"),
        electron_temperature=_scalar(3.0, "eV"),
        chlorine_molecule_feed=_scalar(
            standard_volume_flow_molecules_s(
                35.0,
                standard_temperature_K=273.15,
                standard_pressure_Pa=101325.0,
            ),
            "molecule s^-1",
        ),
        source_power=_scalar(1000.0, "W"),
    )


def _densities(*, electronegativity):
    electron_density = 1.0e16
    negative_density = electronegativity * electron_density
    positive_total = electron_density + negative_density
    return {
        "e": electron_density,
        "Cl2": 1.5e20,
        "Cl": 4.0e19,
        "Cl2+": 0.6 * positive_total,
        "Cl+": 0.4 * positive_total,
        "Cl-": negative_density,
    }


def _positive_mobilities():
    mobilities = (
        lymberopoulos_economou_1995_chlorine_reduced_ion_mobilities())
    return {species: mobilities[species] for species in ("Cl2+", "Cl+")}


def test_published_reduced_mobilities_preserve_values_units_and_conflict():
    source_1995 = (
        lymberopoulos_economou_1995_chlorine_reduced_ion_mobilities())
    source_2002 = ramamurthi_economou_2002_chlorine_reduced_ion_mobilities()
    expected = {
        "Cl2+": (
            LYMBEROPOULOS_1995_REDUCED_CL2PLUS_MOBILITY_M_INV_V_INV_S_INV),
        "Cl+": (
            LYMBEROPOULOS_1995_REDUCED_CLPLUS_MOBILITY_M_INV_V_INV_S_INV),
        "Cl-": (
            LYMBEROPOULOS_1995_REDUCED_CLMINUS_MOBILITY_M_INV_V_INV_S_INV),
    }
    expected_300_K_eV = 300.0 * 1.380649e-23 / E_CHARGE_C

    assert set(source_1995) == set(expected)
    assert set(source_2002) == set(expected)
    for species, value in expected.items():
        assert source_1995[species].reduced_mobility_m_inv_V_inv_s_inv == value
        assert source_2002[species].reduced_mobility_m_inv_V_inv_s_inv == value
        assert source_1995[species].provenance[
            "source_value_cm-1_V-1_s-1"] == value / 100.0
        assert source_1995[species].reference_ion_temperature_eV == 0.12
        assert source_2002[species].reference_ion_temperature_eV == pytest.approx(
            expected_300_K_eV, rel=1.0e-15)
        assert not source_1995[species].supports_prediction
        assert not source_2002[species].supports_prediction

    with pytest.raises(ValueError, match="ion temperature"):
        source_1995["Cl+"].evaluate(
            total_neutral_density_m3=2.0e20,
            ion_temperature_eV=expected_300_K_eV,
        )


def test_reduced_mobility_evaluates_exact_density_scaling():
    model = (
        lymberopoulos_economou_1995_chlorine_reduced_ion_mobilities()[
            "Cl2+"])
    density = 2.5e20
    state = model.evaluate(
        total_neutral_density_m3=density,
        ion_temperature_eV=0.12,
    )
    assert state.mobility_m2_V_s == pytest.approx(
        LYMBEROPOULOS_1995_REDUCED_CL2PLUS_MOBILITY_M_INV_V_INV_S_INV
        / density,
        rel=1.0e-15,
    )
    assert state.total_neutral_density_m3 == density
    assert not state.supports_prediction


def test_provider_exposes_every_mobility_and_edge_intermediate():
    condition = _condition()
    provider = LeeEconomouChlorineChargedTransportProvider(
        reduced_mobilities=_positive_mobilities(),
        ion_temperature=_scalar(0.12, "eV"),
    )
    state = provider.predict(condition, _densities(electronegativity=1.0))

    assert state.geometry == condition.geometry
    assert not state.supports_prediction
    for species, mass_amu in {"Cl2+": 70.906, "Cl+": 35.453}.items():
        transport = state.positive_ion_transport[species]
        provenance = transport.provenance
        mass_kg = mass_amu * ATOMIC_MASS_UNIT_KG
        mobility = provenance["mobility_m2_V_s"]
        expected_collision_frequency = E_CHARGE_C / (mass_kg * mobility)
        expected_mean_speed = np.sqrt(
            8.0 * E_CHARGE_C * 0.12 / (np.pi * mass_kg))
        assert provenance[
            "momentum_collision_frequency_s_inv"] == pytest.approx(
                expected_collision_frequency, rel=1.0e-14)
        assert provenance["mean_ion_speed_m_s"] == pytest.approx(
            expected_mean_speed, rel=1.0e-14)
        assert provenance["momentum_mean_free_path_m"] == pytest.approx(
            expected_mean_speed / expected_collision_frequency,
            rel=1.0e-14,
        )
        assert provenance["ambipolar_diffusivity_m2_s"] == pytest.approx(
            mobility * condition.electron_temperature.value,
            rel=1.0e-14,
        )
        assert provenance["total_neutral_density_m3"] == pytest.approx(1.9e20)
        assert provenance["electronegativity"] == 1.0
        assert transport.axial_flux_velocity_m_s > 0.0
        assert transport.radial_flux_velocity_m_s > 0.0


def test_electronegative_edge_transport_is_recomputed_from_density_state():
    condition = _condition()
    provider = LeeEconomouChlorineChargedTransportProvider(
        reduced_mobilities=_positive_mobilities(),
        ion_temperature=_scalar(0.12, "eV"),
    )
    low_alpha = provider.predict(
        condition, _densities(electronegativity=0.1))
    high_alpha = provider.predict(
        condition, _densities(electronegativity=10.0))

    for species in ("Cl2+", "Cl+"):
        low = low_alpha.positive_ion_transport[species]
        high = high_alpha.positive_ion_transport[species]
        assert high.axial_flux_velocity_m_s < low.axial_flux_velocity_m_s
        assert high.radial_flux_velocity_m_s < low.radial_flux_velocity_m_s
        assert high.provenance["electronegativity"] == 10.0


def test_provider_rejects_cross_temperature_mobility_reuse():
    condition = _condition()
    provider = LeeEconomouChlorineChargedTransportProvider(
        reduced_mobilities=_positive_mobilities(),
        ion_temperature=_scalar(
            300.0 * 1.380649e-23 / E_CHARGE_C,
            "eV",
        ),
    )
    with pytest.raises(ValueError, match="ion temperature"):
        provider.predict(condition, _densities(electronegativity=1.0))
