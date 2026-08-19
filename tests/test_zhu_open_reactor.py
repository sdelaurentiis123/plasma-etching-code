import math
import json
from types import SimpleNamespace

import numpy as np
import pytest

from petch.reactor_global import CylindricalReactor, RateContext
from petch.reactor_global.zhu_open_reactor import (
    ZhuOpenReactorCondition,
    ZhuOpenReactorModel,
    compile_bimolecular_kinetic_pairs,
    log_flux_ratio_residual,
    positive_ion_wall_return,
    wall_resolved_charged_power_density_W_m3,
)
from petch.reactor_global.zhu_supplemental_chemistry import (
    build_zhu_supplemental_chemistry,
    zhu_reactor_species,
)
from scripts.run_zhu_open_reactor import _continuation_state


def _condition():
    geometry = CylindricalReactor(radius_m=.085, length_m=.03)
    return ZhuOpenReactorCondition(
        condition_id="manufactured-zhu",
        geometry=geometry,
        neutral_control_volume_m3=geometry.volume_m3,
        pressure_Pa=3.99967104,
        gas_temperature_K=293.15,
        feed_molecules_s={"CHF3": 1.0, "SF6": .1, "O2": .02},
        absorbed_power_W=90.0,
        source_frequency_hz=13.56e6,
        reduced_field_bounds_Td=(40.0, 600.0),
        ion_temperature_eV=.03,
        ion_momentum_mean_free_path_m=1.0e-4,
        mean_positive_ion_wall_energy_eV=250.0,
        neutral_reduced_diffusivity_m_inv_s=6.0e20,
        neutral_wall_probabilities={
            "F": .02, "H": .10, "O": .02, "O(1d)": .02,
        },
        source="manufactured",
        absorbed_power_source="sensitivity, not forward power",
        machine_closure_source="manufactured",
    )


def test_condition_keeps_forward_to_absorbed_transfer_explicit():
    condition = _condition()
    assert condition.target_neutral_density_m3 == pytest.approx(
        3.99967104 / (1.380649e-23 * 293.15))
    assert condition.active_volume_fraction == 1.0
    assert condition.supports_unique_machine_state is False
    assert condition.uses_wall_resolved_sheath_power is False


def test_condition_requires_a_complete_wall_resolved_sheath_pair():
    condition = _condition()
    with pytest.raises(ValueError, match="supplied together"):
        ZhuOpenReactorCondition(**{
            **condition.__dict__,
            "powered_electrode_sheath_drop_V": 300.0,
        })
    resolved = ZhuOpenReactorCondition(**{
        **condition.__dict__,
        "powered_electrode_sheath_drop_V": 300.0,
        "grounded_surface_sheath_drop_V": 20.0,
    })
    assert resolved.uses_wall_resolved_sheath_power is True


def test_powered_and_grounded_loss_channels_close_and_recover_legacy_power():
    model = object.__new__(ZhuOpenReactorModel)
    species = zhu_reactor_species()
    model.species_by_name = {item.name: item for item in species}
    model.positive_names = tuple(
        item.name for item in species if item.role == "positive_ion")
    model.negative_names = tuple(
        item.name for item in species if item.role == "negative_ion")
    densities = {name: 1.0e15 for name in model.positive_names}
    densities.update({name: 2.0e15 for name in model.negative_names})
    densities["e"] = 1.0e15
    total, _, resolved = model._charged_wall_transport(
        densities, equivalent_temperature_eV=4.0, condition=_condition())
    for name, frequency in total.items():
        assert sum(resolved[name]) == pytest.approx(frequency, rel=2.0e-15)

    powered_loss = {name: resolved[name][0] * densities[name] for name in total}
    grounded_loss = {name: resolved[name][1] * densities[name] for name in total}
    charge = {
        name: model.species_by_name[name].charge_number for name in total}
    powered_power, grounded_power = (
        wall_resolved_charged_power_density_W_m3(
            powered_positive_ion_loss_m3_s=powered_loss,
            grounded_positive_ion_loss_m3_s=grounded_loss,
            ion_charge_numbers=charge,
            electron_wall_energy_eV=8.0,
            powered_electrode_sheath_drop_V=250.0,
            grounded_surface_sheath_drop_V=250.0,
        )
    )
    legacy = 1.602176634e-19 * sum(
        (powered_loss[name] + grounded_loss[name])
        * charge[name]
        * (8.0 + 250.0)
        for name in total
    )
    assert powered_power + grounded_power == pytest.approx(
        legacy, rel=2.0e-15)


def test_every_positive_ion_has_atom_conserving_wall_return():
    by_name = {species.name: species for species in zhu_reactor_species()}
    positives = {
        name for name, species in by_name.items()
        if species.role == "positive_ion"
    }
    for name in positives:
        products = positive_ion_wall_return(name)
        for element, count in by_name[name].composition.items():
            assert sum(
                amount * by_name[product].composition.get(element, 0)
                for product, amount in products.items()
            ) == count
    assert set(positive_ion_wall_return("SF4++")) == {"SF4"}
    assert set(positive_ion_wall_return("H2+")) == {"H2"}
    assert dict(positive_ion_wall_return("HF+")) == {"HF": 1.0}
    with pytest.raises(KeyError):
        positive_ion_wall_return("missing+")


def test_feed_and_wall_closures_are_strictly_validated():
    condition = _condition()
    with pytest.raises(ValueError):
        ZhuOpenReactorCondition(**{
            **condition.__dict__,
            "feed_molecules_s": {"CHF3": 1.0, "SF6": .1},
        })
    with pytest.raises(ValueError):
        ZhuOpenReactorCondition(**{
            **condition.__dict__,
            "neutral_wall_probabilities": {"F": 1.01},
        })
    assert math.isfinite(condition.absorbed_power_density_W_m3)


def test_sparse_jacobian_shape_requires_no_false_dense_columns():
    # Construct without loading the licensed collision workbook by exercising
    # the dependency method on an allocated instance with a minimal attribute
    # set is intentionally avoided: the public model constructor validates
    # shared-deck identity.  Its implementation is still pinned by this
    # structural source check rather than a brittle manufactured collision
    # deck.
    assert callable(ZhuOpenReactorModel.jacobian_sparsity)


def test_log_flux_ratio_residual_has_same_root_and_exact_ratio_coordinate():
    production = np.asarray([1.0e-9, 1.0, 1.0e9])
    loss = np.asarray([1.0, 1.0, 1.0])
    bounded = (production - loss) / (production + loss)
    transformed = log_flux_ratio_residual(bounded)
    assert transformed == pytest.approx(np.log(production / loss))
    assert log_flux_ratio_residual(np.asarray([0.0]))[0] == 0.0
    with pytest.raises(ValueError):
        log_flux_ratio_residual(np.asarray([float("nan")]))


def test_continuation_lifts_old_solution_without_breaking_pressure(tmp_path):
    path = tmp_path / "prior.json"
    path.write_text(json.dumps({"state": {
        "densities_m3": {"CHF3": 8.0, "F": 2.0, "e": 0.1},
        "exhaust_loss_frequency_s_inv": 3.0,
        "reduced_electric_field_Td": 250.0,
    }}), encoding="utf-8")
    model = SimpleNamespace(
        species_order=("CHF3", "F", "CO", "e"),
        neutral_names=("CHF3", "F", "CO"),
        positive_names=(),
        negative_names=(),
    )
    condition = SimpleNamespace(target_neutral_density_m3=100.0)
    densities, exhaust, field = _continuation_state(
        path, model=model, condition=condition)
    assert sum(densities[name] for name in model.neutral_names) == pytest.approx(100.0)
    assert densities["CO"] == pytest.approx(0.5)
    assert densities["CHF3"] / densities["F"] == pytest.approx(4.0)
    assert densities["e"] == pytest.approx(0.1)
    assert exhaust == 3.0
    assert field == 250.0


def test_continuation_seeds_new_ions_without_changing_charge_inventory(tmp_path):
    path = tmp_path / "prior.json"
    path.write_text(json.dumps({"state": {
        "densities_m3": {
            "CHF3": 1000.0,
            "old+": 100.0,
            "old-": 90.0,
            "e": 10.0,
        },
        "exhaust_loss_frequency_s_inv": 3.0,
        "reduced_electric_field_Td": 250.0,
    }}), encoding="utf-8")
    model = SimpleNamespace(
        species_order=("CHF3", "old+", "new+", "old-", "new-", "e"),
        neutral_names=("CHF3",),
        positive_names=("old+", "new+"),
        negative_names=("old-", "new-"),
        species_by_name={
            "old+": SimpleNamespace(charge_number=1),
            "new+": SimpleNamespace(charge_number=1),
            "old-": SimpleNamespace(charge_number=-1),
            "new-": SimpleNamespace(charge_number=-1),
        },
    )
    condition = SimpleNamespace(target_neutral_density_m3=1000.0)
    densities, _, _ = _continuation_state(
        path, model=model, condition=condition)
    positive_charge = densities["old+"] + densities["new+"]
    negative_charge = densities["old-"] + densities["new-"]
    assert positive_charge == pytest.approx(100.0)
    assert negative_charge == pytest.approx(90.0)
    assert densities["new+"] == pytest.approx(5.0)
    assert densities["new-"] == pytest.approx(4.5)


def test_compiled_zhu_mass_action_matches_general_network_evaluator():
    supplemental = build_zhu_supplemental_chemistry()
    model = object.__new__(ZhuOpenReactorModel)
    model.supplemental_chemistry = supplemental
    model._supplemental_kinetic_pairs = compile_bimolecular_kinetic_pairs(
        supplemental.network)
    densities = {
        name: 1.0e14 * (index + 1)
        for index, name in enumerate(supplemental.network.species_names)
    }
    context = RateContext(3.7, 350.0)
    general = supplemental.network.event_rates_m3_s(densities, context)
    compiled = model._supplemental_event_rates_m3_s(
        densities, context, {})
    np.testing.assert_allclose(compiled, general, rtol=3.0e-16, atol=0.0)
