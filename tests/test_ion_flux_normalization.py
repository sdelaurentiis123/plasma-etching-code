"""Gates for the explicit, target-fitted ion-flux normalization.

The aggregate positive-ion flux is the one quantity in the removal chain with no
measurement behind it: the source publishes it as HPEM model output
(``evidence_type: HPEM_simulation``) and its reactor model carries no
experimental validation. A normalization may be explored there without altering
the other published boundary values, but it is not independently inferred:
Krueger publishes neither ion-species composition nor a blanket etch datum.

These gates pin the three properties that keep the calibration honest:
  1. default 1.0 is bitwise inert -- no existing result can move;
  2. it scales ions ONLY, leaving every neutral flux untouched;
  3. a calibrated flux carries the ``declared_calibration`` evidence kind, so it
     can never be audited as a measured or model-published value.
"""
import sys

import numpy as np
import pytest

sys.path.insert(0, "src")

from petch.reactor_boundary import (
    _normalize_ion_flux, load_krueger_2024_reactor_flux_deck)

DATA = "data/experimental/krueger_2024"

# krueger-2024.txt L298, Table I: "Ions  1.2 x 10^16" cm-2 s-1.
_ION_PUBLISHED_M2_S = 1.2e20


def _by_name(deck):
    return {item.name: item for item in deck.species_fluxes}


def test_published_ion_flux_is_the_table_i_value():
    deck = load_krueger_2024_reactor_flux_deck(DATA)
    ion = _by_name(deck)["ions"]
    assert ion.flux_m2_s == pytest.approx(_ION_PUBLISHED_M2_S, rel=1e-12)
    assert ion.evidence_kind == "HPEM_simulation"


def test_unit_normalization_is_inert():
    """Default 1.0 must leave the deck bitwise identical -- no silent drift."""
    deck = load_krueger_2024_reactor_flux_deck(DATA)
    same = _normalize_ion_flux(deck, 1.0)
    for name, record in _by_name(deck).items():
        assert _by_name(same)[name].flux_m2_s == record.flux_m2_s
        assert _by_name(same)[name].evidence_kind == record.evidence_kind


@pytest.mark.parametrize("factor", [1.5, 2.4, 2.8, 4.0])
def test_scales_ions_only_and_leaves_neutrals_bitwise(factor):
    deck = load_krueger_2024_reactor_flux_deck(DATA)
    scaled = _normalize_ion_flux(deck, factor)
    base, new = _by_name(deck), _by_name(scaled)
    assert new["ions"].flux_m2_s == pytest.approx(
        base["ions"].flux_m2_s * factor, rel=1e-12)
    for name, record in base.items():
        if record.role == "positive_ion_mixture":
            continue
        assert new[name].flux_m2_s == record.flux_m2_s, name


def test_calibrated_flux_is_labelled_as_calibration():
    """An audit must never mistake a calibrated flux for a published one."""
    scaled = _normalize_ion_flux(load_krueger_2024_reactor_flux_deck(DATA), 2.8)
    ion = _by_name(scaled)["ions"]
    assert ion.evidence_kind == "declared_calibration"
    assert scaled.provenance["ion_flux_normalization"] == pytest.approx(2.8)
    assert "target fit" in scaled.provenance["ion_flux_normalization_basis"]


@pytest.mark.parametrize("bad", [0.0, -1.0, np.nan, np.inf])
def test_invalid_normalization_refused(bad):
    deck = load_krueger_2024_reactor_flux_deck(DATA)
    with pytest.raises(ValueError):
        _normalize_ion_flux(deck, bad)


def test_karahashi_1000ev_marker_is_not_a_universal_ion_yield_ceiling():
    """The former 2.40x "measurement-only bound" used a false premise.

    Figure 4's same CF3+ series rises above the rounded 1.5 value at 1000 eV.
    More importantly, it ends at 2000 eV and contains no simultaneous molecule
    flux, while Krueger's aggregate IEAD extends to 4800 eV. It cannot supply a
    universal per-ion ceiling or an ion-flux lower bound for that reactor.
    """
    from pathlib import Path
    from petch.experimental_data import load_karahashi_2007_reactive_ion_yields

    path = (Path(__file__).parents[1] / "data" / "experimental"
            / "karahashi_2007" / "figure4_reactive_ion_yields.csv")
    rows = load_karahashi_2007_reactive_ion_yields(path)
    cf3 = [row for row in rows if row.species == "CF3+"]
    at_1000 = next(row for row in cf3 if row.energy_eV == 1000.0)
    assert at_1000.yield_sio2_per_ion == pytest.approx(1.4703)
    assert max(row.yield_sio2_per_ion for row in cf3) == pytest.approx(1.8736)


@pytest.mark.parametrize("extra", [
    dict(low_frequency_power_kw=6.0, oxygen_to_fluorocarbon_ratio=0.5),
    dict(low_frequency_power_kw=6.0, oxygen_to_fluorocarbon_ratio=1.5),
    dict(low_frequency_power_kw=6.0, oxygen_to_fluorocarbon_ratio=2.5),
    dict(low_frequency_power_kw=4.0),
    dict(low_frequency_power_kw=8.0),
])
def test_every_transfer_condition_carries_the_same_calibration(extra):
    """Regression: the Figure-16 oxygen table supplies its own fluxes and once
    discarded the normalization, so the oxygen conditions ran uncalibrated while
    the power conditions ran calibrated -- which would have silently corrupted
    the out-of-sample scorecard.  Every condition must scale identically.
    """
    from petch.reactor_boundary import build_krueger_2024_transfer_boundary

    kw = dict(reference_plane_m=3e-6, neutral_direction_polar_order=12,
              neutral_direction_azimuthal_order=24,
              ion_azimuthal_closure="axisymmetric_uniform",
              ion_azimuthal_order=16)
    plain = build_krueger_2024_transfer_boundary(DATA, **extra, **kw)
    scaled = build_krueger_2024_transfer_boundary(
        DATA, ion_flux_normalization=2.8, **extra, **kw)

    def ion(boundary):
        return [s.flux_m2_s for s in boundary.species if s.charge_number != 0][0]

    def neutrals(boundary):
        return [s.flux_m2_s for s in boundary.species if s.charge_number == 0]

    assert ion(scaled) == pytest.approx(ion(plain) * 2.8, rel=1e-12)
    assert neutrals(scaled) == neutrals(plain)
