"""Gates for the DECLARED ion-flux calibration (RESULTS_BLANKET_ANCHOR_2026-08-06).

The aggregate positive-ion flux is the one quantity in the removal chain with no
measurement behind it: the source publishes it as HPEM model output
(``evidence_type: HPEM_simulation``) and its reactor model carries no
experimental validation.  A normalization therefore acts there and never on the
beam-measured yields, which two independent experiments corroborate to 4.7%
(gray1993_mit, karahashi2007_hyomen).

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
    assert "beam-measured" in scaled.provenance["ion_flux_normalization_basis"]


@pytest.mark.parametrize("bad", [0.0, -1.0, np.nan, np.inf])
def test_invalid_normalization_refused(bad):
    deck = load_krueger_2024_reactor_flux_deck(DATA)
    with pytest.raises(ValueError):
        _normalize_ion_flux(deck, bad)


def test_measurement_only_lower_bound_on_the_ion_flux():
    """The bound that motivates the calibration, recomputed from measurements.

    825 nm / 60 s (SEM, Krueger Fig. 7b) over 2.2e28 m^-3 needs 3.025e20
    units/m2/s at the floor.  The measured per-ion ceiling for the most reactive
    single ion at >=1 keV, saturating, is 1.5 molecules/ion (Karahashi 2007
    Fig. 3).  Floor ion delivery at AR~21 is 0.70 (cascade scan).  Those three
    measured numbers alone put the wafer-plane ion flux at >= 2.4x the published
    HPEM value -- independent of any petch parameter.
    """
    removal = (825e-9 / 60.0) * 2.2e28
    assert removal == pytest.approx(3.025e20, rel=1e-3)
    floor_bound = removal / 1.5
    wafer_bound = floor_bound / 0.70
    assert wafer_bound / _ION_PUBLISHED_M2_S == pytest.approx(2.40, abs=0.02)


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
