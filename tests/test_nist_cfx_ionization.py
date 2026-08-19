import csv

import pytest

from petch.reactor_global.cfx_electron_ionization import (
    DATA_DIRECTORY,
    build_nist_1996_cfx_ionization_network,
    load_nist_1996_cfx_ionization_curves,
)
from petch.reactor_global.network import RateContext
from scripts.extract_nist_cfx_ionization import main as extraction_main


def test_six_measured_curves_close_atoms_charge_and_have_rates():
    curves = load_nist_1996_cfx_ionization_curves()
    network = build_nist_1996_cfx_ionization_network()

    assert len(curves) == len(network.reactions) == 6
    network.assert_closed_conservation()
    context = RateContext(electron_temperature_eV=5.0)
    assert all(curve.rate_coefficient().coefficient_si(context) > 0 for curve in curves)
    assert {curve.relative_uncertainty for curve in curves} == {.15, .16, .18, .20}


def test_70eV_source_values_and_single_energy_firewall():
    curves = {
        (curve.target_neutral, curve.product_ion): curve
        for curve in load_nist_1996_cfx_ionization_curves()
    }
    at_70 = {}
    for key, curve in curves.items():
        at_70[key] = dict(zip(
            curve.electron_energy_eV, curve.cross_section_m2
        ))[70.0]

    assert at_70[("CF3", "CF3+")] == pytest.approx(.376e-20)
    assert at_70[("CF3", "CF2+")] == pytest.approx(.76e-20)
    assert at_70[("CF3", "CF+")] == pytest.approx(.68e-20)
    assert at_70[("CF2", "CF2+")] == pytest.approx(1.03e-20)
    assert at_70[("CF2", "CF+")] == pytest.approx(1.19e-20)
    assert at_70[("CF", "CF+")] == pytest.approx(1.25e-20)

    with (DATA_DIRECTORY / "table33_cf2_dissociative_ionization.csv").open(
        newline="", encoding="utf-8"
    ) as stream:
        anchors = [
            row for row in csv.DictReader(stream)
            if row["evidence_class"].startswith("single_energy")
        ]
    assert len(anchors) == 1
    assert ("CF2", "F+") not in curves


def test_committed_tables_replay(monkeypatch):
    monkeypatch.setattr("sys.argv", ["extract_nist_cfx_ionization.py", "--check"])
    extraction_main()
