from pathlib import Path

import numpy as np
import pytest

from petch.interaction_data import load_kounis_melas_2024_tables
from petch.surface_interaction_table import SurfaceInteractionDomainError


DATA = (
    Path(__file__).parents[1] / "data" / "surface_interactions" / "kounis_melas_2024")


def test_deepmd_si_cl_ar_tables_replay_archived_nodes_and_uncertainty():
    tables = load_kounis_melas_2024_tables(DATA)
    sputter = tables.sputtering.evaluate({"ion_energy": np.array([50.0, 100.0, 200.0])})
    rie = tables.reactive_ion_etch.evaluate({
        "cl2_to_ar_flux_ratio": np.array([10.0, 50.0, 100.0, 200.0])})

    assert np.allclose(
        sputter.values["physical_sputter_yield"],
        [0.009433962264150945, 0.0625, 0.20754716981132074])
    assert np.allclose(
        sputter.standard_uncertainty["amorphous_layer_thickness"],
        [0.8181278173348, 1.083936648716531, 1.586888329896573])
    assert np.allclose(
        rie.values["reactive_etch_yield"],
        [0.24182079610957588, 0.37330997562376655,
         0.5007335313599666, 0.6061765299147771])
    assert tables.sputtering.provenance["evidence_type"] == "DeepMD_molecular_dynamics"


def test_deepmd_ale_product_table_retains_species_resolved_yields():
    table = load_kounis_melas_2024_tables(DATA).ale_products
    dosage = table.axes[0].values
    evaluated = table.evaluate({"ar_ion_dosage": dosage})

    assert np.isclose(evaluated.values["sicl_yield"][0], 0.07833333333333332)
    assert np.isclose(evaluated.values["sicl2_yield"][0], 0.065)
    assert np.all(evaluated.values["sicl2_yield"][3:] == 0.0)
    assert evaluated.standard_uncertainty["cl_yield"].shape == dosage.shape


def test_deepmd_tables_refuse_unvalidated_energy_and_flux_ratio_extrapolation():
    tables = load_kounis_melas_2024_tables(DATA)
    with pytest.raises(SurfaceInteractionDomainError, match="ion_energy"):
        tables.sputtering.evaluate({"ion_energy": 300.0})
    with pytest.raises(SurfaceInteractionDomainError, match="cl2_to_ar_flux_ratio"):
        tables.reactive_ion_etch.evaluate({"cl2_to_ar_flux_ratio": 5.0})


def test_deepmd_leave_one_out_error_is_reported_separately_from_md_uncertainty():
    tables = load_kounis_melas_2024_tables(DATA)
    sputter = tables.sputtering.leave_one_out_interpolation_audit(
        "physical_sputter_yield")
    rie = tables.reactive_ion_etch.leave_one_out_interpolation_audit(
        "reactive_etch_yield")

    assert np.array_equal(sputter.coordinates, [100.0])
    assert sputter.absolute_error[0] > 0.0
    assert np.array_equal(rie.coordinates, [50.0, 100.0])
    assert np.all(rie.absolute_error > 0.0)
    source_uncertainty = tables.reactive_ion_etch.standard_uncertainty[
        "reactive_etch_yield"][1:-1]
    assert not np.array_equal(rie.absolute_error, source_uncertainty)


def test_deepmd_tables_reject_modified_primary_data(tmp_path):
    target = tmp_path / "tables"; target.mkdir()
    for source in DATA.glob("*.csv"):
        (target / source.name).write_bytes(source.read_bytes())
    with (target / "RIE.csv").open("ab") as stream:
        stream.write(b"\n")
    with pytest.raises(ValueError, match="checksum mismatch"):
        load_kounis_melas_2024_tables(target)
