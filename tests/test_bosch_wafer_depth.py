from pathlib import Path

import numpy as np
import pytest

from petch.bosch_process_data import load_bosch_process_traces
from petch.bosch_wafer_depth import (
    build_bosch_reference_surface_mechanisms,
    predict_bosch_wafer_depth,
)
from petch.reactor_global.bosch_spts_reduced import (
    BoschSPTSReducedParameters,
    DeterministicBoschSPTSReactorToWafer,
)


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "experimental" / "zenodo_17122442"


@pytest.fixture(scope="module")
def default_prediction():
    trace = load_bosch_process_traces(
        DATA / "Process_data.nc", DATA / "Dictionary_process.nc")[0]
    boundary = DeterministicBoschSPTSReactorToWafer(
        BoschSPTSReducedParameters()).solve(trace)
    silicon, oxide = build_bosch_reference_surface_mechanisms()
    return predict_bosch_wafer_depth(boundary, silicon, oxide)


def test_default_measured_waveform_path_produces_a_conserved_radial_prediction(
        default_prediction):
    prediction = default_prediction

    assert prediction.wafer_mean_silicon_depth_m > 0.0
    assert prediction.wafer_mean_oxide_loss_m > 0.0
    assert prediction.silicon_to_oxide_selectivity > 0.0
    assert prediction.maximum_material_ledger_relative_residual < 1.0e-12
    assert np.any(prediction.remaining_film_units_m2 > 0.0)
    assert prediction.provenance["target_depth_used"] is False
    assert prediction.provenance["surface_parameters_shared_across_every_wafer"]


def test_reported_wafer_means_are_exact_annular_area_means(default_prediction):
    prediction = default_prediction
    area = prediction.radial_area_m2

    assert prediction.wafer_mean_silicon_depth_m == pytest.approx(
        np.dot(prediction.silicon_depth_m, area) / np.sum(area), rel=2.0e-15)
    assert prediction.wafer_mean_oxide_loss_m == pytest.approx(
        np.dot(prediction.oxide_loss_m, area) / np.sum(area), rel=2.0e-15)
    assert np.any(area == 0.0)
