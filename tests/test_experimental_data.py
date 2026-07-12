from pathlib import Path

import numpy as np
import pytest

from petch.experimental_data import (
    load_bosch_wafer_measurements,
    load_krueger_2024_evidence,
)


DATA = (
    Path(__file__).parents[1]
    / "data"
    / "experimental"
    / "zenodo_17122442"
    / "Si_Oxide_etch_9_points.csv"
)
KRUEGER_DATA = Path(__file__).parents[1] / "data" / "experimental" / "krueger_2024"


def test_bosch_wafer_measurements_have_verified_provenance_and_units():
    rows = load_bosch_wafer_measurements(DATA)

    assert len(rows) == 684
    identified = [row for row in rows if row.experiment_key is not None]
    unidentified = [row for row in rows if row.experiment_key is None]
    assert len({(row.experiment_key, row.wafer_number) for row in identified}) == 75
    assert len(unidentified) == 9
    assert np.isclose(min(row.silicon_etch_um for row in rows), 38.2659)
    assert np.isclose(max(row.silicon_etch_um for row in rows), 42.8646)
    assert np.isclose(min(row.oxide_etch_um for row in rows), 0.5351)
    assert np.isclose(max(row.oxide_etch_um for row in rows), 0.7417)


def test_bosch_wafer_measurements_reject_unverified_content(tmp_path):
    altered = tmp_path / "measurements.csv"
    altered.write_bytes(DATA.read_bytes() + b"\n")

    with pytest.raises(ValueError, match="checksum mismatch"):
        load_bosch_wafer_measurements(altered)


def test_krueger_2024_keeps_calibration_measurements_and_simulated_inputs_separate():
    evidence = load_krueger_2024_evidence(KRUEGER_DATA)

    assert {metric.symbol for metric in evidence.calibration_metrics} == {
        "wm", "wt", "wf", "hf", "hm", "ah",
    }
    assert all(metric.evidence_type == "experiment" for metric in evidence.calibration_metrics)
    assert all(flux.evidence_type == "HPEM_simulation" for flux in evidence.boundary_fluxes)
    assert np.isclose(sum(flux.value_cm2_s for flux in evidence.boundary_fluxes), 3.984e17)
    assert all(item.split != "calibration" for item in evidence.transfer_observations)


def test_krueger_2024_exposes_held_out_experimental_transfer_trends():
    evidence = load_krueger_2024_evidence(KRUEGER_DATA)
    experimental = {
        (item.family, item.control, item.observable): item.value
        for item in evidence.transfer_observations
        if item.evidence_type in {"experiment", "experiment_and_simulation"}
    }

    assert experimental[("oxygen_ratio", "0.5", "feature_clogged")] == "true"
    assert experimental[("oxygen_ratio", "1.5", "etch_depth_rank")] == "maximum"
    assert experimental[("oxygen_ratio", "1.5_to_2.5", "etch_depth_increase")] == "false"
    assert experimental[("low_frequency_power_kw", "4_to_8", "final_profile_difference")] == (
        "few_differences")


def test_krueger_2024_rejects_unverified_transcription(tmp_path):
    target = tmp_path / "krueger_2024"
    target.mkdir()
    for source in KRUEGER_DATA.glob("*.csv"):
        (target / source.name).write_bytes(source.read_bytes())
    with (target / "transfer_observations.csv").open("ab") as stream:
        stream.write(b"\n")

    with pytest.raises(ValueError, match="checksum mismatch"):
        load_krueger_2024_evidence(target)
