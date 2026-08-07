import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
AUDIT = (
    ROOT / "results" / "curated" / "metzler_2016_surface_closure"
    / "finite_range_audit.json"
)


def _report():
    return json.loads(AUDIT.read_text(encoding="utf-8"))


def test_finite_range_replacement_is_no_fit_and_still_falsified():
    report = _report()
    assert report["audit_id"] == (
        "METZLER-2016-FINITE-RANGE-SURFACE-CLOSURE-R1"
    )
    assert report["transport_model"]["finite_range"]
    assert report["transport_model"]["free_attenuation_parameters"] == []
    assert not report["verdict"][
        "finite_range_transport_repairs_surface_response"
    ]
    assert report["verdict"][
        "all_maximum_energy_replays_still_overpredict"
    ]
    assert report["verdict"][
        "smallest_maximum_energy_overprediction_factor"
    ] > 5.0


def test_spectrum_sensitivity_does_not_manufacture_an_iedf():
    report = _report()
    assert not report["verdict"][
        "iedf_sensitivity_identifies_one_boundary"
    ]
    cases = {
        row["spectrum_sensitivity_case"] for row in report["metrics"]
    }
    assert cases == {
        "delta_at_reported_maximum",
        "arcsine_0_to_reported_maximum",
        "arcsine_14eV_to_reported_maximum",
    }
    assert report["integration_convergence"][
        "maximum_relative_prediction_change"
    ] < 5e-3
