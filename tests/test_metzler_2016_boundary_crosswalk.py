import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
AUDIT = (
    ROOT / "results" / "curated" / "metzler_2016_boundary_crosswalk"
    / "audit.json"
)


def _report():
    return json.loads(AUDIT.read_text(encoding="utf-8"))


def test_crosswalk_recovers_six_sio2_fluences_without_calling_them_measured():
    report = _report()
    assert report["audit_id"] == "METZLER-2016-BOUNDARY-CROSSWALK-R1"
    sio2 = [row for row in report["points"] if row["substrate"] == "SiO2"]
    assert len(sio2) == 6
    assert report["gates"][
        "six_SiO2_conditions_recover_duration_linear_fluence"
    ]["passes"]
    assert any(
        "independently measured ion current" in claim
        for claim in report["claim_boundary"]["not_established"]
    )


def test_crosswalk_preserves_duplicate_si_marker_and_density_conditioning():
    report = _report()
    duplicate = [
        row for row in report["points"]
        if row["substrate"] == "Si"
        and row["energy_eV"] == 30
        and row["etch_step_s"] == 40.0
    ]
    assert {row["replicate"] for row in duplicate} == {1, 2}
    assert report["density_conventions"]["SiO2_formula_units_m3"] == 2.2e28
    assert report["material_metrics"]["SiO2"][
        "mean_author_normalized_flux_m2_s"
    ] > 1.7e20
    assert report["material_metrics"]["SiO2"][
        "mean_author_normalized_flux_m2_s"
    ] < 1.9e20
