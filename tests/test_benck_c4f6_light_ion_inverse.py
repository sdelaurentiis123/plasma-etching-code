from hashlib import sha256
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
AUDIT = (
    ROOT / "results" / "curated" / "benck_c4f6_light_ion_inverse_v1"
    / "audit.json"
)


def _load(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def test_benck_light_ion_inverse_is_target_free_and_checksum_bound():
    audit = _load(AUDIT)

    assert audit["calibration_firewall"] == {
        "absolute_current_scale_used": False,
        "feature_depth_used": False,
        "krueger_825_nm_used": False,
        "reaction_parameter_fitted": False,
    }
    assert audit["certification"]["supports_absolute_reactor_flux"] is False
    assert audit["certification"]["supports_feature_depth"] is False
    sources = audit["sources"]
    assert sha256((
        ROOT / "data" / "experimental" / "benck_2003_c4f6"
        / "figure9_mass_resolved_ion_current.csv"
    ).read_bytes()).hexdigest() == sources["benck_feed_csv_sha256"]
    assert sha256((
        ROOT / "data" / "experimental" / "benck_2003_c4f6"
        / "figure14a_cf2_cf_feed_ratio.csv"
    ).read_bytes()).hexdigest() == sources["benck_neutral_ratio_csv_sha256"]
    assert sha256((
        ROOT / "data" / "experimental" / "benck_2003_c4f6"
        / "figure14a_digitization_manifest.json"
    ).read_bytes()).hexdigest() == (
        sources["benck_neutral_ratio_manifest_sha256"]
    )


def test_common_loss_closure_fails_without_clipping_negative_density():
    audit = _load(AUDIT)
    decision = audit["physics_decision"]
    diagnostics = audit["diagnostics"]

    assert decision["common_loss_source_only_closure_is_adequate"] is False
    assert decision["negative_density_is_accepted_as_physical"] is False
    assert decision["negative_density_is_clipped"] is False
    assert diagnostics["ar_containing_condition_count"] == 3
    assert diagnostics["ar_containing_conditions_with_nonnegative_solution"] == 0
    assert diagnostics[
        "ar_containing_conditions_with_nonnegative_solution_when_"
        "CF2_CF_ratio_is_halved"
    ] == 0
    assert diagnostics["maximum_ratio_replay_absolute_error"] < 2.0e-15


def test_inverse_declares_the_required_reactor_rung_without_overclaim():
    audit = _load(AUDIT)

    required = audit["physics_decision"]["required_next_operators"]
    assert "species- and mass-dependent Bohm/wall/exhaust loss" in required
    assert "ion-neutral conversion and charge exchange" in required
    assert audit["certification"]["identifies_unique_missing_operator"] is False
    assert audit["certification"]["supports_steady_reactor_composition"] is False
    assert audit["certification"]["supports_krueger_boundary"] is False
