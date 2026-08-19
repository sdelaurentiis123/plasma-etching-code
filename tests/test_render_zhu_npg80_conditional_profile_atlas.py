import json

from scripts import render_zhu_npg80_conditional_profile_atlas as atlas


def test_profile_atlas_exactly_replays_all_frozen_endpoints():
    payload = atlas._load_source()
    svg = atlas._svg(payload)

    assert atlas.SVG_PATH.read_text(encoding="utf-8") == svg
    assert atlas.MANIFEST_PATH.read_text(encoding="utf-8") == (
        atlas._render_manifest(payload, svg)
    )


def test_profile_atlas_preserves_blind_claim_boundary():
    manifest = json.loads(atlas.MANIFEST_PATH.read_text(encoding="utf-8"))

    assert manifest["target_sem_used"] is False
    assert manifest["target_depth_used"] is False
    assert manifest["profile_count"] == 56
    assert manifest["cleared_profile_count"] == 50
    assert manifest["not_cleared_profile_count"] == 6
    assert manifest["axes"]["width_nm"] == list(atlas.WIDTHS)
    assert "not a validated Oxford surface law" in manifest["claim_boundary"]
