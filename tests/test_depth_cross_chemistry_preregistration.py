from scripts.validate_depth_cross_chemistry_preregistration import (
    load_preregistration,
)


def test_depth_campaign_preregistration_is_checksum_bound_and_cross_chemistry():
    protocol, expected = load_preregistration()

    assert protocol.commit_sha256 == expected
    assert len(protocol.targets) == 86
    assert sum(
        item.split == "held_out_transfer" for item in protocol.targets
    ) == 71
    assert {
        item.chemistry_family
        for item in protocol.targets
        if item.split == "held_out_transfer"
    } == {
        "fluorocarbon_oxide",
        "cyclic_sf6_c4f8_silicon",
        "chlorine_silicon",
    }


def test_yoshie_feature_grid_is_frozen_before_value_reveal():
    protocol, _ = load_preregistration()
    yoshie = tuple(
        item for item in protocol.targets
        if item.observation_id.startswith("yoshie_fig")
    )

    assert len(yoshie) == 49
    assert all(item.split == "held_out_transfer" for item in yoshie)
    assert all(not hasattr(item, "value") for item in yoshie)
