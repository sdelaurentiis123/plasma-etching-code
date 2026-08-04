"""Gates for the community-standard CD triple (Top CD / Neck CD / neck depth).

Campaign 5 showed the legacy ``mask_opening_nm`` cannot separate an aperture
that is too narrow from one that necks at the wrong depth: it minimises over
the mask band and reports a single scalar.  These gates pin the new metrics on
synthetic interfaces whose neck width and depth are prescribed analytically.
"""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from petch.feature_geometry_state_3d import FeatureGeometry3D

from krueger_2024_trench_pilot import measure_krueger_metrics

DX = 0.01
NX, NY, NZ = 14, 3, 60
CENTER_UM = 0.065
SUBSTRATE_TOP = 0.20
MASK_TOP = 0.50
FLOOR = 0.12


def _build(mask_width, feature_width):
    """Assemble a two-material trench whose aperture follows ``mask_width(z)``."""
    x = np.arange(NX) * DX
    z = np.arange(NZ) * DX
    radial = np.abs(x - CENTER_UM)[:, None, None]
    z_grid = z[None, None, :]
    half_mask = 0.5 * np.vectorize(mask_width)(z)[None, None, :]
    half_feature = 0.5 * np.vectorize(feature_width)(z)[None, None, :]
    mask = np.minimum(
        np.minimum(radial - half_mask, z_grid - SUBSTRATE_TOP),
        MASK_TOP - z_grid)
    substrate = np.minimum(
        SUBSTRATE_TOP - z_grid,
        np.maximum(radial - half_feature, FLOOR - z_grid))
    mask = np.broadcast_to(mask, (NX, NY, NZ)).copy()
    substrate = np.broadcast_to(substrate, (NX, NY, NZ)).copy()
    phi = np.maximum(mask, substrate)
    material = np.where(phi < 0.0, 0, np.where(mask >= substrate, 2, 1))
    return FeatureGeometry3D(
        phi=phi, material_id=material, dx=DX, mesh_length_unit_m=1e-6,
        material_levelsets={1: substrate, 2: mask})


def _v_profile(top, neck, z_neck, ramp=0.06):
    def width(z):
        return min(top, neck + (top - neck) * abs(z - z_neck) / ramp)
    return width


def test_prescribed_neck_width_and_depth_recovered():
    """A V-shaped mask opening must report its analytic neck size and depth."""
    z_neck = 0.35
    geometry = _build(_v_profile(0.09, 0.04, z_neck), lambda z: 0.06)
    metrics = measure_krueger_metrics(
        geometry, substrate_top_um=SUBSTRATE_TOP, opening_center_um=CENTER_UM)
    assert metrics["mask_top_z_um"] == pytest.approx(MASK_TOP, abs=0.5 * DX)
    assert metrics["neck_cd_nm"] == pytest.approx(40.0, abs=1.0)
    assert metrics["neck_z_um"] == pytest.approx(z_neck, abs=DX)
    assert metrics["neck_depth_from_mask_top_nm"] == pytest.approx(
        (MASK_TOP - z_neck) * 1e3, abs=10.0)
    # Top CD samples the highest resolved mask plane, where the V is flat.
    assert metrics["top_cd_nm"] == pytest.approx(90.0, abs=1.0)


def test_monotone_taper_necks_at_the_bottom():
    """With no mask neck the minimum must be the deepest resolved aperture."""
    geometry = _build(
        lambda z: 0.05 + 0.08 * (z - SUBSTRATE_TOP) / (MASK_TOP - SUBSTRATE_TOP),
        lambda z: 0.02 + 0.03 * (z - FLOOR) / (SUBSTRATE_TOP - FLOOR))
    metrics = measure_krueger_metrics(
        geometry, substrate_top_um=SUBSTRATE_TOP, opening_center_um=CENTER_UM)
    profile = metrics["aperture_profile"]
    assert profile, "aperture profile must span the feature"
    assert metrics["neck_z_um"] == pytest.approx(
        min(item["z_um"] for item in profile), abs=1e-12)
    assert metrics["neck_z_um"] < SUBSTRATE_TOP
    assert metrics["top_cd_nm"] > metrics["neck_cd_nm"]


def test_legacy_fields_unchanged_and_consistent():
    """The legacy scalar stays the MASK-band minimum, not the global neck."""
    z_neck = 0.35
    geometry = _build(_v_profile(0.09, 0.04, z_neck), lambda z: 0.02)
    metrics = measure_krueger_metrics(
        geometry, substrate_top_um=SUBSTRATE_TOP, opening_center_um=CENTER_UM)
    mask_band = [
        item for item in metrics["aperture_profile"]
        if SUBSTRATE_TOP <= item["z_um"] <= metrics["mask_top_z_um"]]
    legacy = min(item["width_nm"] for item in mask_band)
    assert metrics["mask_opening_nm"] == pytest.approx(legacy, abs=1e-9)
    assert metrics["mask_opening_throat_z_um"] == pytest.approx(z_neck, abs=DX)
    # The narrower substrate trench owns the GLOBAL neck: the two disagree,
    # which is exactly the size-versus-location ambiguity the triple removes.
    assert metrics["neck_cd_nm"] < metrics["mask_opening_nm"]
    assert metrics["neck_z_um"] < SUBSTRATE_TOP
    for name in ("etch_depth_nm", "remaining_mask_thickness_nm",
                 "top_feature_width_nm", "maximum_feature_width_nm"):
        assert np.isfinite(metrics[name])
