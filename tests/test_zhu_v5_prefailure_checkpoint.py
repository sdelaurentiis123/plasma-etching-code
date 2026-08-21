from hashlib import sha256
import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
CHECKPOINT = ROOT / "tests" / "data" / "zhu_v5_prefailure_8ccd278.npz"
CHECKPOINT_SHA256 = (
    "95d43ec8df8372ef14df57cd28df93ac90fe5c2e72da84cceaf66965f7c50998"
)


def test_zhu_v5_prefailure_checkpoint_is_pinned_and_self_describing():
    assert sha256(CHECKPOINT.read_bytes()).hexdigest() == CHECKPOINT_SHA256

    with np.load(CHECKPOINT, allow_pickle=False) as stored:
        metadata = json.loads(str(stored["metadata_json"]))
        assert stored["phi"].shape == (41, 41, 106)
        assert stored["material_id"].shape == stored["phi"].shape
        assert metadata["material_levelset_ids"] == [1, 2, 3]
        assert metadata["elapsed_s"] == 719.8619631901854
        assert metadata["accepted_steps"] == 149
        assert metadata["step_duration_s"] == 4.831288343558282
        assert metadata["surface_state_names"] == [
            "m1__removed_material_units_m2",
            "m2__removed_material_units_m2",
        ]
