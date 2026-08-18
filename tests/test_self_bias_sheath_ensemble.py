import math

import numpy as np

from petch.reactor_global.oxford80_self_bias import (
    Oxford80RIECondition,
    build_oxford80_self_bias_transfer,
    load_oxford80_self_bias_evidence,
)
from petch.reactor_global.self_bias_sheath_ensemble import (
    build_collisionless_self_bias_sheath_ensemble,
)
from petch.reactor_global.wafer_sheath_transfer import (
    DiagnosticConditionedRFSheathTransfer,
)
from scripts.audit_zhu_npg80_self_bias import DEFAULT_EVIDENCE


def test_all_bias_histories_project_to_closed_quasisteady_sheath_snapshots():
    target = Oxford80RIECondition(
        tool_model="Oxford PlasmaPro NPG80 RIE",
        rf_power_W=150.0,
        pressure_mTorr=30.0,
        gas_flows_sccm={"CHF3": 55.0, "SF6": 5.0, "O2": 1.0},
        electrode_temperature_C=20.0,
        duration_s=1200.0,
    )
    bias = build_oxford80_self_bias_transfer(
        target, load_oxford80_self_bias_evidence(DEFAULT_EVIDENCE))
    sheath = DiagnosticConditionedRFSheathTransfer(
        ion_mass_amu={"CF3+": 69.0, "SF5+": 127.0, "O2+": 32.0},
        electrode_area_m2=0.0314,
        plasma_potential_eV=18.0,
        frequency_hz=13.56e6,
        phase_count=24,
        steps_per_period=48,
        steps_per_transit=48,
        source="Oxford-80 self-bias evidence sensitivity",
    )

    ensemble = build_collisionless_self_bias_sheath_ensemble(
        bias_transfer=bias,
        sheath_transfer=sheath,
        positive_ion_flux_m2_s={
            "CF3+": 5.0e18,
            "SF5+": 2.0e18,
            "O2+": 1.0e18,
        },
        electron_temperature_eV=3.0,
        electron_density_m3=2.0e15,
        neutral_gas_temperature_K=350.0,
        normalized_time_nodes=np.array([0.0, 0.5, 1.0]),
    )

    assert len(ensemble.snapshots_by_history) == 5
    assert all(len(values) == 3 for values in ensemble.snapshots_by_history.values())
    drift = ensemble.snapshots_by_history[
        "exact-NGP80 conditioning thresholds"]
    assert [item.bias_magnitude_V for item in drift] == [300.0, 250.0, 200.0]
    assert all(
        abs(item.projection.power_closure_relative_residual) < 1.0e-12
        for values in ensemble.snapshots_by_history.values()
        for item in values
    )
    expected_density = (
        30.0 * 0.133322368 / (1.380649e-23 * 350.0)
    )
    assert math.isclose(ensemble.neutral_number_density_m3, expected_density)
    assert ensemble.time_scale_separation_ratio > 1.0e10
    assert ensemble.molecular_collision_cross_sections_supplied is False
    assert ensemble.supports_collisional_target_iead is False
    assert ensemble.supports_absolute_depth_prediction is False
