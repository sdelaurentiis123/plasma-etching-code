import numpy as np
import pytest

from petch.product_ion_chlorine_si import (
    LeeChangProductIonSiSurfaceSensitivity,
    PRODUCT_ION_MASS_AMU,
)
from petch.reactor_global.wafer_sheath_transfer import (
    SpeciesResolvedIonEnergyDistribution,
)


def _distributions():
    return {
        name: SpeciesResolvedIonEnergyDistribution(
            species=name,
            ion_mass_amu=mass,
            flux_m2_s=1.0e18,
            energy_eV=np.full(16, 60.0),
            weight=np.full(16, 1.0 / 16.0),
        )
        for name, mass in PRODUCT_ION_MASS_AMU.items()
    }


def test_product_ion_limits_close_gross_minus_deposition_exactly():
    model = LeeChangProductIonSiSurfaceSensitivity()
    reflective = model.evaluate(
        _distributions(), chlorination_fraction=0.8, wall_limit="reflective")
    reactive = model.evaluate(
        _distributions(), chlorination_fraction=0.8, wall_limit="reactive")

    assert reflective.total_deposition_rate_si_m2_s == 0.0
    assert reflective.net_removal_rate_si_m2_s == pytest.approx(
        reflective.total_gross_removal_rate_si_m2_s)
    assert reactive.total_deposition_rate_si_m2_s == pytest.approx(1.0e18)
    assert reactive.net_removal_rate_si_m2_s == pytest.approx(
        reactive.total_gross_removal_rate_si_m2_s
        - reactive.total_deposition_rate_si_m2_s)
    assert reactive.net_removal_rate_si_m2_s < (
        reflective.net_removal_rate_si_m2_s)
    assert not reactive.supports_prediction
