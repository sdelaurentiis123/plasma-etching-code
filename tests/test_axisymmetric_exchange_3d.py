import numpy as np
import pytest

from petch.axisymmetric_exchange_3d import (
    build_cylinder_band_exchange,
    cylinder_clausing_transmission,
    santeler_transmission,
)


def test_cylinder_band_operator_is_reciprocal_closed_and_self_viewing():
    operator = build_cylinder_band_exchange(0.5, np.linspace(0.0, 3.0, 13))
    transfer = operator["transfer_fraction"]
    area = operator["band_area"]
    exchange = transfer * area[:, None]
    assert np.allclose(exchange, exchange.T, rtol=0.0, atol=1e-12)
    closure = (transfer.sum(axis=1) + operator["escape_bottom"]
               + operator["escape_top"])
    assert np.allclose(closure, 1.0, rtol=0.0, atol=1e-9)
    # Interior cylinder bands see themselves: strictly positive diagonal.
    assert np.all(np.diag(transfer) > 0.0)


@pytest.mark.parametrize("aspect_ratio, tolerance", [
    (1.0, 0.005), (10.0, 0.005), (50.0, 0.01), (100.0, 0.01)])
def test_clausing_transmission_matches_santeler(aspect_ratio, tolerance):
    # Santeler's closed form is itself only good to ~0.7 percent; the gate allows the
    # union of both error budgets.  The exact-algebra value at 100:1 reproduces the
    # industry figure of ~1.3 percent bottom flux.
    tau = cylinder_clausing_transmission(aspect_ratio)
    reference = santeler_transmission(aspect_ratio)
    assert abs(tau - reference) / reference < tolerance
    if aspect_ratio == 100.0:
        assert tau == pytest.approx(0.013, abs=0.0008)


def test_clausing_transmission_band_convergence():
    coarse = cylinder_clausing_transmission(10.0, bands=120)
    fine = cylinder_clausing_transmission(10.0, bands=480)
    assert abs(coarse - fine) / fine < 2e-3
