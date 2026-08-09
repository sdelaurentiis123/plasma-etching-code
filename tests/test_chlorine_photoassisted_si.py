import pytest

from petch import (
    Du2022ShortwavePhotoEtchYield,
    Hirsch2020PulsedDCAntiSynergySensitivity,
)


def test_du_yield_round_trip_and_measured_range():
    photon_flux = 7.6e16  # 7.6e12 cm^-2 s^-1
    for card in (
        Du2022ShortwavePhotoEtchYield.measured_lower_yield(),
        Du2022ShortwavePhotoEtchYield.measured_upper_yield(),
    ):
        velocity = card.etch_velocity_m_s(
            photon_flux, photon_wavelength_nm=106.0)
        assert card.required_photon_flux_m2_s(
            velocity, photon_wavelength_nm=106.0) == pytest.approx(photon_flux)
        assert card.supports_prediction


def test_du_yield_refuses_139nm_transfer_and_bad_card():
    card = Du2022ShortwavePhotoEtchYield.measured_lower_yield()
    with pytest.raises(ValueError, match="wavelength-specific"):
        card.etch_velocity_m_s(1.0e17, photon_wavelength_nm=139.0)
    with pytest.raises(ValueError):
        Du2022ShortwavePhotoEtchYield(silicon_atoms_per_photon=89.0)


def test_hirsch_digitized_curve_replays_knots_and_is_monotone():
    card = Hirsch2020PulsedDCAntiSynergySensitivity()
    replay = card.relative_yield(card.duty_cycle_percent)
    assert replay == pytest.approx(card.relative_pae_yield, abs=2e-15)
    fine = card.relative_yield([value / 10.0 for value in range(901)])
    assert all(left >= right for left, right in zip(fine[:-1], fine[1:]))
    assert not card.supports_prediction
    with pytest.raises(ValueError, match="outside digitized support"):
        card.relative_yield(95.0)
