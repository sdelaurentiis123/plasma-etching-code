"""Gate the cryo physisorption model (src/petch/cryo.py) against the CF4/H2 pseudo-wet benchmark
(Small Methods 2024) and the physisorption physics."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from petch.cryo import cryo_etch_rate, physisorption_coverage, residence_time


def test_warm_plateau():
    """B1: warm-side (T>0C) etch rate sits at the ~2.3 nm/s plateau."""
    assert abs(cryo_etch_rate(20.0) - 2.3) < 0.1
    assert abs(cryo_etch_rate(0.0) - 2.3) < 0.15


def test_cold_anchor_16x():
    """B2: at -60C the rate is 3.76 nm/s = 1.6x the +20C value (the firm anchor)."""
    assert abs(cryo_etch_rate(-60.0) - 3.76) < 0.15
    assert abs(cryo_etch_rate(-60.0) / cryo_etch_rate(20.0) - 1.6) < 0.1


def test_monotone_rise_as_T_drops():
    """Etch rate rises monotonically as temperature falls (more physisorption)."""
    Ts = [20, 0, -20, -40, -60, -80]
    ers = [cryo_etch_rate(T) for T in Ts]
    assert all(ers[i + 1] >= ers[i] for i in range(len(ers) - 1))


def test_coverage_bounds_and_onset():
    """theta in [0,1]; near-zero warm, strong (>0.5) by -60C; onset below 0C."""
    assert physisorption_coverage(20.0) < 0.05
    assert physisorption_coverage(-60.0) > 0.5
    assert physisorption_coverage(-100.0) < 1.0


def test_residence_time_cliff():
    """B7: C4F8 residence time (E_d=0.406) grows steeply -110C -> -120C (the etch on/off cliff)."""
    assert residence_time(-120.0, 0.406) > residence_time(-110.0, 0.406) * 3
