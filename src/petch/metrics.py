"""Profile / depth metrics. Reused from run_benchmark.py (dimension-agnostic for the 2D phase)."""
import numpy as np


def ours_profile(seg, sub_top, W):
    """Etched-substrate profile from segments: x centered on the trench, depth positive-down."""
    cx = np.r_[seg[:, 0], seg[:, 2]]
    cy = np.r_[seg[:, 1], seg[:, 3]]
    sel = cy < sub_top - 0.02
    return cx[sel] - W / 2, sub_top - cy[sel]


def depth_centre(xc, dep, half=1.5):
    """Max etch depth within +/- half microns of the trench centre."""
    c = np.abs(xc) < half
    return dep[c].max() if c.any() else 0.0


def center_depth(result, half=1.5):
    """Center etch depth straight from a run_etch result dict."""
    seg = result['segs']
    if len(seg) == 0:
        return 0.0
    W = (result['xs'][-1] - result['xs'][0]) + result['dx']
    xc, dep = ours_profile(seg, result['sub_top'], W)
    return depth_centre(xc, dep, half=half)
