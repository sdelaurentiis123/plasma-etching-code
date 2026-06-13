"""petch — 2D feature-scale plasma-etch simulator (Phase 0).

Dimension-agnostic physics (chemistry, params, flags, harness) is shared with the eventual
3D port; only geometry / transport-traversal / levelset get a 3D sibling later.
"""
from .params import PAR, Flags, DEFAULT_FLAGS
from .driver import run_etch
from .geometry import (make_trench, extract_surface, orient_normals,
                       seg_in_mask, profile_bottom)
from .transport import mc_flux, _trace
from .chemistry import surface_rate, surface_rate_langmuir
from .levelset import advect, extend_velocity, reinit
from .metrics import ours_profile, depth_centre, center_depth

__all__ = [
    "PAR", "Flags", "DEFAULT_FLAGS", "run_etch",
    "make_trench", "extract_surface", "orient_normals", "seg_in_mask", "profile_bottom",
    "mc_flux", "_trace", "surface_rate", "surface_rate_langmuir",
    "advect", "extend_velocity", "reinit",
    "ours_profile", "depth_centre", "center_depth",
]

__version__ = "0.1.0"
