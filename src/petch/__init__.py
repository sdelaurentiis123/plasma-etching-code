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
from .surface_kinetics import (
    EnergeticFlux, EnergeticYield, FaceResolvedEnergeticFlux, ParameterEvidence,
    ReducedSiO2FluorocarbonMechanism, ReducedSiO2FluorocarbonParameters,
    SiO2SurfaceState, SurfaceFluxes,
)
from .boundary_transport_3d import (
    BoundaryTransport3DResult, merge_boundary_transport_results_3d,
    trace_boundary_state_field_3d,
    trace_boundary_state_first_hit_3d,
)
from .charging_poisson_3d import (
    NodalPoissonSystem3D, PoissonDiagnostics3D, assemble_q1_stiffness_3d,
    lump_triangle_sheet_charge_3d,
)
from .charging_coupled_3d import (
    DielectricChargingConvergenceError, DielectricChargingStep3DResult,
    SteadyDielectricCharging3DResult, advance_dielectric_charging_3d,
    solve_dielectric_charging_steady_3d,
)
from .surface_interaction_table import (
    InteractionAxis, SurfaceInteractionDomainError, SurfaceInteractionEvaluation,
    SurfaceInteractionInterpolationAudit, SurfaceInteractionTable,
)
from .interaction_data import KounisMelas2024Tables, load_kounis_melas_2024_tables
from .feature_step_3d import (
    FeatureGeometry3D, FeatureSolve3DResult, FeatureStep3DResult, FeatureStepValidity,
    advance_feature_step_3d, conservative_remap_surface_state, solve_feature_3d,
)

# High-level 3D API (ViennaPS-shaped). Importing api pulls in the 3D engine (threed/warp).
from .api import Domain, SF6O2, Process, Result

__all__ = [
    # high-level 3D API (the main public interface)
    "Domain", "SF6O2", "Process", "Result",
    # config + low-level (full control)
    "PAR", "Flags", "DEFAULT_FLAGS", "run_etch",
    "make_trench", "extract_surface", "orient_normals", "seg_in_mask", "profile_bottom",
    "mc_flux", "_trace", "surface_rate", "surface_rate_langmuir",
    "advect", "extend_velocity", "reinit",
    "ours_profile", "depth_centre", "center_depth",
    "EnergeticFlux", "EnergeticYield", "FaceResolvedEnergeticFlux", "ParameterEvidence",
    "ReducedSiO2FluorocarbonMechanism", "ReducedSiO2FluorocarbonParameters",
    "SiO2SurfaceState", "SurfaceFluxes",
    "BoundaryTransport3DResult", "merge_boundary_transport_results_3d",
    "trace_boundary_state_field_3d",
    "trace_boundary_state_first_hit_3d",
    "NodalPoissonSystem3D", "PoissonDiagnostics3D", "assemble_q1_stiffness_3d",
    "lump_triangle_sheet_charge_3d",
    "DielectricChargingConvergenceError", "DielectricChargingStep3DResult",
    "SteadyDielectricCharging3DResult", "advance_dielectric_charging_3d",
    "solve_dielectric_charging_steady_3d",
    "InteractionAxis", "SurfaceInteractionDomainError", "SurfaceInteractionEvaluation",
    "SurfaceInteractionInterpolationAudit", "SurfaceInteractionTable",
    "KounisMelas2024Tables", "load_kounis_melas_2024_tables",
    "FeatureGeometry3D", "FeatureSolve3DResult", "FeatureStep3DResult", "FeatureStepValidity",
    "advance_feature_step_3d", "conservative_remap_surface_state", "solve_feature_3d",
]

__version__ = "0.2.0"
