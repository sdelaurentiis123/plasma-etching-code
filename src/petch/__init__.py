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
    BoundaryTransport3DResult, ChargedSurfaceReimpactPopulation3D,
    estimate_diffuse_form_factors_3d,
    merge_boundary_transport_results_3d,
    trace_charged_surface_events_field_3d,
    trace_boundary_state_field_3d,
    trace_boundary_state_first_hit_3d,
)
from .charging_poisson_3d import (
    NodalPoissonSystem3D, PoissonDiagnostics3D, assemble_q1_stiffness_3d,
    lump_triangle_sheet_charge_3d,
)
from .charging_coupled_3d import (
    CurrentBalanceMetrics3D,
    DielectricChargingConvergenceError, DielectricChargingStep3DResult,
    PhysicalTimeChargingIntegrationError, PhysicalTimeDielectricCharging3DResult,
    SteadyDielectricCharging3DResult,
    advance_dielectric_charging_3d, current_balance_metrics_3d,
    integrate_dielectric_charging_transient_3d, solve_dielectric_charging_steady_3d,
)
from .charged_surface_response_3d import (
    ChargedSurfaceContext3D, ChargedSurfaceResponse3D, ChargedSurfaceTransfer3D,
    OutgoingChargedParticleEvents3D, PerfectAbsorberChargedSurfaceResponse3D,
    GrazingSpecularIonReflection3D, Sobolewski2021ArKineticSEE3D,
    account_charged_surface_transfer_3d, perfect_absorber_surface_transfer_3d,
    sobolewski_2021_ar_kinetic_see_yield,
)
from .charged_surface_cascade_3d import (
    ChargedSurfaceCascade3DResult, augment_transport_with_charged_reimpacts_3d,
    solve_charged_surface_cascade_3d,
)
from .surface_charge_remap_3d import (
    SurfaceChargeRemap3DResult, remap_surface_charge_3d,
)
from .surface_interaction_table import (
    InteractionAxis, SurfaceInteractionDomainError, SurfaceInteractionEvaluation,
    SurfaceInteractionInterpolationAudit, SurfaceInteractionTable,
)
from .interaction_data import KounisMelas2024Tables, load_kounis_melas_2024_tables
from .tabulated_chemistry import (
    TabulatedSiClArMechanism, TabulatedSiPhysicalSputterMechanism,
    TabulatedSiPhysicalSputterStepResult, TabulatedSiSurfaceState,
    TabulatedSiSurfaceStepResult,
)
from .feature_step_3d import (
    FeatureGeometry3D, FeatureSolve3DResult, FeatureStep3DResult, FeatureStepValidity,
    advance_feature_step_3d, conservative_remap_surface_state,
    make_rectangular_trench_geometry_3d, solve_feature_3d,
)
from .validation_demo import (
    JEON_2022_DEMO_VERSION, Jeon2022DemoScore, Jeon2022DemoThresholds,
    Jeon2022Prediction, score_jeon_2022_demo,
)
from .neutral_radiosity_3d import (
    DiffuseFormFactors3D, DiffuseNeutralSolve3D, DiffuseSurfaceEmissionSolve3D,
    solve_diffuse_neutral_radiosity_3d, transport_diffuse_surface_emission_3d,
    transport_surface_product_population_3d,
)
from .experimental_boundary import (
    Jeon2022BoundaryClosure, build_jeon_2022_boundary_state,
)
from .physical_api import COMMON_FEATURE_ENGINE, PhysicalProcess, PhysicalResult
from .surface_exchange import (
    SurfaceMaterialExchange, SurfaceProductPopulation, unresolved_surface_exchange,
    validate_surface_product_routing,
)
from .physical_sputtering import (
    PhysicalSputterMechanism, PhysicalSputterParameters, PhysicalSputterState,
    PhysicalSputterStepResult,
)

# High-level 3D API (ViennaPS-shaped). Importing api pulls in the 3D engine (threed/warp).
from .api import Domain, SF6O2, Process, Result

__all__ = [
    # Explicit dimensional common engine and legacy compatibility interface.
    "COMMON_FEATURE_ENGINE", "PhysicalProcess", "PhysicalResult",
    "SurfaceMaterialExchange", "SurfaceProductPopulation", "unresolved_surface_exchange",
    "validate_surface_product_routing",
    "PhysicalSputterMechanism", "PhysicalSputterParameters", "PhysicalSputterState",
    "PhysicalSputterStepResult",
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
    "BoundaryTransport3DResult", "ChargedSurfaceReimpactPopulation3D",
    "merge_boundary_transport_results_3d", "trace_charged_surface_events_field_3d",
    "estimate_diffuse_form_factors_3d",
    "trace_boundary_state_field_3d",
    "trace_boundary_state_first_hit_3d",
    "NodalPoissonSystem3D", "PoissonDiagnostics3D", "assemble_q1_stiffness_3d",
    "lump_triangle_sheet_charge_3d",
    "CurrentBalanceMetrics3D", "DielectricChargingConvergenceError",
    "DielectricChargingStep3DResult",
    "PhysicalTimeChargingIntegrationError", "PhysicalTimeDielectricCharging3DResult",
    "SteadyDielectricCharging3DResult",
    "advance_dielectric_charging_3d", "current_balance_metrics_3d",
    "integrate_dielectric_charging_transient_3d", "solve_dielectric_charging_steady_3d",
    "ChargedSurfaceContext3D", "ChargedSurfaceResponse3D", "ChargedSurfaceTransfer3D",
    "OutgoingChargedParticleEvents3D", "PerfectAbsorberChargedSurfaceResponse3D",
    "GrazingSpecularIonReflection3D", "Sobolewski2021ArKineticSEE3D",
    "account_charged_surface_transfer_3d", "perfect_absorber_surface_transfer_3d",
    "sobolewski_2021_ar_kinetic_see_yield",
    "ChargedSurfaceCascade3DResult", "augment_transport_with_charged_reimpacts_3d",
    "solve_charged_surface_cascade_3d",
    "SurfaceChargeRemap3DResult", "remap_surface_charge_3d",
    "InteractionAxis", "SurfaceInteractionDomainError", "SurfaceInteractionEvaluation",
    "SurfaceInteractionInterpolationAudit", "SurfaceInteractionTable",
    "KounisMelas2024Tables", "load_kounis_melas_2024_tables",
    "TabulatedSiClArMechanism", "TabulatedSiPhysicalSputterMechanism",
    "TabulatedSiPhysicalSputterStepResult", "TabulatedSiSurfaceState",
    "TabulatedSiSurfaceStepResult",
    "FeatureGeometry3D", "FeatureSolve3DResult", "FeatureStep3DResult", "FeatureStepValidity",
    "advance_feature_step_3d", "conservative_remap_surface_state", "solve_feature_3d",
    "make_rectangular_trench_geometry_3d",
    "JEON_2022_DEMO_VERSION", "Jeon2022DemoScore", "Jeon2022DemoThresholds",
    "Jeon2022Prediction", "score_jeon_2022_demo",
    "DiffuseFormFactors3D", "DiffuseNeutralSolve3D", "DiffuseSurfaceEmissionSolve3D",
    "solve_diffuse_neutral_radiosity_3d", "transport_diffuse_surface_emission_3d",
    "transport_surface_product_population_3d",
    "Jeon2022BoundaryClosure", "build_jeon_2022_boundary_state",
]

__version__ = "0.2.0"
