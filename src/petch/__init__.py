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
    EnergeticFlux, EnergeticYield, FaceResolvedEnergeticFlux,
    LowEnergyActivationYield, ParameterEvidence,
    ReducedSiO2FluorocarbonMechanism, ReducedSiO2FluorocarbonParameters,
    SiO2SurfaceState, SteinbruchelYield, SurfaceFluxes,
)
from .fluorocarbon_lamagna import (
    LaMagnaFluorocarbonParameters, LaMagnaFluorocarbonState,
    LaMagnaFluorocarbonStepResult, LaMagnaGarozzoFluorocarbonMechanism,
)
from .boundary_transport_3d import (
    BoundaryTransport3DResult, ChargedSurfaceReimpactPopulation3D,
    average_boundary_transport_results_3d,
    estimate_diffuse_form_factors_3d,
    merge_boundary_transport_results_3d,
    trace_charged_surface_events_field_3d,
    trace_boundary_state_field_3d,
    trace_boundary_state_first_hit_3d,
)
from .charging_poisson_3d import (
    NodalPoissonSystem3D, PoissonDiagnostics3D, assemble_q1_stiffness_3d,
    lump_mixed_surface_density_3d, lump_triangle_sheet_charge_3d,
)
from .conductor_terminal_3d import (
    ConductorTerminalCurrent3D, RemotePadElectronCollector3D,
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
    ChargedSurfaceCascade3DResult, apply_charged_surface_response_to_transport_3d,
    augment_transport_with_charged_reimpacts_3d, derived_tail_bounce_budget_3d,
    solve_charged_surface_cascade_3d,
)
from .surface_charge_remap_3d import (
    SurfaceChargeRemap3DResult, remap_surface_charge_3d,
)
from .charging_coevolution_3d import (
    CHARGING_RUN_MANIFEST_SCHEMA, ChargingCoevolution3DResult,
    ChargingCoevolutionStep3DResult,
    ExperimentalObservableTolerance3D, PhysicalPatchBalance3D, ResolvedBiasSegment3D,
    SurfaceChargingSaturation3DResult, SurfaceChargingSaturationError,
    integrate_surface_charging_to_saturation_3d, physical_surface_patch_groups_3d,
    solve_charging_coevolution_3d,
)
from .charging_stationarity_3d import (
    PROFILE_STATIONARITY_CONTRACT_DRAFT,
    ProfileChargingStationarity3DResult,
    ProfileChargingStationarityBlock3D,
    ProfileChargingStationarityContract3D,
    assess_profile_charging_stationarity_3d,
)
from .charging_checkpoint_3d import (
    CHARGING_CHECKPOINT_SCHEMA, PhysicalChargingCheckpoint3D,
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
    HWANG_GIAPIS_1997_EEDF_SHA256, HWANG_GIAPIS_1997_IEDF_SHA256,
    HWANG_GIAPIS_1997_PDF_SHA256,
    Jeon2022BoundaryClosure, Jeong2023IonBoundaryClosure,
    build_hwang_giapis_1997_boundary_state,
    build_jeon_2022_boundary_state, build_jeong_2023_boundary_state,
)
from .reactor_boundary import (
    PlasmaDiagnosticState, ReactorSpeciesFlux, TabulatedReactorFluxDeck,
    build_diagnostic_virtual_sheath_boundary, build_tabulated_reactor_boundary,
    load_krueger_2024_reactor_flux_deck,
)
from .physical_api import (
    COMMON_CHARGING_ENGINE, COMMON_CHARGING_ENSEMBLE_ENGINE, COMMON_FEATURE_ENGINE,
    COMMON_FEATURE_MANIFEST_SCHEMA,
    PhysicalChargingEnsembleProcess, PhysicalChargingEnsembleResult,
    PhysicalChargingProcess, PhysicalChargingResult, PhysicalProcess, PhysicalResult,
)
from .surface_exchange import (
    SurfaceMaterialExchange, SurfaceProductPopulation, unresolved_surface_exchange,
    validate_surface_product_routing,
)
from .surface_product_redeposition_3d import (
    SurfaceProductRedeposition3DResult, SurfaceProductRedepositionContract3D,
    SurfaceProductRedepositionLaw3D, transport_surface_product_redeposition_3d,
)
from .physical_sputtering import (
    PhysicalSputterMechanism, PhysicalSputterParameters, PhysicalSputterState,
    PhysicalSputterStepResult,
)
from .chlorine_poly_si import (
    HwangGiapisClSiMechanism, HwangGiapisClSiParameters, HwangGiapisClSiState,
    HwangGiapisClSiStepResult, HwangGiapisClSiYield,
)
from .hwang_giapis_scatter_3d import (
    HwangGiapisForwardScatter3DResult, HwangGiapisSiO2ForwardScatter3D,
    NeutralSurfaceFlight3D, OutgoingNeutralParticleEvents3D,
    apply_hwang_giapis_forward_scatter_to_transport_3d,
    trace_neutral_surface_events_3d,
)
from .silicon_sf6o2 import (
    BelenSiliconParameters, BelenSiliconSF6O2Mechanism, BelenSiliconState,
    BelenSiliconStepResult,
)
from .material_mechanism_3d import (
    MaterialMechanismRouter3D, MaterialSurfaceState3D, MaterialSurfaceStepResult3D,
)
from .physical_arrivals_3d import (
    PhysicalArrivalSample3D, sample_physical_poisson_arrivals_3d,
)
from .profile_observables_3d import (
    EnsembleScalarEstimate3D, FeatureCenterline3D, FeatureCenterlineEnsemble3D,
    TrenchProfileEnsemble3D, TrenchProfileObservables3D,
    measure_feature_centerline_3d, measure_feature_centerline_ensemble_3d,
    measure_trench_profile_ensemble_3d, measure_trench_profile_observables_3d,
)
from .twist_campaign_3d import (
    TwistAspectRatioCampaign3DResult, TwistConditionCampaign3DResult,
    TwistEnsembleRefinement3DResult, TwistEnsembleRefinementContract3D,
    assess_twist_aspect_ratio_campaign_3d, assess_twist_ensemble_refinement_3d,
    score_twist_condition_campaign_3d,
)
from .notching_validation_3d import (
    NOTCHING_VALIDATION_PROTOCOL, NOZAWA_1995_NOTCH_CURVES_SHA256,
    NOZAWA_1995_PDF_SHA256,
    Nozawa1995NotchObservation3D,
    NotchingBenchmarkProtocol3D, NotchingBenchmarkScore3D,
    NotchingCalibrationReveal3D, NotchingHeldOutPrediction3D,
    load_nozawa_1995_notch_observations, score_notching_benchmark_3d,
)

# High-level 3D API (ViennaPS-shaped). Importing api pulls in the 3D engine (threed/warp).
from .api import Domain, SF6O2, Process, Result

__all__ = [
    # Explicit dimensional common engine and legacy compatibility interface.
    "COMMON_FEATURE_ENGINE", "COMMON_CHARGING_ENGINE", "COMMON_CHARGING_ENSEMBLE_ENGINE",
    "COMMON_FEATURE_MANIFEST_SCHEMA", "CHARGING_RUN_MANIFEST_SCHEMA",
    "PhysicalProcess", "PhysicalResult", "PhysicalChargingProcess", "PhysicalChargingResult",
    "PhysicalChargingEnsembleProcess", "PhysicalChargingEnsembleResult",
    "SurfaceMaterialExchange", "SurfaceProductPopulation", "unresolved_surface_exchange",
    "validate_surface_product_routing",
    "SurfaceProductRedeposition3DResult", "SurfaceProductRedepositionContract3D",
    "SurfaceProductRedepositionLaw3D", "transport_surface_product_redeposition_3d",
    "PhysicalSputterMechanism", "PhysicalSputterParameters", "PhysicalSputterState",
    "PhysicalSputterStepResult",
    "BelenSiliconParameters", "BelenSiliconSF6O2Mechanism", "BelenSiliconState",
    "BelenSiliconStepResult",
    "MaterialMechanismRouter3D", "MaterialSurfaceState3D",
    "MaterialSurfaceStepResult3D",
    "PhysicalArrivalSample3D", "sample_physical_poisson_arrivals_3d",
    "EnsembleScalarEstimate3D", "FeatureCenterline3D", "FeatureCenterlineEnsemble3D",
    "TrenchProfileEnsemble3D", "TrenchProfileObservables3D",
    "measure_feature_centerline_3d", "measure_feature_centerline_ensemble_3d",
    "measure_trench_profile_ensemble_3d", "measure_trench_profile_observables_3d",
    "TwistAspectRatioCampaign3DResult", "TwistConditionCampaign3DResult",
    "TwistEnsembleRefinement3DResult", "TwistEnsembleRefinementContract3D",
    "assess_twist_aspect_ratio_campaign_3d", "assess_twist_ensemble_refinement_3d",
    "score_twist_condition_campaign_3d",
    "NOTCHING_VALIDATION_PROTOCOL", "NOZAWA_1995_NOTCH_CURVES_SHA256",
    "NOZAWA_1995_PDF_SHA256",
    "Nozawa1995NotchObservation3D",
    "NotchingBenchmarkProtocol3D", "NotchingBenchmarkScore3D",
    "NotchingCalibrationReveal3D", "NotchingHeldOutPrediction3D",
    "load_nozawa_1995_notch_observations", "score_notching_benchmark_3d",
    "PlasmaDiagnosticState", "ReactorSpeciesFlux", "TabulatedReactorFluxDeck",
    "build_diagnostic_virtual_sheath_boundary", "build_tabulated_reactor_boundary",
    "load_krueger_2024_reactor_flux_deck",
    "Domain", "SF6O2", "Process", "Result",
    # config + low-level (full control)
    "PAR", "Flags", "DEFAULT_FLAGS", "run_etch",
    "make_trench", "extract_surface", "orient_normals", "seg_in_mask", "profile_bottom",
    "mc_flux", "_trace", "surface_rate", "surface_rate_langmuir",
    "advect", "extend_velocity", "reinit",
    "ours_profile", "depth_centre", "center_depth",
    "EnergeticFlux", "EnergeticYield", "FaceResolvedEnergeticFlux",
    "LowEnergyActivationYield", "ParameterEvidence",
    "SteinbruchelYield",
    "LaMagnaFluorocarbonParameters", "LaMagnaFluorocarbonState",
    "LaMagnaFluorocarbonStepResult", "LaMagnaGarozzoFluorocarbonMechanism",
    "ReducedSiO2FluorocarbonMechanism", "ReducedSiO2FluorocarbonParameters",
    "SiO2SurfaceState", "SurfaceFluxes",
    "BoundaryTransport3DResult", "ChargedSurfaceReimpactPopulation3D",
    "average_boundary_transport_results_3d",
    "merge_boundary_transport_results_3d", "trace_charged_surface_events_field_3d",
    "estimate_diffuse_form_factors_3d",
    "trace_boundary_state_field_3d",
    "trace_boundary_state_first_hit_3d",
    "NodalPoissonSystem3D", "PoissonDiagnostics3D", "assemble_q1_stiffness_3d",
    "lump_mixed_surface_density_3d", "lump_triangle_sheet_charge_3d",
    "ConductorTerminalCurrent3D", "RemotePadElectronCollector3D",
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
    "ChargedSurfaceCascade3DResult", "apply_charged_surface_response_to_transport_3d",
    "augment_transport_with_charged_reimpacts_3d", "derived_tail_bounce_budget_3d",
    "solve_charged_surface_cascade_3d",
    "SurfaceChargeRemap3DResult", "remap_surface_charge_3d",
    "ChargingCoevolution3DResult", "ChargingCoevolutionStep3DResult",
    "ExperimentalObservableTolerance3D", "PhysicalPatchBalance3D",
    "ResolvedBiasSegment3D", "SurfaceChargingSaturation3DResult",
    "SurfaceChargingSaturationError", "integrate_surface_charging_to_saturation_3d",
    "physical_surface_patch_groups_3d", "solve_charging_coevolution_3d",
    "PROFILE_STATIONARITY_CONTRACT_DRAFT", "ProfileChargingStationarity3DResult",
    "ProfileChargingStationarityBlock3D", "ProfileChargingStationarityContract3D",
    "assess_profile_charging_stationarity_3d",
    "CHARGING_CHECKPOINT_SCHEMA", "PhysicalChargingCheckpoint3D",
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
    "Jeon2022BoundaryClosure", "Jeong2023IonBoundaryClosure",
    "build_jeon_2022_boundary_state",
    "build_jeong_2023_boundary_state",
    "HWANG_GIAPIS_1997_EEDF_SHA256", "HWANG_GIAPIS_1997_IEDF_SHA256",
    "HWANG_GIAPIS_1997_PDF_SHA256",
    "build_hwang_giapis_1997_boundary_state",
    "HwangGiapisClSiMechanism", "HwangGiapisClSiParameters",
    "HwangGiapisClSiState", "HwangGiapisClSiStepResult", "HwangGiapisClSiYield",
    "HwangGiapisForwardScatter3DResult", "HwangGiapisSiO2ForwardScatter3D",
    "NeutralSurfaceFlight3D", "OutgoingNeutralParticleEvents3D",
    "apply_hwang_giapis_forward_scatter_to_transport_3d",
    "trace_neutral_surface_events_3d",
]

__version__ = "0.3.0"
