"""Profile-relevant stochastic stationarity contract for feature charging.

The signed CCA-R2 B1/B2 contract remains the default charging-convergence contract.  This module
defines an explicit, opt-in *draft* contract for a different question: whether another independent
physical-time charging block changes the field, kinetic currents, delivered surface fluxes, or the
profile increment by a physically relevant amount.  It never changes the kinetic operator and it
cannot authorize an experimental claim while the revision remains a draft.

Two consecutive blocks are required.  Their scoring sample epochs must be disjoint, and every block
must use the exact hard-visibility operator.  Standard errors enter as conservative upper bounds on
the measured changes; noisy agreement therefore cannot earn stationarity.
"""
from __future__ import annotations

from dataclasses import dataclass
from types import MappingProxyType
from typing import Mapping

import numpy as np


PROFILE_STATIONARITY_CONTRACT_DRAFT = "CCA-PROFILE-STATIONARY-2026-07-15-DRAFT"


def _readonly_vector(value, name, *, nonnegative=False):
    array = np.asarray(value, dtype=float).copy()
    if (array.ndim != 1 or array.size == 0 or np.any(~np.isfinite(array))
            or (nonnegative and np.any(array < 0.0))):
        qualifier = " nonnegative" if nonnegative else ""
        raise ValueError(f"{name} must be a nonempty finite{qualifier} vector")
    array.setflags(write=False)
    return array


def _readonly_array(value, name):
    array = np.asarray(value, dtype=float).copy()
    if array.size == 0 or np.any(~np.isfinite(array)):
        raise ValueError(f"{name} must be a nonempty finite array")
    array.setflags(write=False)
    return array


def _readonly_mapping(value, name, size, *, nonnegative=True):
    output = {}
    for species, supplied in dict(value).items():
        if not isinstance(species, str) or not species:
            raise ValueError(f"{name} species names must be nonempty strings")
        array = _readonly_vector(supplied, f"{name}[{species!r}]", nonnegative=nonnegative)
        if array.shape != (size,):
            raise ValueError(f"{name}[{species!r}] must match the face count")
        output[species] = array
    if not output:
        raise ValueError(f"{name} must contain at least one species")
    return MappingProxyType(output)


@dataclass(frozen=True)
class ProfileChargingStationarityContract3D:
    """Declared tolerances for two-block, profile-relevant charging stationarity.

    Absolute potential and profile-increment tolerances carry units.  Current and transported-flux
    tolerances are L1 changes normalized by their corresponding mean throughput.  Profile velocity
    uses a symmetric relative L1 change.  ``confidence_multiplier`` expands each stochastic metric
    by its independently estimated standard error before it is compared with the tolerance.
    """

    potential_drift_tolerance_v: float
    current_relative_l1_tolerance: float
    transported_flux_relative_l1_tolerance: float
    profile_velocity_relative_l1_tolerance: float
    profile_increment_tolerance_m: float
    minimum_independent_replicates: int = 2
    confidence_multiplier: float = 2.0
    revision: str = PROFILE_STATIONARITY_CONTRACT_DRAFT

    def __post_init__(self):
        values = np.asarray([
            self.potential_drift_tolerance_v,
            self.current_relative_l1_tolerance,
            self.transported_flux_relative_l1_tolerance,
            self.profile_velocity_relative_l1_tolerance,
            self.profile_increment_tolerance_m,
            self.confidence_multiplier], dtype=float)
        if (np.any(~np.isfinite(values)) or np.any(values < 0.0)
                or self.potential_drift_tolerance_v <= 0.0
                or self.current_relative_l1_tolerance <= 0.0
                or self.transported_flux_relative_l1_tolerance <= 0.0
                or self.profile_velocity_relative_l1_tolerance <= 0.0
                or self.profile_increment_tolerance_m <= 0.0
                or int(self.minimum_independent_replicates)
                != self.minimum_independent_replicates
                or self.minimum_independent_replicates < 2
                or not isinstance(self.revision, str) or not self.revision):
            raise ValueError("invalid profile-charging stationarity contract")
        object.__setattr__(
            self, "minimum_independent_replicates", int(self.minimum_independent_replicates))
        object.__setattr__(self, "confidence_multiplier", float(self.confidence_multiplier))

    @property
    def authorizes_experimental_claims(self):
        """Draft revisions are numerical acceptance paths, never validation signatures."""
        return False


@dataclass(frozen=True)
class ProfileChargingStationarityBlock3D:
    """One physical-time block plus independent exact-operator endpoint scoring.

    All current and flux fields are face-resolved.  Standard-error arrays describe the replicate
    mean at the block endpoint, not fluctuations reused from the transient that created the state.
    ``profile_increment_m`` is the signed normal displacement predicted from that independently
    averaged endpoint transport over the same declared profile step.
    """

    potential_start_v: np.ndarray
    potential_end_v: np.ndarray
    positive_face_current_density_a_m2: np.ndarray
    negative_face_current_density_a_m2: np.ndarray
    net_face_current_standard_error_a_m2: np.ndarray
    face_area_m2: np.ndarray
    species_face_flux_m2_s: Mapping[str, np.ndarray]
    species_face_flux_standard_error_m2_s: Mapping[str, np.ndarray]
    profile_velocity_m_s: np.ndarray
    profile_velocity_standard_error_m_s: np.ndarray
    profile_increment_m: np.ndarray
    independent_replicates: int
    scoring_sampling_epochs: tuple[int, ...]
    duration_s: float
    exact_hard_visibility: bool = True

    def __post_init__(self):
        area = _readonly_vector(self.face_area_m2, "face_area_m2", nonnegative=True)
        if np.any(area <= 0.0):
            raise ValueError("face_area_m2 must be strictly positive")
        size = len(area)
        potential_start = _readonly_array(self.potential_start_v, "potential_start_v")
        potential_end = _readonly_array(self.potential_end_v, "potential_end_v")
        if potential_start.shape != potential_end.shape:
            raise ValueError("block endpoint potentials must share one Poisson-grid shape")
        vectors = {}
        for name, nonnegative in (
                ("positive_face_current_density_a_m2", True),
                ("negative_face_current_density_a_m2", True),
                ("net_face_current_standard_error_a_m2", True),
                ("profile_velocity_m_s", False),
                ("profile_velocity_standard_error_m_s", True),
                ("profile_increment_m", False)):
            vectors[name] = _readonly_vector(getattr(self, name), name, nonnegative=nonnegative)
            if vectors[name].shape != (size,):
                raise ValueError(f"{name} must match the face count")
        flux = _readonly_mapping(
            self.species_face_flux_m2_s, "species_face_flux_m2_s", size)
        flux_error = _readonly_mapping(
            self.species_face_flux_standard_error_m2_s,
            "species_face_flux_standard_error_m2_s", size)
        epochs = tuple(int(value) for value in self.scoring_sampling_epochs)
        if (set(flux) != set(flux_error)
                or int(self.independent_replicates) != self.independent_replicates
                or self.independent_replicates < 2
                or len(epochs) != self.independent_replicates
                or len(set(epochs)) != len(epochs) or any(value < 0 for value in epochs)
                or not np.isfinite(self.duration_s) or self.duration_s <= 0.0
                or not isinstance(self.exact_hard_visibility, (bool, np.bool_))):
            raise ValueError("invalid independent stationarity-block scoring contract")
        object.__setattr__(self, "face_area_m2", area)
        object.__setattr__(self, "potential_start_v", potential_start)
        object.__setattr__(self, "potential_end_v", potential_end)
        for name, value in vectors.items():
            object.__setattr__(self, name, value)
        object.__setattr__(self, "species_face_flux_m2_s", flux)
        object.__setattr__(self, "species_face_flux_standard_error_m2_s", flux_error)
        object.__setattr__(self, "independent_replicates", int(self.independent_replicates))
        object.__setattr__(self, "scoring_sampling_epochs", epochs)
        object.__setattr__(self, "duration_s", float(self.duration_s))
        object.__setattr__(self, "exact_hard_visibility", bool(self.exact_hard_visibility))


@dataclass(frozen=True)
class ProfileChargingStationarity3DResult:
    """Gate result for two consecutive independently scored charging blocks."""

    passed: bool
    potential_drift_upper_v: float
    current_relative_l1_upper: float
    transported_flux_relative_l1_upper: float
    profile_velocity_relative_l1_upper: float
    profile_increment_difference_upper_m: float
    reasons: tuple[str, ...]
    diagnostics: Mapping[str, object]
    contract_revision: str

    def __post_init__(self):
        values = np.asarray([
            self.potential_drift_upper_v, self.current_relative_l1_upper,
            self.transported_flux_relative_l1_upper,
            self.profile_velocity_relative_l1_upper,
            self.profile_increment_difference_upper_m], dtype=float)
        if (np.any(np.isnan(values)) or np.any(values < 0.0)
                or not isinstance(self.contract_revision, str)
                or not self.contract_revision):
            raise ValueError("invalid profile-charging stationarity result")
        object.__setattr__(self, "passed", bool(self.passed))
        object.__setattr__(self, "reasons", tuple(self.reasons))
        object.__setattr__(self, "diagnostics", MappingProxyType(dict(self.diagnostics)))


def _relative_l1_upper(first, second, first_error, second_error, weight, confidence):
    first = np.asarray(first, dtype=float)
    second = np.asarray(second, dtype=float)
    uncertainty = confidence * np.sqrt(
        np.asarray(first_error, dtype=float) ** 2
        + np.asarray(second_error, dtype=float) ** 2)
    numerator = float(np.sum((np.abs(second - first) + uncertainty) * weight))
    denominator = float(np.sum(0.5 * (np.abs(first) + np.abs(second)) * weight))
    if denominator <= np.finfo(float).tiny:
        return 0.0 if numerator <= np.finfo(float).tiny else float("inf")
    return numerator / denominator


def assess_profile_charging_stationarity_3d(
        first: ProfileChargingStationarityBlock3D,
        second: ProfileChargingStationarityBlock3D,
        contract: ProfileChargingStationarityContract3D):
    """Assess two consecutive physical-time blocks without changing their operator.

    The second block's potential drift is the field-stationarity test.  Current changes are
    normalized by charged-particle throughput, not the small ion current of one shadowed patch.
    Delivered flux and profile-motion tests ensure that any residual charging motion is irrelevant
    to the downstream feature evolution being requested.
    """
    if (not isinstance(first, ProfileChargingStationarityBlock3D)
            or not isinstance(second, ProfileChargingStationarityBlock3D)
            or not isinstance(contract, ProfileChargingStationarityContract3D)):
        raise TypeError("stationarity assessment requires two blocks and one contract")
    size = len(first.face_area_m2)
    if (len(second.face_area_m2) != size
            or first.potential_end_v.shape != second.potential_start_v.shape
            or set(first.species_face_flux_m2_s) != set(second.species_face_flux_m2_s)
            or not np.allclose(first.face_area_m2, second.face_area_m2, rtol=2e-13, atol=0.0)):
        raise ValueError("stationarity blocks must score the same fixed surface mesh")
    if np.max(np.abs(first.potential_end_v - second.potential_start_v)) > 1e-10 * max(
            1.0, float(np.max(np.abs(first.potential_end_v)))):
        raise ValueError("stationarity blocks must be consecutive in physical state")
    if set(first.scoring_sampling_epochs) & set(second.scoring_sampling_epochs):
        raise ValueError("stationarity blocks must use disjoint endpoint-scoring sample epochs")
    if (first.independent_replicates < contract.minimum_independent_replicates
            or second.independent_replicates < contract.minimum_independent_replicates):
        raise ValueError("stationarity blocks do not meet the replicate-count contract")

    confidence = contract.confidence_multiplier
    potential_drift = float(np.max(np.abs(
        second.potential_end_v - second.potential_start_v)))

    first_net = (first.positive_face_current_density_a_m2
                 - first.negative_face_current_density_a_m2)
    second_net = (second.positive_face_current_density_a_m2
                  - second.negative_face_current_density_a_m2)
    current_uncertainty = confidence * np.sqrt(
        first.net_face_current_standard_error_a_m2 ** 2
        + second.net_face_current_standard_error_a_m2 ** 2)
    current_numerator = float(np.sum(
        (np.abs(second_net - first_net) + current_uncertainty) * first.face_area_m2))
    current_denominator = float(np.sum(0.25 * (
        first.positive_face_current_density_a_m2
        + first.negative_face_current_density_a_m2
        + second.positive_face_current_density_a_m2
        + second.negative_face_current_density_a_m2) * first.face_area_m2))
    current_relative = (
        current_numerator / current_denominator
        if current_denominator > np.finfo(float).tiny
        else (0.0 if current_numerator <= np.finfo(float).tiny else float("inf")))

    flux_numerator = 0.0
    flux_denominator = 0.0
    for species in sorted(first.species_face_flux_m2_s):
        first_flux = first.species_face_flux_m2_s[species]
        second_flux = second.species_face_flux_m2_s[species]
        uncertainty = confidence * np.sqrt(
            first.species_face_flux_standard_error_m2_s[species] ** 2
            + second.species_face_flux_standard_error_m2_s[species] ** 2)
        flux_numerator += float(np.sum(
            (np.abs(second_flux - first_flux) + uncertainty) * first.face_area_m2))
        flux_denominator += float(np.sum(
            0.5 * (first_flux + second_flux) * first.face_area_m2))
    flux_relative = (
        flux_numerator / flux_denominator
        if flux_denominator > np.finfo(float).tiny
        else (0.0 if flux_numerator <= np.finfo(float).tiny else float("inf")))

    velocity_relative = _relative_l1_upper(
        first.profile_velocity_m_s, second.profile_velocity_m_s,
        first.profile_velocity_standard_error_m_s,
        second.profile_velocity_standard_error_m_s,
        first.face_area_m2, confidence)
    increment_uncertainty = confidence * np.sqrt(
        (first.profile_velocity_standard_error_m_s * first.duration_s) ** 2
        + (second.profile_velocity_standard_error_m_s * second.duration_s) ** 2)
    increment_difference = float(np.max(
        np.abs(second.profile_increment_m - first.profile_increment_m)
        + increment_uncertainty))

    reasons = []
    if not first.exact_hard_visibility or not second.exact_hard_visibility:
        reasons.append("exact hard-visibility operator was not used in both blocks")
    if potential_drift > contract.potential_drift_tolerance_v:
        reasons.append("second-block potential drift exceeds tolerance")
    if current_relative > contract.current_relative_l1_tolerance:
        reasons.append("independent kinetic-current change exceeds tolerance")
    if flux_relative > contract.transported_flux_relative_l1_tolerance:
        reasons.append("delivered species-flux change exceeds tolerance")
    if velocity_relative > contract.profile_velocity_relative_l1_tolerance:
        reasons.append("predicted profile-velocity change exceeds tolerance")
    if increment_difference > contract.profile_increment_tolerance_m:
        reasons.append("predicted profile-increment change exceeds tolerance")
    diagnostics = dict(
        contract_kind="profile-relevant two-block stochastic stationarity",
        experimental_claim_authorized=contract.authorizes_experimental_claims,
        exact_operator_statement="independent scoring; exact hard visibility; no residual smoothing",
        first_scoring_sampling_epochs=first.scoring_sampling_epochs,
        second_scoring_sampling_epochs=second.scoring_sampling_epochs,
        first_independent_replicates=first.independent_replicates,
        second_independent_replicates=second.independent_replicates,
        confidence_multiplier=confidence,
        potential_drift_tolerance_v=contract.potential_drift_tolerance_v,
        current_relative_l1_tolerance=contract.current_relative_l1_tolerance,
        transported_flux_relative_l1_tolerance=(
            contract.transported_flux_relative_l1_tolerance),
        profile_velocity_relative_l1_tolerance=(
            contract.profile_velocity_relative_l1_tolerance),
        profile_increment_tolerance_m=contract.profile_increment_tolerance_m)
    return ProfileChargingStationarity3DResult(
        passed=not reasons,
        potential_drift_upper_v=potential_drift,
        current_relative_l1_upper=current_relative,
        transported_flux_relative_l1_upper=flux_relative,
        profile_velocity_relative_l1_upper=velocity_relative,
        profile_increment_difference_upper_m=increment_difference,
        reasons=tuple(reasons), diagnostics=diagnostics,
        contract_revision=contract.revision)
