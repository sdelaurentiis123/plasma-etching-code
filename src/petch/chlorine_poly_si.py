"""Source-backed Cl+/poly-Si removal law for the Nozawa notching replay.

The model is the deliberately small profile law used by Hwang and Giapis,
JVST B 15, 70 (1997), Eq. (4.1), DOI 10.1116/1.589258.  It is not a
universal chlorine chemistry: it represents energetic Cl+ (and explicitly
identified fast neutral Cl) striking a chlorinated poly-Si surface during
overetch.  Spontaneous neutral-Cl etching, evolving chlorination, reaction
products, implantation, and mask/oxide removal are outside this law and are
reported rather than hidden in a fitted rate.
"""
from __future__ import annotations

from dataclasses import dataclass
from types import MappingProxyType
from typing import Mapping

import numpy as np

from .surface_exchange import SurfaceMaterialExchange, unresolved_surface_exchange
from .surface_kinetics import MechanismValidity, ParameterEvidence, SurfaceFluxes


@dataclass(frozen=True)
class HwangGiapisClSiYield:
    """Hwang--Giapis Eq. (4.1), including its critical-angle cap.

    ``Y = C max(sqrt(E)-sqrt(E_th), 0)`` below the critical angle.  At
    larger incidence angles the paper multiplies this value by
    ``cos(theta) / cos(theta_critical)``, so the yield decreases toward
    grazing incidence instead of growing there.  The paper chose
    ``C=0.1 / sqrt(eV)``, ``E_th=10 eV``, and
    ``theta_critical=45 deg``.  The result is a yield and is not clipped to
    one.
    """

    prefactor_per_sqrt_eV: float = 0.1
    threshold_energy_eV: float = 10.0
    critical_angle_deg: float = 45.0

    def __post_init__(self):
        if (not np.isfinite(self.prefactor_per_sqrt_eV)
                or self.prefactor_per_sqrt_eV < 0.0
                or not np.isfinite(self.threshold_energy_eV)
                or self.threshold_energy_eV < 0.0
                or not np.isfinite(self.critical_angle_deg)
                or not 0.0 < self.critical_angle_deg < 90.0):
            raise ValueError("invalid Hwang--Giapis Cl/poly-Si yield parameters")

    def evaluate(self, energy_eV, cosine_incidence):
        energy = np.asarray(energy_eV, dtype=float)
        cosine = np.asarray(cosine_incidence, dtype=float)
        if (np.any(~np.isfinite(energy)) or np.any(energy < 0.0)
                or np.any(~np.isfinite(cosine))
                or np.any((cosine < 0.0) | (cosine > 1.0))):
            raise ValueError("incident energies/cosines are outside the physical domain")
        energy_factor = np.maximum(
            np.sqrt(energy) - np.sqrt(self.threshold_energy_eV), 0.0)
        critical_cosine = np.cos(np.deg2rad(self.critical_angle_deg))
        angular_factor = np.minimum(cosine / critical_cosine, 1.0)
        return self.prefactor_per_sqrt_eV * energy_factor * angular_factor


@dataclass(frozen=True)
class HwangGiapisClSiState:
    removed_si_atoms_m2: np.ndarray | float = 0.0

    def __post_init__(self):
        value = np.asarray(self.removed_si_atoms_m2, dtype=float).copy()
        if np.any(~np.isfinite(value)) or np.any(value < 0.0):
            raise ValueError("removed silicon inventory must be finite and nonnegative")
        value.setflags(write=False)
        object.__setattr__(self, "removed_si_atoms_m2", value)

    @classmethod
    def bare(cls, shape=()):
        return cls(np.zeros(shape))

    def conservative_surface_fields(self):
        return {"removed_si_atoms_m2": self.removed_si_atoms_m2}

    def conservative_surface_upper_bounds(self):
        return {"removed_si_atoms_m2": None}

    def with_conservative_surface_fields(self, fields):
        fields = dict(fields)
        if set(fields) != {"removed_si_atoms_m2"}:
            raise ValueError("Cl/poly-Si remap fields do not match its state contract")
        return type(self)(fields["removed_si_atoms_m2"])


@dataclass(frozen=True)
class HwangGiapisClSiParameters:
    bulk_si_atom_density_m3: float
    projectile_species: tuple[str, ...]
    yield_law: HwangGiapisClSiYield
    evidence: Mapping[str, ParameterEvidence]
    known_omissions: tuple[str, ...] = (
        "spontaneous neutral-Cl etching is omitted; Hwang--Giapis estimated up to 60 nm over the full process",
        "surface chlorination and its history are not evolved",
        "removed-Si reaction-product identity and redeposition are unresolved",
        "photoresist and SiO2 removal are excluded by the benchmark model",
    )

    def __post_init__(self):
        projectiles = tuple(self.projectile_species)
        evidence = dict(self.evidence)
        required = {
            "bulk_si_atom_density_m3", "yield_prefactor_per_sqrt_eV",
            "threshold_energy_eV", "critical_angle_deg",
        }
        if (not np.isfinite(self.bulk_si_atom_density_m3)
                or self.bulk_si_atom_density_m3 <= 0.0
                or not projectiles or any(not name for name in projectiles)
                or len(set(projectiles)) != len(projectiles)
                or not isinstance(self.yield_law, HwangGiapisClSiYield)
                or set(evidence) != required
                or any(not isinstance(item, ParameterEvidence)
                       for item in evidence.values())):
            raise ValueError("invalid Hwang--Giapis Cl/poly-Si parameters or evidence")
        object.__setattr__(self, "projectile_species", projectiles)
        object.__setattr__(self, "evidence", MappingProxyType(evidence))
        object.__setattr__(self, "known_omissions", tuple(self.known_omissions))

    @classmethod
    def hwang_giapis_1997(cls):
        paper = (
            "Hwang & Giapis, JVST B 15, 70 (1997), "
            "DOI 10.1116/1.589258, Sec. IV A, Eq. (4.1)")
        return cls(
            bulk_si_atom_density_m3=5.0e28,
            projectile_species=("Cl+", "Cl_fast_neutral"),
            yield_law=HwangGiapisClSiYield(),
            evidence={
                "bulk_si_atom_density_m3": ParameterEvidence(
                    "derived from crystalline-Si density and molar mass; petch ALE inventory uses 5.0e28 m^-3",
                    "derived physical constant", supports_prediction_within_declared_domain=True),
                "yield_prefactor_per_sqrt_eV": ParameterEvidence(
                    paper, "published computational time-scale choice", note=(
                        "C=0.1 makes a normally incident 177 eV ion approximately unit yield; "
                        "it does not independently validate absolute etch time")),
                "threshold_energy_eV": ParameterEvidence(
                    paper, "published externally measured input", note="E_th=10 eV",
                    supports_prediction_within_declared_domain=True),
                "critical_angle_deg": ParameterEvidence(
                    paper, "published model assumption", note=(
                        "45 deg central value; the paper brackets 30--60 deg and reports weak notch sensitivity")),
            })


@dataclass(frozen=True)
class HwangGiapisClSiStepResult:
    state: HwangGiapisClSiState
    etch_velocity_m_s: np.ndarray
    removed_si_atoms_m2: np.ndarray
    material_exchange: SurfaceMaterialExchange
    product_populations: tuple = ()
    validity: MechanismValidity | None = None

    def __post_init__(self):
        for name in ("etch_velocity_m_s", "removed_si_atoms_m2"):
            value = np.asarray(getattr(self, name), dtype=float).copy()
            if np.any(~np.isfinite(value)) or np.any(value < 0.0):
                raise ValueError("invalid Cl/poly-Si step result")
            value.setflags(write=False)
            object.__setattr__(self, name, value)
        object.__setattr__(self, "product_populations", tuple(self.product_populations))
        if (not isinstance(self.state, HwangGiapisClSiState)
                or not isinstance(self.material_exchange, SurfaceMaterialExchange)
                or not isinstance(self.validity, MechanismValidity)
                or self.product_populations):
            raise ValueError("invalid Cl/poly-Si result contract")


class HwangGiapisClSiMechanism:
    """Conservative energetic chlorine removal on chlorinated poly-Si."""

    def __init__(self, parameters: HwangGiapisClSiParameters | None = None):
        self.parameters = (
            HwangGiapisClSiParameters.hwang_giapis_1997()
            if parameters is None else parameters)
        if not isinstance(self.parameters, HwangGiapisClSiParameters):
            raise TypeError("parameters must be HwangGiapisClSiParameters")

    @property
    def provenance(self):
        law = self.parameters.yield_law
        return MappingProxyType({
            "model": "Hwang--Giapis Cl+/chlorinated-poly-Si Eq. (4.1)",
            "projectile_species": self.parameters.projectile_species,
            "parameters": {
                "bulk_si_atom_density_m3": self.parameters.bulk_si_atom_density_m3,
                "yield_prefactor_per_sqrt_eV": law.prefactor_per_sqrt_eV,
                "threshold_energy_eV": law.threshold_energy_eV,
                "critical_angle_deg": law.critical_angle_deg,
            },
            "bounds": {
                "bulk_si_atom_density_m3": (4.9e28, 5.1e28),
                "yield_prefactor_per_sqrt_eV": (0.1, 0.1),
                "threshold_energy_eV": (10.0, 20.0),
                "critical_angle_deg": (30.0, 60.0),
            },
            "evidence": {
                name: item.source for name, item in self.parameters.evidence.items()},
            "known_omissions": self.parameters.known_omissions,
            "claim": (
                "published benchmark model; C fixes the computational etch-time scale and "
                "does not independently validate absolute rate"),
        })

    @staticmethod
    def initial_state(shape=()):
        return HwangGiapisClSiState.bare(shape)

    def validity(self, fluxes: SurfaceFluxes):
        unsupported_neutral = tuple(sorted(
            name for name, value in fluxes.neutral_flux_m2_s.items()
            if np.any(np.asarray(value) > 0.0)))
        unsupported_energetic = tuple(sorted({
            population.name for population in fluxes.energetic_fluxes
            if population.name not in self.parameters.projectile_species
            and np.any(np.asarray(population.flux_m2_s) > 0.0)}))
        reasons = () if not (unsupported_neutral or unsupported_energetic) else (
            "positive incident flux has no declared Cl/poly-Si reaction channel",)
        nonpredictive = tuple(sorted(
            name for name, item in self.parameters.evidence.items()
            if not item.supports_prediction_within_declared_domain))
        return MechanismValidity(
            within_declared_scope=not reasons, reasons=reasons,
            unsupported_neutral_species=unsupported_neutral,
            known_model_form_omissions=self.parameters.known_omissions,
            parameter_evidence_supports_prediction=not nonpredictive,
            nonpredictive_parameters=nonpredictive)

    def advance(self, state, fluxes: SurfaceFluxes, duration_s: float, *, strict=True):
        if not isinstance(state, HwangGiapisClSiState):
            raise TypeError("Cl/poly-Si mechanism requires HwangGiapisClSiState")
        if not np.isfinite(duration_s) or duration_s < 0.0:
            raise ValueError("duration_s must be finite and nonnegative")
        validity = self.validity(fluxes)
        if strict and not validity.within_declared_scope:
            raise ValueError("surface mechanism outside declared scope: " + "; ".join(validity.reasons))
        shape = state.removed_si_atoms_m2.shape
        removal_rate = np.zeros(shape)
        for population in fluxes.energetic_fluxes:
            if population.name in self.parameters.projectile_species:
                removal_rate = removal_rate + np.broadcast_to(
                    population.yield_rate_m2_s(self.parameters.yield_law), shape)
        removed = removal_rate * float(duration_s)
        exchange = unresolved_surface_exchange(
            removed_units_m2={"poly_si_atoms": removed}, limitations=(
                "removed poly-Si product identity is not declared by the benchmark law",
                "unresolved products are not transported or redeposited",
            ))
        return HwangGiapisClSiStepResult(
            HwangGiapisClSiState(state.removed_si_atoms_m2 + removed),
            removal_rate / self.parameters.bulk_si_atom_density_m3,
            removed, exchange, (), validity)
