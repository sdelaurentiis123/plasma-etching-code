"""Fail-closed TiO2/fluorocarbon surface-deck contract.

The common fluorinated-oxide kernel is executable, conservative, and
material-labelled, but published SiO2 coefficients do not become TiO2 data by
changing that label.  This module makes every TiO2-specific numeric input and
its evidence explicit.  Choi's nonmonotonic oxygen response is represented by
a separate oxygen-blocked TiO2 site fraction with analytic adsorption and
energetic-cleanup updates; its two coefficients remain mandatory data inputs.

Consequently a complete instance can be built for controlled sensitivity
studies, while target prediction remains fail-closed until the remaining
roughness model-form gap and parameter evidence are resolved. There are
deliberately no default surface probabilities or yield curves in this module.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from types import MappingProxyType
from typing import Mapping

from .surface_kinetics import (
    EnergeticYield,
    ParameterEvidence,
    ReducedFluorinatedOxideMechanism,
    ReducedFluorinatedOxideParameters,
)
from .tio2_ion_dose import tio2_formula_unit_density_m3


TIO2_REDUCED_SURFACE_REQUIRED_EVIDENCE = (
    "site_density_m2",
    "bulk_formula_density_m3",
    "polymer_monolayer_density_m2",
    "polymer_bulk_unit_density_m3",
    "complex_formation_probability",
    "polymer_deposition_probability_on_substrate",
    "polymer_deposition_probability_on_polymer",
    "oxygen_polymer_etch_probability",
    "oxygen_blocking_probability",
    "oxygen_blocker_ion_removal_yield",
    "bare_sio2_yield",
    "complex_sio2_yield",
    "polymer_sputter_yield",
)

TIO2_TARGET_MODEL_FORM_GAPS = (
    "chemistry_dependent_roughness_evolution",
)


@dataclass(frozen=True)
class Tio2SurfaceDeckReadiness:
    """Evidence and model-form status for one numerical TiO2 deck."""

    missing_parameter_evidence: tuple[str, ...]
    nonpredictive_parameter_evidence: tuple[str, ...]
    unresolved_model_form: tuple[str, ...]
    supports_reduced_sensitivity: bool
    supports_absolute_target_prediction: bool


@dataclass(frozen=True)
class Tio2ReducedSurfaceDeck:
    """All numerical inputs needed by the current reduced TiO2 kernel.

    The historical ``bare_sio2_yield`` evidence key is retained because it is
    part of the common-kernel validity schema.  In this deck it always means
    the bare-TiO2 energetic removal law; material labels and provenance make
    that interpretation explicit.
    """

    mass_density_kg_m3: float
    site_density_m2: float
    passivation_monolayer_density_m2: float
    passivation_bulk_unit_density_m3: float
    fluorination_probability: Mapping[str, float]
    passivation_deposition_probability_on_tio2: Mapping[str, float]
    passivation_deposition_probability_on_passivation: Mapping[str, float]
    oxygen_species: str
    oxygen_passivation_removal_probability: float
    oxygen_blocking_probability: float
    oxygen_blocker_ion_removal_yield: EnergeticYield
    bare_tio2_yield: EnergeticYield
    fluorinated_tio2_yield: EnergeticYield
    passivation_sputter_yield: EnergeticYield
    evidence: Mapping[str, ParameterEvidence]
    declared_inert_neutral_species: tuple[str, ...] = ()
    additional_known_omissions: tuple[str, ...] = ()

    def __post_init__(self):
        evidence = dict(self.evidence)
        if any(not isinstance(item, ParameterEvidence) for item in evidence.values()):
            raise TypeError("TiO2 deck evidence values must be ParameterEvidence objects")
        object.__setattr__(self, "evidence", MappingProxyType(evidence))
        object.__setattr__(
            self, "declared_inert_neutral_species",
            tuple(self.declared_inert_neutral_species),
        )
        object.__setattr__(
            self, "additional_known_omissions",
            tuple(self.additional_known_omissions),
        )

    def readiness(self) -> Tio2SurfaceDeckReadiness:
        required = set(TIO2_REDUCED_SURFACE_REQUIRED_EVIDENCE)
        missing = tuple(sorted(required - set(self.evidence)))
        nonpredictive = tuple(sorted(
            name for name in required
            if name not in self.evidence
            or not self.evidence[name].supports_prediction_within_declared_domain
        ))
        model_form = tuple(TIO2_TARGET_MODEL_FORM_GAPS)
        reduced_ready = not missing
        return Tio2SurfaceDeckReadiness(
            missing_parameter_evidence=missing,
            nonpredictive_parameter_evidence=nonpredictive,
            unresolved_model_form=model_form,
            supports_reduced_sensitivity=reduced_ready,
            supports_absolute_target_prediction=(
                reduced_ready and not nonpredictive and not model_form
            ),
        )

    def build_parameters(
        self, *, allow_reduced_sensitivity: bool = False,
    ) -> ReducedFluorinatedOxideParameters:
        """Build the common-kernel parameters without hiding claim limits."""
        status = self.readiness()
        if status.missing_parameter_evidence:
            raise ValueError(
                "TiO2 deck is missing parameter evidence: "
                + ", ".join(status.missing_parameter_evidence)
            )
        if not allow_reduced_sensitivity and status.unresolved_model_form:
            raise ValueError(
                "TiO2 target prediction is blocked by unresolved model form: "
                + ", ".join(status.unresolved_model_form)
            )
        omissions = tuple(dict.fromkeys(
            TIO2_TARGET_MODEL_FORM_GAPS + self.additional_known_omissions
        ))
        return ReducedFluorinatedOxideParameters(
            site_density_m2=self.site_density_m2,
            bulk_formula_density_m3=tio2_formula_unit_density_m3(
                self.mass_density_kg_m3
            ),
            polymer_monolayer_density_m2=self.passivation_monolayer_density_m2,
            polymer_bulk_unit_density_m3=self.passivation_bulk_unit_density_m3,
            complex_formation_probability=dict(self.fluorination_probability),
            polymer_deposition_probability_on_substrate=dict(
                self.passivation_deposition_probability_on_tio2
            ),
            polymer_deposition_probability_on_polymer=dict(
                self.passivation_deposition_probability_on_passivation
            ),
            oxygen_species=self.oxygen_species,
            oxygen_polymer_etch_probability=(
                self.oxygen_passivation_removal_probability
            ),
            oxygen_blocking_probability=self.oxygen_blocking_probability,
            oxygen_blocker_ion_removal_yield=(
                self.oxygen_blocker_ion_removal_yield
            ),
            bare_sio2_yield=self.bare_tio2_yield,
            complex_sio2_yield=self.fluorinated_tio2_yield,
            polymer_sputter_yield=self.passivation_sputter_yield,
            material_name="ALD TiO2",
            material_inventory_name="TiO2_formula_unit",
            declared_inert_neutral_species=self.declared_inert_neutral_species,
            evidence=self.evidence,
            known_omissions=omissions,
        )

    def build_mechanism(
        self, *, allow_reduced_sensitivity: bool = False,
    ) -> ReducedFluorinatedOxideMechanism:
        return ReducedFluorinatedOxideMechanism(self.build_parameters(
            allow_reduced_sensitivity=allow_reduced_sensitivity
        ))


__all__ = [
    "TIO2_REDUCED_SURFACE_REQUIRED_EVIDENCE",
    "TIO2_TARGET_MODEL_FORM_GAPS",
    "Tio2ReducedSurfaceDeck",
    "Tio2SurfaceDeckReadiness",
]
