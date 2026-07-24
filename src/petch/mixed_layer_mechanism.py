"""Feature-engine adapter for the element-resolved mixed-layer chemistry.

Wraps `petch.mixed_layer` behind the same interface the reduced SiO2
fluorocarbon mechanism exposes to the 3-D feature engine: `initial_state`,
`advance(state, fluxes, duration_s)` returning a `SurfaceStepResult`,
`neutral_reaction_probability`, `validity`, and the conservative surface-state
remap contract. Declared omissions are reported, never silently absorbed:

- Ion spectra are compressed to a flux-weighted mean energy and cosine per
  face before the stopping-table lookups (Jensen bias; exact per-event
  integration is the declared follow-up).
- Product routing (SiF4 / CO / COF2 emission laws) is reported as an
  unresolved ledger, exactly like the incumbent reduced mechanism.
- No spontaneous (ion-free) chemical etch channel yet.
"""

from __future__ import annotations

from dataclasses import dataclass
from types import MappingProxyType

import numpy as np

from petch.mixed_layer import (
    _MONOLAYER_AREAL_M2,
    MixedLayerParams,
    MixedLayerState,
    SurfaceFluxes as ModuleFluxes,
    step as mixed_layer_step,
)
from petch.surface_exchange import unresolved_surface_exchange
from petch.surface_kinetics import (
    EnergeticFlux,
    FaceResolvedEnergeticFlux,
    MechanismValidity,
    SurfaceFluxes,
    SurfaceStepResult,
)

_SIO2_FORMULA_DENSITY_M3 = 2.2e28
_CARBON_ATOM_DENSITY_M3 = 1.0e29

_STATE_FIELDS = ("n_c_film", "n_f_film", "n_si", "n_o", "n_c", "n_f",
                 "removed_formula_units_m2")


@dataclass(frozen=True)
class MixedLayerSurfaceState:
    """Per-face mixed-layer state under the engine's conservative remap contract."""

    n_c_film: np.ndarray | float
    n_f_film: np.ndarray | float
    n_si: np.ndarray | float
    n_o: np.ndarray | float
    n_c: np.ndarray | float
    n_f: np.ndarray | float
    removed_formula_units_m2: np.ndarray | float = 0.0

    def __post_init__(self):
        arrays = np.broadcast_arrays(*[
            np.asarray(getattr(self, name), dtype=float)
            for name in _STATE_FIELDS])
        for name, array in zip(_STATE_FIELDS, arrays):
            array = np.array(array, copy=True)
            if np.any(~np.isfinite(array)) or np.any(array < 0.0):
                raise ValueError(f"invalid mixed-layer state field {name}")
            array.setflags(write=False)
            object.__setattr__(self, name, array)

    @classmethod
    def bare(cls, shape=()):
        zero = np.zeros(shape)
        return cls(zero, zero, zero, zero, zero, zero, zero)

    def conservative_surface_fields(self):
        return {name: getattr(self, name) for name in _STATE_FIELDS}

    def conservative_surface_upper_bounds(self):
        return {name: None for name in _STATE_FIELDS}

    def with_conservative_surface_fields(self, fields):
        fields = dict(fields)
        if set(fields) != set(_STATE_FIELDS):
            raise ValueError("mixed-layer remap fields do not match its state contract")
        return type(self)(*[fields[name] for name in _STATE_FIELDS])

    def _module_state(self):
        return MixedLayerState(self.n_c_film, self.n_f_film, self.n_si,
                               self.n_o, self.n_c, self.n_f)


class MixedLayerMechanism:
    """Mixed-layer chemistry behind the reduced-mechanism engine interface."""

    def __init__(self, parameters: MixedLayerParams, *,
                 precursor_species=("CFx",), fluorine_species=("F",),
                 oxygen_species=("O",), inert_species=(),
                 chemisorption_probability=None,
                 default_max_step_s=0.01):
        if not isinstance(parameters, MixedLayerParams):
            raise TypeError("parameters must be MixedLayerParams")
        self.parameters = parameters
        # precursor_species: iterable of names (1 C, params ratio F each) or a
        # mapping {name: (carbon_atoms, fluorine_atoms)} for stoichiometric
        # multi-species chemistries (e.g. Krueger's CF/CF2/CF3/C2F3/...).
        if isinstance(precursor_species, dict):
            self.precursor_stoichiometry = {
                str(name): (float(c), float(f))
                for name, (c, f) in precursor_species.items()}
        else:
            self.precursor_stoichiometry = {
                str(name): (1.0, float(parameters.precursor_fc_ratio))
                for name in precursor_species}
        self.precursor_species = tuple(self.precursor_stoichiometry)
        self.chemisorption_probability = {
            str(name): float(value)
            for name, value in dict(chemisorption_probability or {}).items()}
        unknown = set(self.chemisorption_probability) - set(self.precursor_species)
        if unknown:
            raise ValueError(
                f"chemisorption probabilities for unmapped species: {sorted(unknown)}")
        self.fluorine_species = tuple(fluorine_species)
        self.oxygen_species = tuple(oxygen_species)
        self.inert_species = tuple(inert_species)
        self.default_max_step_s = float(default_max_step_s)
        overlap = (set(self.precursor_species) & set(self.fluorine_species)
                   | set(self.precursor_species) & set(self.oxygen_species)
                   | set(self.fluorine_species) & set(self.oxygen_species))
        if overlap:
            raise ValueError(f"species mapped to multiple channels: {sorted(overlap)}")
        self.provenance = MappingProxyType({
            "model": "mixed-layer-two-reservoir-v1",
            "substrate": parameters.substrate,
            "parameters": {
                "sticking_probability": float(parameters.sticking_probability),
                "fluorine_film_sticking": float(parameters.fluorine_film_sticking),
                "oxidation_probability": float(parameters.oxidation_probability),
                "mixing_efficiency": float(parameters.mixing_efficiency),
                "volatilization_yield": float(parameters.volatilization_yield),
                "film_sputter_yield": float(parameters.film_sputter_yield),
                "reference_energy_eV": float(parameters.reference_energy_eV),
                "precursor_fc_ratio": float(parameters.precursor_fc_ratio),
                "ion": [int(parameters.ion_atomic_number),
                        float(parameters.ion_mass_amu)],
            },
            "energy_model": "zbl-deposited-in-layer (stopping tables, no fitted law)",
            "chemisorption_probability": dict(self.chemisorption_probability),
            "species_channels": {
                "precursor": list(self.precursor_species),
                "fluorine": list(self.fluorine_species),
                "oxygen": list(self.oxygen_species),
                "inert": list(self.inert_species),
            },
        })

    # -- engine contract -------------------------------------------------

    def initial_state(self, shape=()):
        return MixedLayerSurfaceState.bare(shape)

    def validity(self, fluxes: SurfaceFluxes) -> MechanismValidity:
        mapped = (set(self.precursor_species) | set(self.fluorine_species)
                  | set(self.oxygen_species) | set(self.inert_species))
        unsupported = tuple(sorted(
            name for name, value in fluxes.neutral_flux_m2_s.items()
            if name not in mapped and np.any(np.asarray(value) > 0.0)))
        reasons = ()
        if unsupported:
            reasons = (f"unmapped nonzero neutral species: {', '.join(unsupported)}",)
        return MechanismValidity(
            within_declared_scope=not unsupported,
            reasons=reasons,
            unsupported_neutral_species=unsupported,
            known_model_form_omissions=(
                "ion spectrum compressed to flux-weighted mean energy/cosine",
                "no spontaneous (ion-free) chemical etch channel",
                "chemisorption channel uses substrate-agnostic published probabilities",
                "single projectile species for stopping tables",
                "product routing unresolved (SiF4/CO/COF2 emission laws undeclared)",
            ),
            parameter_evidence_supports_prediction=False,
            nonpredictive_parameters=(
                "sticking_probability", "fluorine_film_sticking",
                "oxidation_probability", "mixing_efficiency",
                "volatilization_yield", "film_sputter_yield"),
        )

    def neutral_reaction_probability(self, state: MixedLayerSurfaceState):
        """Per-collision loss probability for each mapped neutral species."""
        if not isinstance(state, MixedLayerSurfaceState):
            raise TypeError("neutral reaction probabilities require MixedLayerSurfaceState")
        par = self.parameters
        module = state._module_state()
        d_fc = np.asarray(module.film_thickness_nm(), dtype=float)
        theta_film = 1.0 - np.exp(-np.minimum(
            d_fc * 1e-9 * 7.5e28 / _MONOLAYER_AREAL_M2, 700.0))
        theta_f_layer = np.minimum(
            np.asarray(state.n_f, dtype=float) / _MONOLAYER_AREAL_M2, 1.0)
        film_total = np.asarray(state.n_c_film + state.n_f_film, dtype=float)
        with np.errstate(divide="ignore", invalid="ignore"):
            x_c = np.where(film_total > 0.0,
                           state.n_c_film / np.maximum(film_total, 1e-300), 0.0)
        layer_c_fraction = (np.asarray(state.n_c, dtype=float)
                            / (np.asarray(state.n_c, dtype=float) + _MONOLAYER_AREAL_M2))
        probability = {}
        for name in self.precursor_species:
            probability[name] = np.full_like(theta_film, par.sticking_probability)
        for name in self.fluorine_species:
            probability[name] = (par.fluorine_film_sticking * theta_film
                                 + (1.0 - theta_film) * (1.0 - theta_f_layer))
        for name in self.oxygen_species:
            probability[name] = par.oxidation_probability * (
                theta_film * x_c + (1.0 - theta_film) * layer_c_fraction)
        for name in self.inert_species:
            probability[name] = np.zeros_like(theta_film)
        return probability

    def advance(self, state: MixedLayerSurfaceState, fluxes: SurfaceFluxes,
                duration_s: float, *, max_step_s=None, strict=True):
        if not isinstance(state, MixedLayerSurfaceState):
            raise TypeError("advance requires MixedLayerSurfaceState")
        if not np.isfinite(duration_s) or duration_s < 0.0:
            raise ValueError("duration_s must be finite and nonnegative")
        validity = self.validity(fluxes)
        if strict and not validity.within_declared_scope:
            raise ValueError("surface mechanism outside declared scope: "
                             + "; ".join(validity.reasons))
        shape = np.asarray(state.n_f).shape
        module_fluxes = self._module_fluxes(fluxes, shape)

        if duration_s == 0.0:
            zeros = np.zeros(shape)
            return SurfaceStepResult(
                state=state,
                etch_velocity_m_s=zeros,
                formed_complex_units_m2=zeros,
                removed_complex_units_m2=zeros,
                removed_bare_formula_units_m2=zeros,
                deposited_polymer_units_m2=zeros,
                removed_polymer_units_m2=zeros,
                material_exchange=unresolved_surface_exchange(
                    removed_units_m2={}, deposited_units_m2={},
                    limitations=validity.known_model_form_omissions),
                validity=validity)

        cap = self.default_max_step_s if max_step_s is None else float(max_step_s)
        n_steps = max(1, int(np.ceil(duration_s / cap)))
        dt = duration_s / n_steps
        module_state = state._module_state()
        removed = np.zeros(shape)
        film_deposited = np.zeros(shape)
        film_removed_before = np.asarray(
            module_state.n_c_film + module_state.n_f_film, dtype=float)
        deposited_total = np.zeros(shape)
        for _ in range(n_steps):
            result = mixed_layer_step(module_state, module_fluxes, dt,
                                      self.parameters)
            removed = removed + np.asarray(result.substrate_removal_rate) * dt
            dep_rate = (self.parameters.sticking_probability
                        * np.asarray(module_fluxes.precursor_flux)
                        * (1.0 + np.asarray(module_fluxes.precursor_fc_ratio)))
            deposited_total = deposited_total + dep_rate * dt
            module_state = result.state
        film_after = np.asarray(
            module_state.n_c_film + module_state.n_f_film, dtype=float)
        film_deposited = deposited_total
        film_removed = np.maximum(
            film_removed_before + deposited_total - film_after, 0.0)

        bulk_density = (_CARBON_ATOM_DENSITY_M3
                        if self.parameters.substrate == "carbon"
                        else _SIO2_FORMULA_DENSITY_M3)
        velocity = removed / bulk_density / duration_s
        zeros = np.zeros(shape)
        def representational_floor(value):
            array = np.asarray(value, dtype=float)
            floor = -1e-6 * _MONOLAYER_AREAL_M2
            if np.any(array < floor):
                raise ValueError("mixed-layer reservoir went materially negative "
                                 "(integrator bug, not rounding)")
            return np.maximum(array, 0.0)

        new_state = MixedLayerSurfaceState(
            representational_floor(module_state.n_c_film),
            representational_floor(module_state.n_f_film),
            representational_floor(module_state.n_si),
            representational_floor(module_state.n_o),
            representational_floor(module_state.n_c),
            representational_floor(module_state.n_f),
            np.asarray(state.removed_formula_units_m2) + removed)
        inventory_name = ("carbon_atom" if self.parameters.substrate == "carbon"
                          else "sio2_formula")
        exchange = unresolved_surface_exchange(
            removed_units_m2={inventory_name: removed,
                              "fluorocarbon_film_atom": film_removed},
            deposited_units_m2={"fluorocarbon_film_atom": film_deposited},
            limitations=validity.known_model_form_omissions)
        return SurfaceStepResult(
            state=new_state,
            etch_velocity_m_s=velocity,
            formed_complex_units_m2=zeros,
            removed_complex_units_m2=zeros,
            removed_bare_formula_units_m2=removed,
            deposited_polymer_units_m2=film_deposited,
            removed_polymer_units_m2=film_removed,
            material_exchange=exchange,
            validity=validity)

    # -- helpers ---------------------------------------------------------

    def _module_fluxes(self, fluxes: SurfaceFluxes, shape) -> ModuleFluxes:
        def channel(names):
            total = np.zeros(shape)
            for name in names:
                value = fluxes.neutral_flux_m2_s.get(name)
                if value is not None:
                    total = total + np.broadcast_to(
                        np.asarray(value, dtype=float), shape)
            return total

        ion_flux = np.zeros(shape)
        energy_weighted = np.zeros(shape)
        cosine_weighted = np.zeros(shape)
        for population in fluxes.energetic_fluxes:
            if isinstance(population, FaceResolvedEnergeticFlux):
                flux = np.zeros(int(population.face_count))
                e_sum = np.zeros(int(population.face_count))
                c_sum = np.zeros(int(population.face_count))
                np.add.at(flux, population.event_face,
                          population.event_flux_m2_s)
                np.add.at(e_sum, population.event_face,
                          population.event_flux_m2_s * population.event_energy_eV)
                np.add.at(c_sum, population.event_face,
                          population.event_flux_m2_s
                          * population.event_cosine_incidence)
                flux = np.broadcast_to(flux, shape)
                e_sum = np.broadcast_to(e_sum, shape)
                c_sum = np.broadcast_to(c_sum, shape)
            elif isinstance(population, EnergeticFlux):
                base = np.broadcast_to(
                    np.asarray(population.flux_m2_s, dtype=float), shape)
                mean_e = float(np.dot(population.weight, population.energy_eV))
                mean_c = float(np.dot(population.weight,
                                      population.cosine_incidence))
                flux = base
                e_sum = base * mean_e
                c_sum = base * mean_c
            else:
                raise TypeError("unsupported energetic flux population")
            ion_flux = ion_flux + flux
            energy_weighted = energy_weighted + e_sum
            cosine_weighted = cosine_weighted + c_sum
        with np.errstate(divide="ignore", invalid="ignore"):
            safe = np.maximum(ion_flux, 1e-300)
            mean_energy = np.where(ion_flux > 0.0, energy_weighted / safe, 0.0)
            mean_cosine = np.where(ion_flux > 0.0, cosine_weighted / safe, 1.0)
        carbon_flux = np.zeros(shape)
        fluorine_bound = np.zeros(shape)
        chem_c = np.zeros(shape)
        chem_f = np.zeros(shape)
        for name, (c_atoms, f_atoms) in self.precursor_stoichiometry.items():
            value = fluxes.neutral_flux_m2_s.get(name)
            if value is None:
                continue
            base = np.broadcast_to(np.asarray(value, dtype=float), shape)
            carbon_flux = carbon_flux + base * c_atoms
            fluorine_bound = fluorine_bound + base * f_atoms
            probability = self.chemisorption_probability.get(name, 0.0)
            chem_c = chem_c + probability * base * c_atoms
            chem_f = chem_f + probability * base * f_atoms
        with np.errstate(divide="ignore", invalid="ignore"):
            fc_ratio = np.where(
                carbon_flux > 0.0,
                fluorine_bound / np.maximum(carbon_flux, 1e-300),
                self.parameters.precursor_fc_ratio)
        return ModuleFluxes(
            precursor_flux=carbon_flux,
            fluorine_flux=channel(self.fluorine_species),
            oxygen_flux=channel(self.oxygen_species),
            ion_flux=ion_flux,
            ion_energy_eV=mean_energy,
            cosine_incidence=np.clip(mean_cosine, 0.0, 1.0),
            precursor_fc_ratio=fc_ratio,
            chemisorption_carbon_flux=chem_c,
            chemisorption_fluorine_flux=chem_f)


# Krueger 2024 species channels: fluorocarbon radicals with explicit C/F
# stoichiometry (the boundary publishes NO atomic-F flux — fluorine reaches
# the layer through the film AND by direct complex-formation chemisorption).
KRUEGER_2024_PRECURSOR_STOICHIOMETRY = {
    "CF": (1.0, 1.0), "CF2": (1.0, 2.0), "CF3": (1.0, 3.0),
    "C2F3": (2.0, 3.0), "C2F4": (2.0, 4.0),
    "C3F5": (3.0, 5.0), "C3F6": (3.0, 6.0),
}

# Krueger's published complex-formation probabilities (thesis Appendix B):
# the direct-chemisorption channel delivering bound fluorine into the mixed
# layer on open oxide. Published constants — data, not knobs.
KRUEGER_2024_CHEMISORPTION_PROBABILITY = {
    "CF": 0.2729, "CF2": 0.2729, "CF3": 0.2, "C2F3": 0.2,
    "C2F4": 0.001, "C3F5": 0.001, "C3F6": 0.001,
}


def build_krueger_2024_mixed_layer_mechanisms(
        *, oxide_parameters=None, mask_parameters=None,
        volatilization_yield=1.0):
    """Oxide + mask mixed-layer mechanisms wired to the Krueger species set.

    ``volatilization_yield`` is the single base-condition-anchored absolute
    rate constant (per-channel reference yield), the same one-constant-on-base
    discipline the K24-DEKNOB-1 study declared for Lambda.
    """
    oxide_parameters = oxide_parameters or MixedLayerParams(
        substrate="sio2", volatilization_yield=float(volatilization_yield))
    mask_parameters = mask_parameters or MixedLayerParams(
        substrate="carbon", volatilization_yield=float(volatilization_yield))
    common = dict(
        precursor_species=dict(KRUEGER_2024_PRECURSOR_STOICHIOMETRY),
        fluorine_species=(),
        oxygen_species=("O",),
        inert_species=("C3F4",),
        chemisorption_probability=dict(KRUEGER_2024_CHEMISORPTION_PROBABILITY))
    return (MixedLayerMechanism(oxide_parameters, **common),
            MixedLayerMechanism(mask_parameters, **common))
