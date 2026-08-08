"""Explicit legacy-SIGLO chlorine source-replay adapter.

The raw LXCat/SIGLO bytes remain user supplied and outside the package.  This
adapter recognizes one hash-locked 2013 deck, reproduces BOLOS/BOLSIG edge
padding as a declared numerical convention, declares double-ionization
multiplicity, and maps every inelastic row into the six-species chlorine
ledger.  It is a source-replay/sensitivity provider, not current collision
evidence or a validated reactor mechanism.
"""
from __future__ import annotations

from dataclasses import dataclass, replace
import hashlib
import json
from pathlib import Path

from .argon import ELECTRON_MASS_AMU
from .chlorine import (
    CHLORINE_ATOM_MASS_AMU,
    lee_lieberman_chlorine_species,
)
from .chlorine_wall import (
    ChlorineWallRecombinationBoundary,
    LogLinearChlorineWallRecombinationProvider,
)
from .electron_collision_chemistry import (
    ElectronCollisionChemistry,
    ElectronCollisionHeavyMapping,
)
from .electron_collision_deck import (
    ElectronCollisionDeck,
    ElectronCollisionProcess,
    load_bolsig_lxcat_file,
)
from .evaluated_chlorine import (
    nist_hayes_atomic_chlorine_ionization_collision_process,
)


LEGACY_SIGLO_CL2_2013_SHA256 = (
    "b8b1ff807d0586795dcabf4511a2a44b8b04828626b45e0fd95a319817c01063"
)
COMSOL_64_ATOMIC_CL_MOMENTUM_SHA256 = (
    "6c1391eda01d90ac504e73e3f5306b510e8a5d721860b03b50cdba39edd29661"
)
_EXPECTED_ROWS = (
    ("ATTACHMENT", "Cl^- + Cl", 0.0),
    ("ELASTIC", None, None),
    ("EXCITATION", "Cl2(v)", 0.069),
    ("EXCITATION", "Cl2(v)", 0.139),
    ("EXCITATION", "Cl2(3PI_u)", 3.36),
    ("EXCITATION", "Cl2(1PI_u)", 4.3),
    ("EXCITATION", "Cl2(3PI_g)", 6.38),
    ("EXCITATION", "Cl2(1PI_g)", 7.01),
    ("EXCITATION", "Cl2(3SIG_u)", 7.02),
    ("EXCITATION", "Cl2(1PI_ub)", 10.54),
    ("EXCITATION", "Cl2(1SIG_ub)", 10.7),
    ("EXCITATION", "Cl^- + Cl^+", 11.0),
    ("IONIZATION", "Cl2^+", 11.49),
    ("IONIZATION", "Cl^+ + Cl", 11.49),
    ("IONIZATION", "Cl2^++", 35.5),
    ("IONIZATION", "Cl^++ + Cl", 43.5),
)


@dataclass(frozen=True)
class LegacySigloChlorineReplay:
    raw_payload_sha256: str
    derived_deck: ElectronCollisionDeck
    collision_chemistry: ElectronCollisionChemistry
    maximum_energy_eV: float
    missing_reactor_channels: tuple[str, ...] = (
        "atomic_chlorine_ionization",
        "electron_detachment_from_Clminus",
        "tracked_vibrational_and_electronic_state_kinetics",
    )
    declared_sensitivity_closures: tuple[str, ...] = (
        "constant edge-value cross-section padding to solver maximum energy",
        "elastic low-energy edge padded with first printed value",
        "Cl2++ and Cl++ channels collapse to two tracked Cl+ products",
    )
    supports_direct_swarm_grade: bool = False
    supports_reactor_state_prediction: bool = False
    supports_wafer_flux: bool = False
    supports_feature_depth: bool = False


@dataclass(frozen=True)
class LegacySigloComsolChlorineReplay:
    """Source replay with calculated Cl momentum and measured ionization."""

    molecular_replay: LegacySigloChlorineReplay
    atomic_momentum_payload_sha256: str
    derived_deck: ElectronCollisionDeck
    collision_chemistry: ElectronCollisionChemistry
    maximum_energy_eV: float
    missing_reactor_channels: tuple[str, ...] = (
        "electron_detachment_from_Clminus",
        "tracked_vibrational_and_electronic_state_kinetics",
    )
    declared_sensitivity_closures: tuple[str, ...] = (
        "official COMSOL calculated atomic-Cl momentum row",
        "constant final atomic-Cl momentum cross section through grid maximum",
        "NIST Hayes atomic-Cl ionization with sub-threshold values zeroed",
    )
    supports_direct_swarm_grade: bool = False
    supports_reactor_state_prediction: bool = False
    supports_wafer_flux: bool = False
    supports_feature_depth: bool = False


@dataclass(frozen=True)
class DeclaredLogLinearWallExtrapolation:
    """Sensitivity-only extension of a direct-data wall regression.

    The wrapped provider retains the measured domain. This adapter permits a
    wider, explicitly printed composition interval only when the same fitted
    law remains a physical probability. It can never support prediction.
    """

    provider: LogLinearChlorineWallRecombinationProvider
    allowed_cl_to_cl2_ratio: tuple[float, float]
    source: str
    name: str = "declared_log_linear_wall_extrapolation"
    version: str = "1"

    def __post_init__(self):
        if not isinstance(
            self.provider, LogLinearChlorineWallRecombinationProvider
        ):
            raise TypeError("a log-linear wall provider is required")
        try:
            lower, upper = (
                float(value) for value in self.allowed_cl_to_cl2_ratio)
        except (TypeError, ValueError):
            raise ValueError(
                "wall extrapolation interval must contain two numbers"
            ) from None
        endpoint_probabilities = tuple(
            10.0 ** (
                self.provider.intercept_log10
                + self.provider.slope_per_ratio * ratio)
            for ratio in (lower, upper)
        )
        if (
            lower < 0.0
            or upper <= lower
            or lower > self.provider.valid_cl_to_cl2_ratio[0]
            or upper < self.provider.valid_cl_to_cl2_ratio[1]
            or any(
                not 0.0 < probability <= 1.0
                for probability in endpoint_probabilities)
            or not str(self.source).strip()
            or not str(self.name).strip()
            or not str(self.version).strip()
        ):
            raise ValueError("invalid physical wall extrapolation")
        object.__setattr__(
            self, "allowed_cl_to_cl2_ratio", (lower, upper))

    @property
    def supports_prediction(self) -> bool:
        return False

    def predict(
        self,
        *,
        cl_to_cl2_ratio: float,
        pressure_Pa: float,
        icp_power_W: float,
        gas_temperature_K: float,
    ) -> ChlorineWallRecombinationBoundary:
        ratio = float(cl_to_cl2_ratio)
        if not self.allowed_cl_to_cl2_ratio[0] <= ratio <= (
            self.allowed_cl_to_cl2_ratio[1]
        ):
            raise ValueError(
                "Cl/Cl2 ratio is outside the declared extrapolation interval: "
                f"ratio={ratio:.8g}, interval={self.allowed_cl_to_cl2_ratio}")
        # Composition is the only extended coordinate. Preserve every other
        # applicability gate from the direct Stafford marker board.
        domains = {
            "pressure": (float(pressure_Pa), self.provider.valid_pressure_Pa),
            "ICP power": (
                float(icp_power_W), self.provider.valid_icp_power_W),
            "gas temperature": (
                float(gas_temperature_K),
                self.provider.valid_gas_temperature_K,
            ),
        }
        for quantity, (value, domain) in domains.items():
            if not domain[0] <= value <= domain[1]:
                raise ValueError(
                    f"{quantity} is outside the source wall-provider domain")
        probability = 10.0 ** (
            self.provider.intercept_log10
            + self.provider.slope_per_ratio * ratio)
        return ChlorineWallRecombinationBoundary(
            recombination_probability=probability,
            surface_state=self.provider.surface_state,
            source=self.source,
            evidence_kind="sensitivity",
            valid_cl_to_cl2_ratio=self.allowed_cl_to_cl2_ratio,
            valid_pressure_Pa=self.provider.valid_pressure_Pa,
            valid_icp_power_W=self.provider.valid_icp_power_W,
            valid_gas_temperature_K=self.provider.valid_gas_temperature_K,
            relative_measurement_uncertainty=None,
            provenance={
                **self.provider.provenance,
                "fit_form": (
                    "log10(gamma) = intercept + slope * nCl/nCl2"),
                "slope_per_ratio": self.provider.slope_per_ratio,
                "intercept_log10": self.provider.intercept_log10,
                "direct_marker_ratio_domain": (
                    self.provider.valid_cl_to_cl2_ratio),
                "declared_extrapolation_ratio_domain": (
                    self.allowed_cl_to_cl2_ratio),
                "coefficient_selection_target": None,
            },
        )


def _verify_rows(deck: ElectronCollisionDeck) -> None:
    actual = tuple(
        (process.kind, process.product, process.energy_loss_eV)
        for process in deck.processes
    )
    if actual != _EXPECTED_ROWS:
        raise RuntimeError(
            "legacy SIGLO chlorine process topology/signatures changed"
        )


def derive_legacy_siglo_cl2_replay(
    raw_deck: ElectronCollisionDeck,
    *,
    maximum_energy_eV: float = 200.0,
) -> LegacySigloChlorineReplay:
    """Build the explicit source-replay deck from an already hash-gated deck."""
    maximum = float(maximum_energy_eV)
    if raw_deck.payload_sha256 != LEGACY_SIGLO_CL2_2013_SHA256:
        raise RuntimeError("legacy SIGLO chlorine raw hash mismatch")
    if maximum <= 100.0:
        raise ValueError("legacy SIGLO replay maximum energy must exceed 100 eV")
    _verify_rows(raw_deck)
    processes = []
    for index, process in enumerate(raw_deck.processes):
        energy = list(process.electron_energy_eV)
        cross_section = list(process.cross_section_m2)
        if process.kind in {"ELASTIC", "MOMENTUM", "EFFECTIVE"} and energy[0] > 0.0:
            energy.insert(0, 0.0)
            cross_section.insert(0, cross_section[0])
        if energy[-1] < maximum:
            energy.append(maximum)
            cross_section.append(cross_section[-1])
        processes.append(replace(
            process,
            electron_energy_eV=tuple(energy),
            cross_section_m2=tuple(cross_section),
            electron_number_change=(
                2 if index in {14, 15} else process.electron_number_change),
        ))
    derivation = {
        "schema": "petch.legacy_siglo_cl2_replay.v1",
        "raw_payload_sha256": raw_deck.payload_sha256,
        "maximum_energy_eV": maximum,
        "low_edge": "constant_first_value_for_momentum",
        "high_edge": "constant_last_value_for_all_processes",
        "double_ionization_rows": [14, 15],
    }
    derived_sha = hashlib.sha256(json.dumps(
        derivation, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")).hexdigest()
    derived = ElectronCollisionDeck(
        processes=tuple(processes),
        payload_sha256=derived_sha,
        source_database="legacy SIGLO Cl2 2013 explicit source replay",
        retrieved_at=raw_deck.retrieved_at,
        source_reference=(
            f"{raw_deck.source_reference}; raw_sha256={raw_deck.payload_sha256}; "
            f"derivation={json.dumps(derivation, sort_keys=True)}"
        ),
        packaged_or_redistributed=False,
    )
    common_source = (
        "legacy SIGLO Cl2 2013 source-replay mapping; no reactor fit"
    )
    mappings = []

    def mapping(index, name, reactants, products, evidence="published_compilation"):
        mappings.append(ElectronCollisionHeavyMapping(
            process_index=index,
            reaction_name=name,
            heavy_reactants=reactants,
            heavy_products=products,
            source=common_source,
            evidence_kind=evidence,
        ))

    mapping(0, "dissociative_attachment", {"Cl2": 1}, {"Cl-": 1, "Cl": 1})
    mapping(2, "vibrational_excitation_069", {"Cl2": 1}, {"Cl2": 1})
    mapping(3, "vibrational_excitation_139", {"Cl2": 1}, {"Cl2": 1})
    for index, state in zip(range(4, 9), (
        "3PI_u", "1PI_u", "3PI_g", "1PI_g", "3SIG_u",
    )):
        mapping(
            index,
            f"dissociative_excitation_{state}",
            {"Cl2": 1},
            {"Cl": 2},
        )
    mapping(9, "rydberg_excitation_1PI_ub", {"Cl2": 1}, {"Cl2": 1})
    mapping(10, "rydberg_excitation_1SIG_ub", {"Cl2": 1}, {"Cl2": 1})
    mapping(11, "ion_pair_formation", {"Cl2": 1}, {"Cl-": 1, "Cl+": 1})
    mapping(12, "molecular_ionization", {"Cl2": 1}, {"Cl2+": 1})
    mapping(13, "dissociative_ionization", {"Cl2": 1}, {"Cl+": 1, "Cl": 1})
    mapping(
        14,
        "double_molecular_ionization_collapsed_to_two_Clplus",
        {"Cl2": 1},
        {"Cl+": 2},
        evidence="sensitivity",
    )
    mapping(
        15,
        "double_atomic_ionization_collapsed_to_two_Clplus",
        {"Cl2": 1},
        {"Cl+": 2},
        evidence="sensitivity",
    )
    chemistry = ElectronCollisionChemistry(
        derived, lee_lieberman_chlorine_species(), tuple(mappings))
    return LegacySigloChlorineReplay(
        raw_payload_sha256=raw_deck.payload_sha256,
        derived_deck=derived,
        collision_chemistry=chemistry,
        maximum_energy_eV=maximum,
    )


def load_legacy_siglo_cl2_replay(
    path: str | Path,
    *,
    maximum_energy_eV: float = 200.0,
) -> LegacySigloChlorineReplay:
    raw = load_bolsig_lxcat_file(
        path,
        source_database="SIGLO LXCat legacy chlorine",
        retrieved_at="2013-06-04",
        source_reference=str(Path(path)),
        target="Cl2",
        expected_sha256=LEGACY_SIGLO_CL2_2013_SHA256,
    )
    return derive_legacy_siglo_cl2_replay(
        raw, maximum_energy_eV=maximum_energy_eV)


def _parse_two_column_cross_section(payload: bytes) -> tuple[
    tuple[float, ...], tuple[float, ...]
]:
    rows = []
    for line in payload.decode("utf-8-sig").splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        fields = stripped.split()
        if len(fields) != 2:
            raise RuntimeError("atomic momentum table is not two-column data")
        try:
            rows.append((float(fields[0]), float(fields[1])))
        except ValueError:
            raise RuntimeError(
                "atomic momentum table contains nonnumeric data") from None
    if (
        len(rows) != 28
        or rows[0][0] != 0.0
        or rows[-1][0] != 25.0
        or any(
            right[0] <= left[0]
            for left, right in zip(rows, rows[1:])
        )
        or any(energy < 0.0 or cross_section < 0.0
               for energy, cross_section in rows)
    ):
        raise RuntimeError("atomic momentum table topology changed")
    return (
        tuple(energy for energy, _ in rows),
        tuple(cross_section for _, cross_section in rows),
    )


def _derive_legacy_siglo_comsol_chlorine_replay(
    molecular_replay: LegacySigloChlorineReplay,
    *,
    atomic_momentum_energy_eV: tuple[float, ...],
    atomic_momentum_cross_section_m2: tuple[float, ...],
    atomic_momentum_payload_sha256: str,
) -> LegacySigloComsolChlorineReplay:
    """Compose already verified source components without reactor fitting."""
    if not isinstance(molecular_replay, LegacySigloChlorineReplay):
        raise TypeError("a molecular chlorine replay is required")
    maximum = molecular_replay.maximum_energy_eV
    energies = tuple(float(value) for value in atomic_momentum_energy_eV)
    cross_sections = tuple(
        float(value) for value in atomic_momentum_cross_section_m2)
    if energies[-1] < maximum:
        energies = (*energies, maximum)
        cross_sections = (*cross_sections, cross_sections[-1])
    atomic_momentum = ElectronCollisionProcess(
        kind="MOMENTUM",
        target="Cl",
        product=None,
        electron_energy_eV=energies,
        cross_section_m2=cross_sections,
        mass_ratio=ELECTRON_MASS_AMU / CHLORINE_ATOM_MASS_AMU,
        comments=(
            "official COMSOL 6.4 chlorine global-model asset",
            "calculated atomic-Cl elastic momentum-transfer source replay",
            "constant final value padded to solver maximum",
            "no reactor or feature fit",
        ),
    )
    atomic_ionization = (
        nist_hayes_atomic_chlorine_ionization_collision_process())
    processes = (
        *molecular_replay.derived_deck.processes,
        atomic_momentum,
        atomic_ionization,
    )
    derivation = {
        "schema": "petch.legacy_siglo_comsol_chlorine_replay.v1",
        "molecular_derived_sha256": (
            molecular_replay.derived_deck.payload_sha256),
        "atomic_momentum_payload_sha256": atomic_momentum_payload_sha256,
        "atomic_momentum_high_edge": "constant_last_value",
        "atomic_ionization": (
            "Christophorou-Olthoff Table 25 Hayes; NIST ASD threshold"),
        "maximum_energy_eV": maximum,
    }
    derived_sha = hashlib.sha256(json.dumps(
        derivation, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")).hexdigest()
    deck = ElectronCollisionDeck(
        processes=processes,
        payload_sha256=derived_sha,
        source_database=(
            "legacy SIGLO Cl2 plus official COMSOL Cl momentum plus "
            "NIST Hayes Cl ionization source replay"),
        retrieved_at="2026-08-08",
        source_reference=json.dumps(derivation, sort_keys=True),
        packaged_or_redistributed=False,
    )
    mappings = tuple(molecular_replay.collision_chemistry.mappings) + (
        ElectronCollisionHeavyMapping(
            process_index=len(processes) - 1,
            reaction_name="atomic_chlorine_ionization",
            heavy_reactants={"Cl": 1},
            heavy_products={"Cl+": 1},
            source=(
                "Christophorou-Olthoff 1999 Table 25 Hayes measurements; "
                "NIST ASD threshold; no reactor fit"),
            evidence_kind="measured",
        ),
    )
    chemistry = ElectronCollisionChemistry(
        deck, lee_lieberman_chlorine_species(), mappings)
    return LegacySigloComsolChlorineReplay(
        molecular_replay=molecular_replay,
        atomic_momentum_payload_sha256=atomic_momentum_payload_sha256,
        derived_deck=deck,
        collision_chemistry=chemistry,
        maximum_energy_eV=maximum,
    )


def load_legacy_siglo_comsol_chlorine_replay(
    molecular_cl2_path: str | Path,
    atomic_cl_momentum_path: str | Path,
    *,
    maximum_energy_eV: float = 200.0,
) -> LegacySigloComsolChlorineReplay:
    """Load a rights-safe two-target source replay from local raw assets."""
    molecular = load_legacy_siglo_cl2_replay(
        molecular_cl2_path, maximum_energy_eV=maximum_energy_eV)
    payload = Path(atomic_cl_momentum_path).read_bytes()
    digest = hashlib.sha256(payload).hexdigest()
    if digest != COMSOL_64_ATOMIC_CL_MOMENTUM_SHA256:
        raise RuntimeError("official COMSOL atomic-Cl momentum hash mismatch")
    energies, cross_sections = _parse_two_column_cross_section(payload)
    return _derive_legacy_siglo_comsol_chlorine_replay(
        molecular,
        atomic_momentum_energy_eV=energies,
        atomic_momentum_cross_section_m2=cross_sections,
        atomic_momentum_payload_sha256=digest,
    )
