"""Transport-definition-safe grading of the direct pure-Cl2 swarm board."""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
import math
import re

from .chlorine_swarm import GonzalezMaganaPureChlorineSwarmBoard


@dataclass(frozen=True)
class ElectronSwarmPrediction:
    """One prediction carrying the deck and observable definition used."""

    observation_id: str
    transport_definition: str
    value_si: float
    si_unit: str
    solver_id: str
    collision_deck_sha256: str

    def __post_init__(self):
        if (
            not str(self.observation_id).strip()
            or not str(self.transport_definition).strip()
            or not math.isfinite(self.value_si)
            or not str(self.si_unit).strip()
            or not str(self.solver_id).strip()
            or not re.fullmatch(r"[0-9a-f]{64}", self.collision_deck_sha256)
        ):
            raise ValueError("invalid electron swarm prediction")


@dataclass(frozen=True)
class ChlorineSwarmResidual:
    observation_id: str
    observable: str
    reduced_field_Td: float
    measured_value_si: float
    predicted_value_si: float
    si_unit: str
    signed_error_si: float
    typical_interval_lower_si: float
    typical_interval_upper_si: float
    inside_source_wide_typical_interval: bool


@dataclass(frozen=True)
class GonzalezMaganaPureChlorineSwarmGrade:
    """A residual board, not automatic reactor or depth certification."""

    residuals: tuple[ChlorineSwarmResidual, ...]
    solver_id: str
    collision_deck_sha256: str

    def __post_init__(self):
        residuals = tuple(self.residuals)
        if (
            len(residuals) != 52
            or len({item.observation_id for item in residuals}) != 52
            or not str(self.solver_id).strip()
            or not re.fullmatch(r"[0-9a-f]{64}", self.collision_deck_sha256)
        ):
            raise ValueError("invalid pure-chlorine swarm grade")
        object.__setattr__(self, "residuals", residuals)

    @property
    def inside_interval_counts(self) -> dict[str, tuple[int, int]]:
        totals = Counter(item.observable for item in self.residuals)
        inside = Counter(
            item.observable
            for item in self.residuals
            if item.inside_source_wide_typical_interval
        )
        return {
            observable: (inside[observable], total)
            for observable, total in sorted(totals.items())
        }

    @property
    def all_inside_source_wide_typical_intervals(self) -> bool:
        return all(
            item.inside_source_wide_typical_interval
            for item in self.residuals
        )

    @property
    def supports_collision_set_validation(self) -> bool:
        return True

    @property
    def supports_reactor_state_prediction(self) -> bool:
        return False

    @property
    def supports_wafer_flux(self) -> bool:
        return False

    @property
    def supports_feature_depth(self) -> bool:
        return False


def grade_gonzalez_magana_pure_cl2_swarm(
    board: GonzalezMaganaPureChlorineSwarmBoard,
    predictions: tuple[ElectronSwarmPrediction, ...],
) -> GonzalezMaganaPureChlorineSwarmGrade:
    """Grade only predictions with measurement-equivalent definitions."""

    predictions = tuple(predictions)
    if len(predictions) != len(board.measurements):
        raise ValueError("swarm grade requires exactly one prediction per marker")
    by_id = {item.observation_id: item for item in predictions}
    if len(by_id) != len(predictions):
        raise ValueError("duplicate swarm prediction observation id")
    expected_ids = {item.observation_id for item in board.measurements}
    if set(by_id) != expected_ids:
        raise ValueError("swarm prediction ids do not match the direct board")
    solver_ids = {item.solver_id for item in predictions}
    deck_hashes = {item.collision_deck_sha256 for item in predictions}
    if len(solver_ids) != 1 or len(deck_hashes) != 1:
        raise ValueError("one swarm grade must use one solver and collision deck")

    residuals = []
    for measurement in board.measurements:
        prediction = by_id[measurement.observation_id]
        if prediction.transport_definition != measurement.transport_definition:
            raise ValueError(
                "swarm transport definition mismatch for "
                f"{measurement.observation_id}: expected "
                f"{measurement.transport_definition!r}, received "
                f"{prediction.transport_definition!r}"
            )
        if prediction.si_unit != measurement.si_unit:
            raise ValueError(
                f"swarm SI unit mismatch for {measurement.observation_id}")
        half_width = (
            abs(measurement.value_si)
            * measurement.relative_uncertainty_max
        )
        lower = measurement.value_si - half_width
        upper = measurement.value_si + half_width
        residuals.append(ChlorineSwarmResidual(
            observation_id=measurement.observation_id,
            observable=measurement.observable,
            reduced_field_Td=measurement.reduced_field_Td,
            measured_value_si=measurement.value_si,
            predicted_value_si=prediction.value_si,
            si_unit=measurement.si_unit,
            signed_error_si=prediction.value_si - measurement.value_si,
            typical_interval_lower_si=lower,
            typical_interval_upper_si=upper,
            inside_source_wide_typical_interval=(
                lower <= prediction.value_si <= upper),
        ))
    return GonzalezMaganaPureChlorineSwarmGrade(
        residuals=tuple(residuals),
        solver_id=next(iter(solver_ids)),
        collision_deck_sha256=next(iter(deck_hashes)),
    )
