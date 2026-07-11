"""Geometry- and species-agnostic adaptive quadrature over surface elements.

The controller knows only element weights, replicate estimates, and tolerances. It has no concepts such
as floor, corner, wall, ion, electron, or aspect ratio. An evaluator may use Sobol, Gaussian quadrature,
deterministic ordinates, or another nested rule, provided increasing ``log2_samples`` refines its estimate.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np


@dataclass(frozen=True)
class AdaptiveQuadratureResult:
    element_mean: np.ndarray
    element_stderr: np.ndarray
    element_replicates: np.ndarray
    log2_samples: np.ndarray
    total_mean: float
    total_stderr: float
    converged: bool
    rounds: int
    evaluations: int


def adaptive_surface_quadrature(
    evaluator: Callable[[np.ndarray, int, int], np.ndarray],
    n_elements: int,
    *,
    weights=None,
    base_log2: int = 10,
    max_log2: int = 16,
    n_replicates: int = 4,
    seed: int = 0,
    absolute_tolerance: float = 1e-3,
    relative_tolerance: float = 5e-3,
    element_absolute_tolerance: float | None = None,
    element_relative_tolerance: float = 0.0,
    refine_fraction: float = 0.5,
    initial_log2_samples=None,
) -> AdaptiveQuadratureResult:
    """Adapt nested per-element quadrature using replicate uncertainty as the sole error indicator.

    ``evaluator(indices, log2_samples, seed)`` returns one estimate per requested element. At each round,
    elements are ranked by ``abs(weight) * replicate_std`` and the largest contributors are refined one
    level. Convergence requires the aggregate replicate standard error to meet its mixed tolerance and,
    when requested, every element standard error to meet the mixed element tolerance
    ``element_absolute_tolerance + element_relative_tolerance*abs(element_mean)``.
    """
    if n_elements <= 0:
        raise ValueError("n_elements must be positive")
    if not 0.0 < refine_fraction <= 1.0:
        raise ValueError("refine_fraction must lie in (0, 1]")
    if n_replicates < 2:
        raise ValueError("at least two replicates are required to estimate uncertainty")
    if element_relative_tolerance < 0.0:
        raise ValueError("element_relative_tolerance must be nonnegative")
    if base_log2 < 0 or max_log2 < base_log2:
        raise ValueError("require 0 <= base_log2 <= max_log2")
    if weights is None:
        weights = np.full(n_elements, 1.0 / n_elements)
    else:
        weights = np.asarray(weights, dtype=float)
        if weights.shape != (n_elements,):
            raise ValueError("weights must have shape (n_elements,)")

    if initial_log2_samples is None:
        levels = np.full(n_elements, int(base_log2), dtype=np.int64)
    else:
        levels = np.asarray(initial_log2_samples, dtype=np.int64).copy()
        if (levels.shape != (n_elements,) or np.any(levels < base_log2)
                or np.any(levels > max_log2)):
            raise ValueError("initial sample levels must lie within the configured level range")
    estimates = np.empty((n_replicates, n_elements), dtype=float)
    evaluations = 0

    def evaluate(indices, level):
        nonlocal evaluations
        for replicate in range(n_replicates):
            values = np.asarray(evaluator(indices, int(level), int(seed + replicate)), dtype=float)
            if values.shape != indices.shape:
                raise ValueError("evaluator returned the wrong shape")
            estimates[replicate, indices] = values
        evaluations += int(indices.size) * n_replicates * (2 ** int(level))

    for initial_level in np.unique(levels):
        indices = np.where(levels == initial_level)[0]
        evaluate(indices, int(initial_level))
    rounds = 0
    converged = False
    while True:
        rounds += 1
        element_mean = estimates.mean(axis=0)
        element_stderr = estimates.std(axis=0, ddof=1) / np.sqrt(n_replicates)
        totals = estimates @ weights
        total_mean = float(totals.mean())
        total_stderr = float(totals.std(ddof=1) / np.sqrt(n_replicates))
        total_tol = float(absolute_tolerance + relative_tolerance * abs(total_mean))
        element_ok = (element_absolute_tolerance is None
                      or bool(np.all(element_stderr <= (
                          element_absolute_tolerance
                          + element_relative_tolerance * np.abs(element_mean)))))
        if total_stderr <= total_tol and element_ok:
            converged = True
            break

        refinable = np.where(levels < max_log2)[0]
        if refinable.size == 0:
            break
        contribution = np.abs(weights[refinable]) * estimates[:, refinable].std(axis=0, ddof=1)
        order = refinable[np.argsort(-contribution, kind="stable")]
        count = max(1, int(np.ceil(refine_fraction * order.size)))
        selected = order[:count]
        # Elements can have different levels. Group them so every evaluator call has one nested rule.
        for old_level in np.unique(levels[selected]):
            group = selected[levels[selected] == old_level]
            new_level = int(old_level + 1)
            evaluate(group, new_level)
            levels[group] = new_level

    return AdaptiveQuadratureResult(
        element_mean=element_mean.copy(),
        element_stderr=element_stderr.copy(),
        element_replicates=estimates.copy(),
        log2_samples=levels.copy(),
        total_mean=total_mean,
        total_stderr=total_stderr,
        converged=converged,
        rounds=rounds,
        evaluations=evaluations,
    )
