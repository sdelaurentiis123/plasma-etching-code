import numpy as np

from petch.adaptive_quadrature import adaptive_surface_quadrature


class SyntheticEvaluator:
    def __init__(self, truth, noise):
        self.truth = np.asarray(truth, dtype=float)
        self.noise = np.asarray(noise, dtype=float)

    def __call__(self, indices, log2_samples, seed):
        # Deterministic nested-error surrogate: noisy elements converge as N^-1/2.
        phase = np.sin((indices + 1) * (seed + 1) * 1.61803398875)
        return self.truth[indices] + self.noise[indices] * phase / np.sqrt(2 ** log2_samples)


def test_adaptive_quadrature_refines_by_error_not_element_identity():
    truth = np.array([0.2, 0.4, 0.6, 0.8])
    evaluator = SyntheticEvaluator(truth, [0.01, 0.7, 0.02, 0.03])
    result = adaptive_surface_quadrature(
        evaluator, 4, base_log2=4, max_log2=12, n_replicates=4,
        absolute_tolerance=2e-3, relative_tolerance=0.0,
        element_absolute_tolerance=4e-3, refine_fraction=0.5,
    )
    assert result.converged
    assert result.log2_samples[1] > result.log2_samples[[0, 2, 3]].min()
    assert np.isclose(result.total_mean, truth.mean(), atol=3e-3)


def test_adaptive_quadrature_is_reproducible():
    evaluator = SyntheticEvaluator([0.1, 0.3, 0.9], [0.2, 0.1, 0.4])
    kwargs = dict(base_log2=3, max_log2=10, n_replicates=5, seed=17,
                  absolute_tolerance=1e-3, relative_tolerance=0.0)
    first = adaptive_surface_quadrature(evaluator, 3, **kwargs)
    second = adaptive_surface_quadrature(evaluator, 3, **kwargs)
    assert first.total_mean == second.total_mean
    assert first.total_stderr == second.total_stderr
    assert np.array_equal(first.log2_samples, second.log2_samples)


def test_adaptive_quadrature_reports_unmet_tolerance_at_budget_limit():
    evaluator = SyntheticEvaluator([1.0, 2.0], [10.0, 10.0])
    result = adaptive_surface_quadrature(
        evaluator, 2, base_log2=2, max_log2=2, n_replicates=3,
        absolute_tolerance=1e-12, relative_tolerance=0.0,
    )
    assert not result.converged
    assert np.all(result.log2_samples == 2)
