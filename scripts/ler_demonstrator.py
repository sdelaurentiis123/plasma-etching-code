"""LER demonstrator — the full measurement chain on engine-generated edges.

Synthesise a rough mask edge ensemble (Palasantzas), transfer it through the
exact ray-shadowing operator of ``ler_gate2_shadowing`` (the same geometry the
exchange operator computes), then recover the transfer function with the H1
cross-spectral estimator and report the roughness triplet on both sides.

Demonstrator scale, stated honestly: the transfer here is *geometric only* --
static mask, no chemistry coupling, no film, no feature evolution -- so the
measured T(k) is the shadowing transfer, not the full etch transfer. The
chemistry-coupled wide-y campaign (LER_DEMONSTRATOR_PLAN sec. Rung B, ~150
solves with the extrusion guard off) remains the next rung; what this run
validates is that the whole chain -- synthesis, physics, metrology, transfer
estimation -- closes on physics-generated data rather than on synthetic
ground truth.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

_trapezoid = getattr(np, "trapezoid", None) or np.trapz

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from petch.ler_metrology import (  # noqa: E402
    EdgeStatistics, averaged_psd_nm3, fit_edge_statistics, synthesize_edge_nm)
from petch.ler_transfer import estimate_transfer  # noqa: E402
from ler_gate2_shadowing import shadowed_substrate_edge_nm  # noqa: E402


def run_case(*, sigma_nm, correlation_length_nm, roughness_exponent,
             sigma_theta_deg, sidewall_angle_deg, mask_height_nm,
             n_points, spacing_nm, seeds):
    statistics = EdgeStatistics(sigma_nm, correlation_length_nm,
                                roughness_exponent)
    inputs, outputs = [], []
    for seed in seeds:
        edge = synthesize_edge_nm(statistics, n_points=n_points,
                                  spacing_nm=spacing_nm, seed=int(seed))
        etched = shadowed_substrate_edge_nm(
            edge, spacing_nm=spacing_nm, mask_height_nm=mask_height_nm,
            sidewall_angle_deg=sidewall_angle_deg,
            sigma_theta_deg=sigma_theta_deg)
        inputs.append(edge)
        outputs.append(etched)
    inputs = np.array(inputs)
    outputs = np.array(outputs)

    estimate = estimate_transfer(inputs, outputs, spacing_nm=spacing_nm)
    frequencies_in, psd_in = averaged_psd_nm3(inputs, spacing_nm=spacing_nm)
    frequencies_out, psd_out = averaged_psd_nm3(outputs, spacing_nm=spacing_nm)
    stats_in = fit_edge_statistics(frequencies_in, psd_in)
    stats_out = fit_edge_statistics(frequencies_out, psd_out)

    magnitude = estimate.magnitude_squared
    coherence = estimate.coherence
    frequencies = estimate.frequencies_per_nm
    positive = frequencies > 0.0
    # Band split at the correlation frequency 1/(2 pi xi).
    corner = 1.0 / (2.0 * np.pi * correlation_length_nm)
    low = positive & (frequencies <= corner)
    high = positive & (frequencies > corner)
    intrinsic_share = (float(_trapezoid(estimate.intrinsic_psd_nm3[positive],
                                        frequencies[positive]))
                       / float(_trapezoid(psd_out[positive],
                                          frequencies[positive])))
    return {
        "declared": {"sigma_nm": sigma_nm,
                     "correlation_length_nm": correlation_length_nm,
                     "roughness_exponent": roughness_exponent},
        "beam_sigma_theta_deg": sigma_theta_deg,
        "sidewall_angle_deg": sidewall_angle_deg,
        "n_realizations": len(seeds),
        "sigma_in_nm": float(np.mean([np.std(row) for row in inputs])),
        "sigma_out_nm": float(np.mean([np.std(row) for row in outputs])),
        "sigma_ratio": float(np.mean([np.std(row) for row in outputs])
                             / np.mean([np.std(row) for row in inputs])),
        "fitted_in": {"sigma_nm": stats_in.sigma_nm,
                      "xi_nm": stats_in.correlation_length_nm,
                      "alpha": stats_in.roughness_exponent},
        "fitted_out": {"sigma_nm": stats_out.sigma_nm,
                       "xi_nm": stats_out.correlation_length_nm,
                       "alpha": stats_out.roughness_exponent},
        "transfer": {
            "corner_frequency_per_nm": float(corner),
            "mean_T2_low_band": float(np.mean(magnitude[low])),
            "mean_T2_high_band": float(np.mean(magnitude[high])),
            "T2_at_nyquist": float(magnitude[-1]),
            "mean_coherence_low": float(np.mean(coherence[low])),
            "mean_coherence_high": float(np.mean(coherence[high])),
            "min_coherence": float(np.min(coherence[positive])),
            "intrinsic_share_of_output_power": float(intrinsic_share),
        },
    }


def main():
    common = dict(roughness_exponent=0.6, sidewall_angle_deg=86.2,
                  mask_height_nm=150.0, n_points=2048, spacing_nm=1.0,
                  seeds=tuple(range(101, 117)))
    cases = [
        dict(sigma_nm=3.0, correlation_length_nm=15.0, sigma_theta_deg=25.0),
        dict(sigma_nm=3.0, correlation_length_nm=30.0, sigma_theta_deg=25.0),
        dict(sigma_nm=3.0, correlation_length_nm=15.0, sigma_theta_deg=2.0),
    ]
    report = {"scale": "demonstrator (geometric transfer, static mask)",
              "cases": []}
    for case in cases:
        result = run_case(**case, **common)
        report["cases"].append(result)
        t = result["transfer"]
        print(f"xi={case['correlation_length_nm']:4.1f} nm  "
              f"beam={case['sigma_theta_deg']:4.1f} deg  "
              f"sigma {result['sigma_in_nm']:.3f} -> {result['sigma_out_nm']:.3f} "
              f"({result['sigma_ratio']:.4f})  "
              f"|T|^2 low={t['mean_T2_low_band']:.4f} "
              f"high={t['mean_T2_high_band']:.4f} "
              f"nyq={t['T2_at_nyquist']:.4f}  "
              f"coh_min={t['min_coherence']:.4f}  "
              f"intrinsic={t['intrinsic_share_of_output_power']:.2e}",
              flush=True)

    destination = Path("results/curated/ler_demonstrator")
    destination.mkdir(parents=True, exist_ok=True)
    with open(destination / "demonstrator.json", "w") as handle:
        json.dump(report, handle, indent=2)
    print(f"wrote {destination / 'demonstrator.json'}")


if __name__ == "__main__":
    main()
