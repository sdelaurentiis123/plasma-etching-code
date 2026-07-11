"""Forward/backward reciprocity audit in a frozen, nonuniform AR4 charging field.

This is a numerical-invariant gate, not an experimental fit.  The field is generated once by the
deterministic charging solver, frozen, and then scored independently by (1) surface-launched adjoint
gathers and (2) plasma-plane forward particles.  Fluxes are normalized to the same incident flux per
unit horizontal length.  Run the W ladder before changing a particle mover::

    python scripts/backward_actual_field_reciprocity.py --width 16
    python scripts/backward_actual_field_reciprocity.py --width 32
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np
from scipy.stats import gamma as gamma_dist, norm, qmc

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from petch.charging2d import _build_edge_array_geometry
from petch.charging_backward import backward_electron_gather, backward_ion_gather, self_consistent_backward
from petch.charging_general import _trace_general


def _forward_floor_flux(result, geometry, species, n_log2, seed, trace_dt, trace_dt_field):
    nx, nz = geometry["nx"], geometry["nz"]
    t0, t1 = geometry["trench0"], geometry["trench1"]
    sampler = qmc.Sobol(d=4 if species == "electron" else 3, scramble=True, seed=seed)
    u = sampler.random_base2(n_log2)
    if species == "electron":
        energy = gamma_dist.ppf(u[:, 0], a=2.0, scale=4.0)
        ct = np.sqrt(u[:, 1])
        vx = np.sqrt(energy) * np.sqrt(1.0 - ct * ct) * np.cos(2.0 * np.pi * u[:, 2])
        vz = np.sqrt(energy) * ct
        xcoord = u[:, 3]
        charge = -1.0
    else:
        phase = 2.0 * np.pi * u[:, 0]
        vx = np.sqrt(0.25) * norm.ppf(np.clip(u[:, 1], 1e-6, 1.0 - 1e-6))
        vz = np.sqrt(2.0 + 37.0 + 30.0 * np.sin(phase))
        xcoord = u[:, 2]
        charge = 1.0
    x = xcoord * nx
    z = np.full_like(x, 0.51)
    hit_x, hit_z, *_ = _trace_general(
        result["Ex"], result["Ez"], geometry["solid"], x, z, vx, vz, charge, nx, nz,
        200 * nz, trace_dt, trace_dt_field,
    )
    floor_hit = (hit_z == nz - 1) & (hit_x >= t0) & (hit_x < t1)
    # Convert probability under a source uniform over nx to flux per trench-opening width.
    return float(floor_hit.mean() * nx / (t1 - t0))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--width", type=int, default=16)
    parser.add_argument("--mouth", type=int, default=None)
    parser.add_argument("--charge-log2", type=int, default=10)
    parser.add_argument("--score-log2", type=int, default=17)
    parser.add_argument("--iterations", type=int, default=10)
    parser.add_argument("--trace-dt", type=float, default=0.15)
    parser.add_argument("--trace-dt-field", type=float, default=0.10)
    parser.add_argument("--exit-energy-mixture", type=float, default=0.2)
    args = parser.parse_args()

    mouth = args.mouth if args.mouth is not None else 5 * args.width
    geometry = _build_edge_array_geometry(4.0, W=args.width, mouth=mouth)
    result = self_consistent_backward(
        geometry, n_iter=args.iterations, n_log2=args.charge_log2, n_scramble=2,
        ion_ied_phase_exponent=0.0,
    )
    t0, t1, fz = geometry["trench0"], geometry["trench1"], geometry["nz"] - 1
    cells = [(x, fz) for x in range(t0, t1)]
    normals = [(0.0, -1.0)] * len(cells)
    backward_e = float(backward_electron_gather(
        geometry["solid"], result["Ex"], result["Ez"], result["Vs"], cells, normals,
        n_log2=args.charge_log2 + 2, n_scramble=3, seed=101,
        trace_dt=args.trace_dt, trace_dt_field=args.trace_dt_field,
    ).mean())
    backward_i = float(backward_ion_gather(
        geometry["solid"], result["Ex"], result["Ez"], result["Vs"], cells, normals,
        n_log2=args.charge_log2 + 2, n_scramble=3, seed=103, ied_phase_exponent=0.0,
        trace_dt=args.trace_dt, trace_dt_field=args.trace_dt_field,
    ).mean())
    backward_i_exit = float(backward_ion_gather(
        geometry["solid"], result["Ex"], result["Ez"], result["Vs"], cells, normals,
        n_log2=args.charge_log2 + 2, n_scramble=3, seed=103, ied_phase_exponent=0.0,
        exit_state_weight=True,
        exit_energy_mixture=args.exit_energy_mixture,
        trace_dt=args.trace_dt, trace_dt_field=args.trace_dt_field,
    ).mean())
    forward_e = _forward_floor_flux(result, geometry, "electron", args.score_log2, 107,
                                    args.trace_dt, args.trace_dt_field)
    forward_i = _forward_floor_flux(result, geometry, "ion", args.score_log2, 109,
                                    args.trace_dt, args.trace_dt_field)

    print(f"W={args.width} mouth={mouth} iterations={result['iterations']} floor={result['floor_mean']:.3f} V "
          f"trace_dt={args.trace_dt:g}/{args.trace_dt_field:g}")
    print(f"field residual rms={result['field_final']['rms']:.3e}; "
          f"charge residual rms={result['balance_preupdate']['rms_log_ratio']:.3e}")
    for name, backward, forward in (("electron", backward_e, forward_e),
                                    ("ion-1d", backward_i, forward_i),
                                    (f"ion-exit-{args.exit_energy_mixture:g}", backward_i_exit, forward_i)):
        rel = (backward / forward - 1.0) if forward else np.nan
        print(f"{name:8s} backward={backward:.6f} forward={forward:.6f} relative={rel:+.2%}")


if __name__ == "__main__":
    main()
