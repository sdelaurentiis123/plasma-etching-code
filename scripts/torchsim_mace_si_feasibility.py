#!/usr/bin/env python3
# /// script
# requires-python = ">=3.12"
# dependencies = [
#   "torch-sim-atomistic[mace,io]==0.6.1",
# ]
# ///
"""Probe the real TorchSim/MACE path on bulk silicon.

The calculation checks batched MACE-MP inference on a small silicon equation of
state and a short deterministic NVE trajectory.  It does not claim coverage of
ions, radicals, surfaces, charge transfer, or bond-breaking etch products.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import platform
import time

from ase.build import bulk
from mace.calculators.foundations_models import (
    download_mace_mp_checkpoint,
    mace_mp,
    mace_mp_urls,
)
import torch
import torch_sim as ts
from torch_sim.models.mace import MaceModel
from torch_sim.units import MetalUnits


DEFAULT_LATTICE_CONSTANTS_ANGSTROM = (5.20, 5.43, 5.70)


def _device(name: str):
    if name == "auto":
        name = "cuda" if torch.cuda.is_available() else "cpu"
    if name == "mps":
        raise RuntimeError(
            "TorchSim 0.6.1's default Warp neighbor list does not support MPS")
    device = torch.device(name)
    if device.type == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("CUDA was requested but is unavailable")
    return device


def run(lattice_constants, *, steps: int, timestep_ps: float, temperature_k: float,
        device: torch.device):
    dtype = torch.float32
    torch.manual_seed(0)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(0)
    torch.use_deterministic_algorithms(True)

    structures = [
        bulk("Si", "diamond", a=lattice, cubic=True).repeat((2, 2, 2))
        for lattice in lattice_constants]
    state = ts.initialize_state(structures, device=device, dtype=dtype)

    load_started = time.perf_counter()
    checkpoint_path = Path(download_mace_mp_checkpoint("small"))
    checkpoint_sha256 = sha256(checkpoint_path.read_bytes()).hexdigest()
    raw_model = mace_mp(
        model=checkpoint_path, return_raw_model=True, default_dtype="float32",
        device=str(device))
    model_load_s = time.perf_counter() - load_started
    model = MaceModel(
        model=raw_model, device=device, dtype=dtype,
        compute_forces=True, compute_stress=False, enable_cueq=False)

    inference_started = time.perf_counter()
    static = model(state)
    if device.type == "cuda":
        torch.cuda.synchronize()
    static_inference_s = time.perf_counter() - inference_started
    counts = torch.bincount(state.system_idx).detach().cpu()
    energy_per_atom = static["energy"].detach().cpu() / counts
    force_norms = torch.linalg.vector_norm(static["forces"], dim=1)
    max_force = torch.segment_reduce(
        force_norms, reduce="max", lengths=torch.bincount(state.system_idx))

    state.rng = 0
    state = ts.nve_init(
        state=state, model=model,
        kT=torch.full(
            (len(lattice_constants),), temperature_k * MetalUnits.temperature,
            device=device, dtype=dtype))
    initial_total = (
        state.energy
        + ts.calc_kinetic_energy(
            masses=state.masses, momenta=state.momenta,
            system_idx=state.system_idx)).detach().cpu()
    dt = torch.tensor(timestep_ps * MetalUnits.time, device=device, dtype=dtype)
    dynamics_started = time.perf_counter()
    for _ in range(steps):
        state = ts.nve_step(state=state, model=model, dt=dt)
    if device.type == "cuda":
        torch.cuda.synchronize()
    dynamics_s = time.perf_counter() - dynamics_started
    final_total = (
        state.energy
        + ts.calc_kinetic_energy(
            masses=state.masses, momenta=state.momenta,
            system_idx=state.system_idx)).detach().cpu()

    cases = []
    for index, lattice in enumerate(lattice_constants):
        drift = float(final_total[index] - initial_total[index])
        cases.append({
            "lattice_constant_angstrom": lattice,
            "atoms": int(counts[index]),
            "static_energy_eV_per_atom": float(energy_per_atom[index]),
            "static_max_force_eV_per_angstrom": float(max_force[index].detach().cpu()),
            "initial_nve_total_energy_eV": float(initial_total[index]),
            "final_nve_total_energy_eV": float(final_total[index]),
            "nve_energy_drift_eV": drift,
            "nve_relative_energy_drift": drift / max(abs(float(initial_total[index])), 1.0e-12),
        })
    minimum_case = min(cases, key=lambda case: case["static_energy_eV_per_atom"])
    return {
        "schema_version": 1,
        "experiment": "torchsim_mace_mp_bulk_si_feasibility",
        "predictive_etch_physics": False,
        "supported_claim": (
            "TorchSim can batch and integrate a real MACE potential for crystalline Si."),
        "unsupported_claims": [
            "ion or radical impact response",
            "surface reaction probabilities",
            "sputter or reactive etch yields",
            "charged or electronically excited states",
            "transfer to SiO2, fluorocarbon films, masks, or reactor conditions",
        ],
        "software": {
            "torch_version": torch.__version__,
            "torchsim_version": ts.__version__,
            "python_version": platform.python_version(),
            "device": str(device),
            "dtype": str(dtype).removeprefix("torch."),
            "potential": "MACE-MP small 2023-12-10 checkpoint",
            "potential_url": mace_mp_urls["small"],
            "potential_sha256": checkpoint_sha256,
        },
        "configuration": {
            "lattice_constants_angstrom": list(lattice_constants),
            "supercell_repetition": [2, 2, 2],
            "temperature_k": temperature_k,
            "steps": steps,
            "timestep_ps": timestep_ps,
            "integrator": "NVE velocity Verlet",
            "seed": 0,
            "deterministic_algorithms": True,
        },
        "performance": {
            "model_load_s": model_load_s,
            "batched_static_inference_s": static_inference_s,
            "batched_dynamics_s": dynamics_s,
            "system_steps_per_s": steps * len(lattice_constants) / dynamics_s,
        },
        "differentiability": {
            "adapter_energy_requires_grad": bool(static["energy"].requires_grad),
            "adapter_forces_require_grad": bool(static["forces"].requires_grad),
            "adapter_autograd_gate_passed": bool(
                static["energy"].requires_grad and static["forces"].requires_grad),
            "note": (
                "MACE 0.3.16's TorchSim adapter detaches returned energy and forces; "
                "this tested adapter is inference-only for autograd purposes."),
        },
        "lowest_sampled_energy_lattice_constant_angstrom": (
            minimum_case["lattice_constant_angstrom"]),
        "cases": cases,
    }


def _parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--device", choices=("auto", "cpu", "mps", "cuda"), default="cpu")
    parser.add_argument("--steps", type=int, default=20)
    parser.add_argument("--timestep-ps", type=float, default=0.001)
    parser.add_argument("--temperature-k", type=float, default=300.0)
    parser.add_argument(
        "--lattice-constants-angstrom", nargs="+", type=float,
        default=DEFAULT_LATTICE_CONSTANTS_ANGSTROM)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main():
    arguments = _parse_args()
    lattice_constants = tuple(sorted(set(arguments.lattice_constants_angstrom)))
    if (len(lattice_constants) < 2 or any(value <= 0.0 for value in lattice_constants)
            or arguments.steps <= 0 or arguments.timestep_ps <= 0.0
            or arguments.temperature_k < 0.0):
        raise ValueError("invalid lattice sweep or dynamics controls")
    payload = run(
        lattice_constants, steps=arguments.steps,
        timestep_ps=arguments.timestep_ps, temperature_k=arguments.temperature_k,
        device=_device(arguments.device))
    encoded = json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n"
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(encoded, encoding="utf-8")
    print(json.dumps({
        "output": str(arguments.output),
        "sha256": sha256(encoded.encode()).hexdigest(),
        "lowest_sampled_energy_lattice_constant_angstrom": (
            payload["lowest_sampled_energy_lattice_constant_angstrom"]),
        "performance": payload["performance"],
        "cases": payload["cases"],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
