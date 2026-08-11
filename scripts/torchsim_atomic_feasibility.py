#!/usr/bin/env python3
# /// script
# requires-python = ">=3.12"
# dependencies = [
#   "torch-sim-atomistic==0.6.1",
# ]
# ///
"""Exercise TorchSim as a deterministic atomistic-to-petch boundary provider.

This is deliberately a software/contract feasibility experiment, not an etch
prediction.  It fires one Lennard-Jones argon projectile at each member of a
batched, zero-temperature argon slab ensemble and exports the observed event
counts in the payload accepted by :class:`petch.SurfaceInteractionTable`.

The analytic Lennard-Jones potential is intentionally incapable of representing
reactive semiconductor etching.  Keeping that limitation in the generated
provenance prevents a successful software spike from becoming a physics claim.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
import math
from pathlib import Path
import platform
import time

import torch
import torch_sim as ts
from torch_sim.constraints import FixAtoms
from torch_sim.models.lennard_jones import LennardJonesModel
from torch_sim.units import MetalUnits


ARGON_ATOMIC_NUMBER = 18
ARGON_MASS_AMU = 39.948
ARGON_SIGMA_ANGSTROM = 3.405
ARGON_EPSILON_EV = 0.0104
ARGON_LATTICE_ANGSTROM = 5.26
DEFAULT_ENERGIES_EV = (0.02, 0.05, 0.10, 0.20)


def _device(name: str) -> torch.device:
    if name == "auto":
        if torch.cuda.is_available():
            name = "cuda"
        else:
            name = "cpu"
    device = torch.device(name)
    if device.type == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("CUDA was requested but is unavailable")
    if device.type == "mps":
        if not torch.backends.mps.is_available():
            raise RuntimeError("MPS was requested but is unavailable")
        raise RuntimeError(
            "TorchSim 0.6.1's default nvalchemiops/Warp neighbor list does not "
            "support PyTorch MPS; use CPU here or CUDA on a supported host")
    return device


def _fcc_slab(*, device: torch.device, dtype: torch.dtype):
    """Return one finite-z, periodic-xy argon slab plus a normal projectile."""
    base = torch.tensor(
        ((0.0, 0.0, 0.0), (0.0, 0.5, 0.5),
         (0.5, 0.0, 0.5), (0.5, 0.5, 0.0)),
        device=device, dtype=dtype)
    repeats = (2, 2, 2)
    positions = []
    for i in range(repeats[0]):
        for j in range(repeats[1]):
            for k in range(repeats[2]):
                offset = torch.tensor((i, j, k), device=device, dtype=dtype)
                positions.extend(base + offset)
    slab = torch.stack(positions) * ARGON_LATTICE_ANGSTROM
    slab[:, 2] += 6.0
    initial_surface_z = float(slab[:, 2].max())

    lateral_extent = torch.tensor(
        (repeats[0] * ARGON_LATTICE_ANGSTROM,
         repeats[1] * ARGON_LATTICE_ANGSTROM),
        device=device, dtype=dtype)
    projectile = torch.tensor(
        ((0.5 * lateral_extent[0], 0.5 * lateral_extent[1],
          initial_surface_z + 4.5),),
        device=device, dtype=dtype)
    all_positions = torch.cat((slab, projectile), dim=0)
    n_atoms = all_positions.shape[0]
    masses = torch.full(
        (n_atoms,), ARGON_MASS_AMU, device=device, dtype=dtype)
    atomic_numbers = torch.full(
        (n_atoms,), ARGON_ATOMIC_NUMBER, device=device, dtype=torch.int64)
    cell = torch.diag(torch.tensor(
        (float(lateral_extent[0]), float(lateral_extent[1]), 50.0),
        device=device, dtype=dtype))
    # Anchor the lowest crystallographic plane. All other target atoms and the
    # projectile evolve under the same deterministic Hamiltonian.
    fixed = torch.where(slab[:, 2] <= slab[:, 2].min() + 1.0e-8)[0]
    state = ts.SimState(
        positions=all_positions,
        masses=masses,
        cell=cell,
        atomic_numbers=atomic_numbers,
        pbc=[True, True, False],
    )
    state.constraints = [FixAtoms(atom_idx=fixed)]
    return state, initial_surface_z


def _differentiability_probe(device: torch.device, dtype: torch.dtype):
    """Separate static differentiability from gradients through NVE dynamics."""
    positions = torch.tensor(
        ((5.0, 5.0, 5.0),
         (5.0 + 2.0 ** (1.0 / 6.0) * ARGON_SIGMA_ANGSTROM, 5.0, 5.0)),
        device=device, dtype=dtype, requires_grad=True)
    state = ts.SimState(
        positions=positions,
        masses=torch.full((2,), ARGON_MASS_AMU, device=device, dtype=dtype),
        cell=torch.eye(3, device=device, dtype=dtype) * 20.0,
        atomic_numbers=torch.full(
            (2,), ARGON_ATOMIC_NUMBER, device=device, dtype=torch.int64),
        pbc=False)
    model = LennardJonesModel(
        sigma=ARGON_SIGMA_ANGSTROM,
        epsilon=ARGON_EPSILON_EV,
        cutoff=2.5 * ARGON_SIGMA_ANGSTROM,
        device=device, dtype=dtype, compute_forces=True, retain_graph=True)
    result = model(state)
    energy_gradient = torch.autograd.grad(result["energy"].sum(), positions)[0]
    error = torch.max(torch.abs(energy_gradient + result["forces"]))
    static_result = {
        "energy_eV": float(result["energy"].detach().cpu()),
        "force_energy_gradient_max_abs_error_eV_per_angstrom": float(
            error.detach().cpu()),
    }

    def final_separation(sigma):
        trajectory_positions = torch.tensor(
            ((5.0, 5.0, 5.0), (9.5, 5.0, 5.0)),
            device=device, dtype=dtype)
        trajectory_state = ts.SimState(
            positions=trajectory_positions,
            masses=torch.full((2,), ARGON_MASS_AMU, device=device, dtype=dtype),
            cell=torch.eye(3, device=device, dtype=dtype) * 20.0,
            atomic_numbers=torch.full(
                (2,), ARGON_ATOMIC_NUMBER, device=device, dtype=torch.int64),
            pbc=False)
        trajectory_model = LennardJonesModel(
            sigma=sigma, epsilon=ARGON_EPSILON_EV,
            cutoff=2.5 * ARGON_SIGMA_ANGSTROM,
            device=device, dtype=dtype, compute_forces=True, retain_graph=True)
        trajectory_state = ts.nve_init(
            state=trajectory_state, model=trajectory_model,
            kT=torch.zeros(1, device=device, dtype=dtype))
        dt = torch.tensor(0.001 * MetalUnits.time, device=device, dtype=dtype)
        for _ in range(50):
            trajectory_state = ts.nve_step(
                state=trajectory_state, model=trajectory_model, dt=dt)
        return torch.linalg.vector_norm(
            trajectory_state.positions[1] - trajectory_state.positions[0])

    sigma = torch.tensor(
        ARGON_SIGMA_ANGSTROM, device=device, dtype=dtype, requires_grad=True)
    separation = final_separation(sigma)
    if separation.requires_grad:
        gradient = torch.autograd.grad(
            separation, sigma, allow_unused=True)[0]
    else:
        gradient = None
    delta = 1.0e-3
    upper = final_separation(ARGON_SIGMA_ANGSTROM + delta)
    lower = final_separation(ARGON_SIGMA_ANGSTROM - delta)
    finite_difference = (upper - lower) / (2.0 * delta)
    static_result.update({
        "nve_parameter_gradient_available": gradient is not None,
        "nve_final_separation_gradient_wrt_sigma_autodiff": (
            None if gradient is None else float(gradient.detach().cpu())),
        "nve_final_separation_gradient_wrt_sigma_finite_difference": float(
            finite_difference.detach().cpu()),
        "nve_parameter_gradient_gate_passed": bool(
            gradient is not None
            and torch.isclose(
                gradient.detach(), finite_difference.detach(), rtol=1.0e-3,
                atol=1.0e-6)),
    })
    return static_result


def _classify(system, initial_surface_z: float):
    positions = system.positions.detach().cpu()
    momenta = system.momenta.detach().cpu()
    target_z = positions[:-1, 2]
    projectile_z = float(positions[-1, 2])
    projectile_pz = float(momenta[-1, 2])
    ejected = int(torch.count_nonzero(target_z > initial_surface_z + 2.0))
    if projectile_z > initial_surface_z + 3.0 and projectile_pz > 0.0:
        outcome = "reflected"
    elif projectile_z <= initial_surface_z + 2.0:
        outcome = "trapped_or_implanted"
    else:
        outcome = "undetermined"
    return {
        "projectile_outcome": outcome,
        "projectile_final_z_angstrom": projectile_z,
        "projectile_final_pz_sqrt_amu_eV": projectile_pz,
        "ejected_target_count": ejected,
        "maximum_target_z_angstrom": float(target_z.max()),
    }


def run(energies_eV, *, steps: int, timestep_ps: float, device: torch.device):
    dtype = torch.float64
    torch.manual_seed(0)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(0)
    torch.use_deterministic_algorithms(True)

    states = []
    surface_heights = []
    for _energy in energies_eV:
        state, surface_z = _fcc_slab(device=device, dtype=dtype)
        states.append(state)
        surface_heights.append(surface_z)
    state = ts.concatenate_states(states)
    model = LennardJonesModel(
        sigma=ARGON_SIGMA_ANGSTROM,
        epsilon=ARGON_EPSILON_EV,
        cutoff=2.5 * ARGON_SIGMA_ANGSTROM,
        device=device, dtype=dtype, compute_forces=True)
    state = ts.nve_init(
        state=state, model=model,
        kT=torch.zeros(len(energies_eV), device=device, dtype=dtype))
    counts = torch.bincount(state.system_idx)
    projectile_indices = torch.cumsum(counts, dim=0) - 1
    prescribed_momenta = torch.zeros_like(state.momenta)
    prescribed_momenta[projectile_indices, 2] = -torch.sqrt(
        2.0 * ARGON_MASS_AMU
        * torch.tensor(energies_eV, device=device, dtype=dtype))
    state.set_constrained_momenta(prescribed_momenta)
    initial_energy = state.energy.detach().cpu()
    initial_kinetic = ts.calc_kinetic_energy(
        masses=state.masses, momenta=state.momenta, system_idx=state.system_idx
    ).detach().cpu()
    initial_total = initial_energy + initial_kinetic
    dt = torch.tensor(timestep_ps * MetalUnits.time, device=device, dtype=dtype)

    started = time.perf_counter()
    for _ in range(steps):
        state = ts.nve_step(state=state, model=model, dt=dt)
    if device.type == "cuda":
        torch.cuda.synchronize()
    elif device.type == "mps":
        torch.mps.synchronize()
    elapsed_s = time.perf_counter() - started

    final_kinetic = ts.calc_kinetic_energy(
        masses=state.masses, momenta=state.momenta, system_idx=state.system_idx
    ).detach().cpu()
    final_total = state.energy.detach().cpu() + final_kinetic
    split = state.split()
    cases = []
    for index, (energy, system, surface_z) in enumerate(
            zip(energies_eV, split, surface_heights, strict=True)):
        drift = float(final_total[index] - initial_total[index])
        scale = max(abs(float(initial_total[index])), 1.0e-12)
        cases.append({
            "incident_energy_eV": energy,
            "initial_total_energy_eV": float(initial_total[index]),
            "final_total_energy_eV": float(final_total[index]),
            "energy_drift_eV": drift,
            "relative_energy_drift": drift / scale,
            "energy_drift_over_incident_energy": drift / energy,
            **_classify(system, surface_z),
        })

    return {
        "schema_version": 1,
        "experiment": "torchsim_deterministic_batched_impact_feasibility",
        "predictive_physics": False,
        "nonpredictive_reason": (
            "Lennard-Jones argon is a software/contract probe and cannot represent "
            "reactive semiconductor plasma-surface chemistry."),
        "software": {
            "torch_version": torch.__version__,
            "torchsim_version": ts.__version__,
            "python_version": platform.python_version(),
            "device": str(device),
            "dtype": str(dtype).removeprefix("torch."),
        },
        "configuration": {
            "energies_eV": list(energies_eV),
            "steps": steps,
            "timestep_ps": timestep_ps,
            "integrator": "NVE velocity Verlet",
            "potential": "analytic Lennard-Jones Ar",
            "sigma_angstrom": ARGON_SIGMA_ANGSTROM,
            "epsilon_eV": ARGON_EPSILON_EV,
            "seed": 0,
            "deterministic_algorithms": True,
            "batched_systems": len(energies_eV),
        },
        "performance": {
            "elapsed_s": elapsed_s,
            "trajectory_steps_per_s": steps * len(energies_eV) / elapsed_s,
        },
        "differentiability": _differentiability_probe(device, dtype),
        "cases": cases,
        "interaction_table": {
            "schema_version": 1,
            "material": "Ar_Lennard_Jones_slab_NONPREDICTIVE",
            "incident_species": ["Ar_projectile"],
            "axes": [{
                "name": "incident_energy",
                "values": list(energies_eV),
                "unit": "eV",
                "interpolation": "linear",
            }],
            "outputs": {
                "ejected_target_count_per_projectile": [
                    case["ejected_target_count"] for case in cases],
            },
            "output_units": {
                "ejected_target_count_per_projectile": "Ar/projectile",
            },
            "provenance": {
                "source": "petch TorchSim 0.6.1 feasibility experiment",
                "evidence_type": "toy_atomistic_model_nonpredictive",
                "supports_prediction_within_declared_domain": False,
                "potential": "analytic Lennard-Jones Ar",
                "purpose": "software_and_exchange_contract_only",
            },
            "standard_uncertainty": {},
            "bounds": {
                "ejected_target_count_per_projectile": [0.0, None],
            },
            "conservation_groups": {},
        },
    }


def _parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--device", choices=("auto", "cpu", "mps", "cuda"), default="cpu")
    parser.add_argument("--steps", type=int, default=1500)
    parser.add_argument("--timestep-ps", type=float, default=0.001)
    parser.add_argument("--energies-eV", nargs="+", type=float, default=DEFAULT_ENERGIES_EV)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main():
    arguments = _parse_args()
    energies = tuple(sorted(set(arguments.energies_eV)))
    if (len(energies) < 2 or any(not math.isfinite(value) or value <= 0.0 for value in energies)
            or arguments.steps <= 0 or arguments.timestep_ps <= 0.0):
        raise ValueError("provide at least two increasing positive energies and positive integration controls")
    payload = run(
        energies, steps=arguments.steps, timestep_ps=arguments.timestep_ps,
        device=_device(arguments.device))
    encoded = json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n"
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(encoded, encoding="utf-8")
    print(json.dumps({
        "output": str(arguments.output),
        "sha256": sha256(encoded.encode()).hexdigest(),
        "device": payload["software"]["device"],
        "elapsed_s": payload["performance"]["elapsed_s"],
        "cases": payload["cases"],
        "differentiability": payload["differentiability"],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
