"""Atom-counted Si/Cl2/Ar+ ALE depth transfer from DeepMD to experiment.

This module deliberately does not use the overloaded ``monolayer`` label as a
conversion between atomistic and experimental fluence.  The released DeepMD
ALE trajectory is converted through its physical cell area and explicit atom
count.  The experiment is converted through its measured positive-ion flux and
bombardment duration.

The first 1000 normal-incidence Ar impacts use the released, chlorinated DeepMD
cycle.  Any remaining experimental fluence is a bare-Si physical-sputter tail
using the independent released DeepMD sputter table.  This is a no-depth-fit
cross-source transfer.  It is not a blind prediction: the atomistic work and
experiment share a research lineage, and the experimental mean energy was
inferred rather than measured as a species-resolved IEAD.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .surface_interaction_table import SurfaceInteractionTable


DIAMOND_SI_LATTICE_CONSTANT_ANGSTROM = 5.43
DIAMOND_SI_ATOMS_PER_CONVENTIONAL_CELL = 8
DEEPMD_CELL_REPEATS_XY = 6
DEEPMD_CELL_LENGTH_ANGSTROM = (
    DEEPMD_CELL_REPEATS_XY * DIAMOND_SI_LATTICE_CONSTANT_ANGSTROM)
DEEPMD_SI_ATOMS_PER_MATERIAL_ML = 72
DEEPMD_AR_IMPACTS_PER_ALE_CYCLE = 1000
DEEPMD_CL2_IMPACTS_PER_ALE_CYCLE = 2255


@dataclass(frozen=True)
class VellaHaoAleBoundary:
    """Facility-conditioned experimental boundary reported by Vella et al."""

    positive_ion_flux_cm2_s: float = 3.7e16
    ion_bombardment_duration_s: float = 3.0

    def __post_init__(self):
        if (not np.isfinite(self.positive_ion_flux_cm2_s)
                or self.positive_ion_flux_cm2_s <= 0.0
                or not np.isfinite(self.ion_bombardment_duration_s)
                or self.ion_bombardment_duration_s <= 0.0):
            raise ValueError("ALE boundary flux and duration must be positive and finite")

    @property
    def ion_fluence_cm2(self):
        return self.positive_ion_flux_cm2_s * self.ion_bombardment_duration_s


@dataclass(frozen=True)
class SiClAleAbsoluteDepth:
    """Dimensional ledger for one no-fit absolute-depth prediction."""

    mean_ion_energy_eV: float
    steady_cycle_si_material_ml: float
    simulated_cell_area_cm2: float
    simulated_ar_fluence_cm2: float
    experimental_positive_ion_fluence_cm2: float
    physical_sputter_tail_fluence_cm2: float
    physical_sputter_yield_si_per_ar: float
    source_reported_sputter_yield_uncertainty: float
    chlorinated_transient_si_atoms_cm2: float
    physical_sputter_tail_si_atoms_cm2: float
    total_removed_si_atoms_cm2: float
    silicon_atomic_density_cm3: float
    chlorinated_transient_depth_nm: float
    physical_sputter_tail_depth_nm: float
    total_depth_nm: float
    source_reported_tail_depth_uncertainty_nm: float
    ale_cycle_table_fingerprint: str
    sputter_table_fingerprint: str

    @property
    def dimensional_atom_balance_residual_cm2(self):
        return (
            self.total_removed_si_atoms_cm2
            - self.chlorinated_transient_si_atoms_cm2
            - self.physical_sputter_tail_si_atoms_cm2
        )


def silicon_atomic_density_cm3(
        lattice_constant_angstrom=DIAMOND_SI_LATTICE_CONSTANT_ANGSTROM):
    """Diamond-cubic Si number density from eight atoms per conventional cell."""
    lattice_constant_cm = float(lattice_constant_angstrom) * 1.0e-8
    if not np.isfinite(lattice_constant_cm) or lattice_constant_cm <= 0.0:
        raise ValueError("lattice constant must be positive and finite")
    return DIAMOND_SI_ATOMS_PER_CONVENTIONAL_CELL / lattice_constant_cm ** 3


def deepmd_cell_area_cm2():
    return (DEEPMD_CELL_LENGTH_ANGSTROM * 1.0e-8) ** 2


def steady_cycle_removal_material_ml(
        ale_cycle_table: SurfaceInteractionTable,
        ion_energy_eV: float,
        *,
        tail_cycles: int = 3,
):
    """Mean of completed-cycle increments, excluding the first transient cycle."""
    if tail_cycles < 1 or tail_cycles > 3:
        raise ValueError("the four-cycle source supports a 1--3 cycle tail mean")
    cycles = np.asarray([1.0, 2.0, 3.0, 4.0])
    evaluated = ale_cycle_table.evaluate({
        "completed_cycle": cycles,
        "ion_energy": np.full(cycles.shape, float(ion_energy_eV)),
    })
    cumulative = evaluated.values["cumulative_si_etched_material_ml"]
    increments = np.diff(cumulative)
    if np.any(increments < 0.0):
        raise RuntimeError("cumulative DeepMD ALE removal decreased between cycles")
    return float(np.mean(increments[-tail_cycles:]))


def predict_vella_hao_ale_depth(
        ion_energy_eV: float,
        ale_cycle_table: SurfaceInteractionTable,
        sputter_table: SurfaceInteractionTable,
        boundary: VellaHaoAleBoundary = VellaHaoAleBoundary(),
):
    """Predict absolute EPC without calibrating to any measured depth.

    The supported energy domain is the intersection of the released ALE
    trajectory (40--100 eV) and physical-sputter table (50--200 eV).
    Table extrapolation remains refused.
    """
    energy = float(ion_energy_eV)
    if not np.isfinite(energy):
        raise ValueError("ion energy must be finite")
    cycle_removal_ml = steady_cycle_removal_material_ml(
        ale_cycle_table, energy)
    sputter = sputter_table.evaluate({"ion_energy": energy})
    sputter_yield = float(sputter.values["physical_sputter_yield"])
    sputter_uncertainty = float(
        sputter.standard_uncertainty["physical_sputter_yield"])

    cell_area = deepmd_cell_area_cm2()
    simulated_fluence = DEEPMD_AR_IMPACTS_PER_ALE_CYCLE / cell_area
    experimental_fluence = boundary.ion_fluence_cm2
    tail_fluence = experimental_fluence - simulated_fluence
    if tail_fluence < 0.0:
        raise ValueError(
            "experimental ion fluence is smaller than the atomistic ALE sequence")

    transient_atoms = (
        cycle_removal_ml
        * DEEPMD_SI_ATOMS_PER_MATERIAL_ML
        / cell_area
    )
    tail_atoms = tail_fluence * sputter_yield
    total_atoms = transient_atoms + tail_atoms
    density = silicon_atomic_density_cm3()
    cm_to_nm = 1.0e7
    transient_depth = transient_atoms / density * cm_to_nm
    tail_depth = tail_atoms / density * cm_to_nm
    tail_depth_uncertainty = (
        tail_fluence * sputter_uncertainty / density * cm_to_nm)
    return SiClAleAbsoluteDepth(
        mean_ion_energy_eV=energy,
        steady_cycle_si_material_ml=cycle_removal_ml,
        simulated_cell_area_cm2=cell_area,
        simulated_ar_fluence_cm2=simulated_fluence,
        experimental_positive_ion_fluence_cm2=experimental_fluence,
        physical_sputter_tail_fluence_cm2=tail_fluence,
        physical_sputter_yield_si_per_ar=sputter_yield,
        source_reported_sputter_yield_uncertainty=sputter_uncertainty,
        chlorinated_transient_si_atoms_cm2=transient_atoms,
        physical_sputter_tail_si_atoms_cm2=tail_atoms,
        total_removed_si_atoms_cm2=total_atoms,
        silicon_atomic_density_cm3=density,
        chlorinated_transient_depth_nm=transient_depth,
        physical_sputter_tail_depth_nm=tail_depth,
        total_depth_nm=transient_depth + tail_depth,
        source_reported_tail_depth_uncertainty_nm=tail_depth_uncertainty,
        ale_cycle_table_fingerprint=ale_cycle_table.fingerprint,
        sputter_table_fingerprint=sputter_table.fingerprint,
    )


def printed_rom_chlorine_creation_per_ar(
        theta_top: float,
        theta_mixed: float,
        sicl2_yield: float,
):
    """Cl atoms created per Ar by the published transient-ROM output equation.

    The coverage equations remove two Cl atoms with a theta-squared SiCl2
    channel, but the printed product equation emits SiCl2 linearly in theta.
    Their difference is an elemental source except at zero/full coverage.
    """
    top = float(theta_top)
    mixed = float(theta_mixed)
    yield_value = float(sicl2_yield)
    if (not np.isfinite(top) or not np.isfinite(mixed)
            or not np.isfinite(yield_value)
            or top < 0.0 or top > 1.0
            or mixed < 0.0 or mixed > 1.0
            or yield_value < 0.0):
        raise ValueError("coverages must be in [0,1] and yield must be nonnegative")
    return 2.0 * yield_value * (
        top + mixed - top ** 2 - mixed ** 2)
