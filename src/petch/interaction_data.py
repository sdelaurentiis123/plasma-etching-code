"""Checksum-gated physical surface-interaction data from public primary sources."""
from __future__ import annotations

import csv
from dataclasses import dataclass
from hashlib import sha256
from pathlib import Path

import numpy as np

from .surface_interaction_table import InteractionAxis, SurfaceInteractionTable


KOUNIS_MELAS_2024_SHA256 = {
    "Sputtering.csv": "80ae627c1cec67258496ee7d22bd130817b678c1fd3288d5141436fcf374ee3c",
    "RIE.csv": "7cc634ae1218ba12d1e30ba7e6b4aefc0f4f0cc6de04ced8120115a60786cc77",
    "Products.csv": "79a7cd3a2618a3fc3d65946d2db5247870d428b58270f78b0ffe46b5116bd9bf",
    "ALE_cycle_endpoints.csv": (
        "4867fbc5c671ac7c9c1efba8371f8d8e3009ec21b77f7f41f10a486414d94a36"),
}
KOUNIS_MELAS_2024_ARCHIVE_SHA256 = (
    "4c9fa0b9268ac314da77b1012906dff4e45c5af79afd7ea674b26ace48e0f269")
KOUNIS_MELAS_2024_ALE_TRAJECTORY_SHA256 = (
    "e76b9895a8191c923a492130452a73a4c929589b5260d9847e1d750fc22be884")


@dataclass(frozen=True)
class KounisMelas2024Tables:
    sputtering: SurfaceInteractionTable
    reactive_ion_etch: SurfaceInteractionTable
    ale_products: SurfaceInteractionTable
    ale_cycles: SurfaceInteractionTable


def _verified_rows(path, expected_fields, verify_checksum):
    path = Path(path); payload = path.read_bytes()
    expected_hash = KOUNIS_MELAS_2024_SHA256.get(path.name)
    if expected_hash is None:
        raise ValueError(f"unrecognized interaction-data file: {path.name}")
    if verify_checksum and sha256(payload).hexdigest() != expected_hash:
        raise ValueError(f"checksum mismatch for interaction data: {path}")
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames != expected_fields:
            raise ValueError(f"unexpected interaction-data schema: {reader.fieldnames}")
        rows = list(reader)
    if not rows:
        raise ValueError(f"interaction data is empty: {path}")
    return rows


def _provenance(table_name, **conditions):
    return {
        "source": (
            "Kounis-Melas et al., Data from Deep Potential Molecular Dynamics Simulations "
            "of Low-Temperature Plasma-Surface Interactions"),
        "evidence_type": "DeepMD_molecular_dynamics",
        "supports_prediction_within_declared_domain": True,
        "dataset_doi": "10.34770/rjv6-2w31",
        "paper_doi": "10.1116/6.0004027",
        "osti_dataset_id": "2589032",
        "license": "CC-BY-4.0",
        "archive_sha256": KOUNIS_MELAS_2024_ARCHIVE_SHA256,
        "source_table": table_name,
        "source_table_sha256": KOUNIS_MELAS_2024_SHA256[table_name],
        "conditions": conditions,
    }


def load_kounis_melas_2024_tables(directory, *, verify_checksum=True):
    """Load the archived Si-Cl2-Ar+ DeepMD summaries into general interaction tables.

    These are MD outputs and retain that evidence label. They are a sourced Si chemistry extension,
    not a substitute for SiO2/fluorocarbon surface data and not experimental validation by themselves.
    """
    directory = Path(directory)
    sputter_rows = _verified_rows(
        directory / "Sputtering.csv",
        ["Energy (eV)", "Yield", "Yield ±", "Thickness (Å)", "Thickness ±"],
        verify_checksum)
    rie_rows = _verified_rows(
        directory / "RIE.csv", ["Flux Ratio", "Yield", "Yield ±"], verify_checksum)
    product_rows = _verified_rows(
        directory / "Products.csv",
        ["Ion Dosage (cm^-2) × 10^15", "Si Yield", "Si Yield ±", "SiCl Yield",
         "SiCl Yield ±", "SiCl2 Yield", "SiCl2 Yield ±", "Cl Yield", "Cl Yield ±"],
        verify_checksum)
    ale_cycle_rows = _verified_rows(
        directory / "ALE_cycle_endpoints.csv",
        ["Completed cycle", "Ion energy (eV)",
         "Cumulative Si etched (material ML)",
         "Post-bombardment Cl (material ML)", "Source row index",
         "Source cycle coordinate"],
        verify_checksum)

    energy = np.asarray([float(row["Energy (eV)"]) for row in sputter_rows])
    sputtering = SurfaceInteractionTable(
        material="Si(100)", incident_species=("Ar+",),
        axes=(InteractionAxis("ion_energy", energy, "eV"),),
        outputs={
            "physical_sputter_yield": [float(row["Yield"]) for row in sputter_rows],
            "amorphous_layer_thickness": [
                float(row["Thickness (Å)"]) for row in sputter_rows],
        },
        output_units={
            "physical_sputter_yield": "Si/Ar+",
            "amorphous_layer_thickness": "angstrom",
        },
        standard_uncertainty={
            "physical_sputter_yield": [float(row["Yield ±"]) for row in sputter_rows],
            "amorphous_layer_thickness": [
                float(row["Thickness ±"]) for row in sputter_rows],
        },
        bounds={
            "physical_sputter_yield": (0.0, None),
            "amorphous_layer_thickness": (0.0, None),
        },
        provenance=_provenance(
            "Sputtering.csv", incidence_angle_deg=0.0, substrate_temperature_k=298.0))

    flux_ratio = np.asarray([float(row["Flux Ratio"]) for row in rie_rows])
    reactive_ion_etch = SurfaceInteractionTable(
        material="Si(100)", incident_species=("Ar+", "Cl2"),
        axes=(InteractionAxis(
            "cl2_to_ar_flux_ratio", flux_ratio, "1", interpolation="log"),),
        outputs={"reactive_etch_yield": [float(row["Yield"]) for row in rie_rows]},
        output_units={"reactive_etch_yield": "Si/Ar+"},
        standard_uncertainty={
            "reactive_etch_yield": [float(row["Yield ±"]) for row in rie_rows]},
        bounds={"reactive_etch_yield": (0.0, None)},
        provenance=_provenance(
            "RIE.csv", ar_ion_energy_eV=100.0, incidence_angle_deg=0.0,
            substrate_temperature_k=298.0))

    dosage = np.asarray([
        float(row["Ion Dosage (cm^-2) × 10^15"]) for row in product_rows])
    product_names = {
        "si_yield": ("Si Yield", "Si Yield ±", "Si/Ar+"),
        "sicl_yield": ("SiCl Yield", "SiCl Yield ±", "SiCl/Ar+"),
        "sicl2_yield": ("SiCl2 Yield", "SiCl2 Yield ±", "SiCl2/Ar+"),
        "cl_yield": ("Cl Yield", "Cl Yield ±", "Cl/Ar+"),
    }
    ale_products = SurfaceInteractionTable(
        material="Si(100)", incident_species=("Ar+", "Cl2"),
        axes=(InteractionAxis("ar_ion_dosage", dosage, "1e15 cm^-2"),),
        outputs={
            name: [float(row[column]) for row in product_rows]
            for name, (column, _, _) in product_names.items()},
        output_units={name: unit for name, (_, _, unit) in product_names.items()},
        standard_uncertainty={
            name: [float(row[uncertainty]) for row in product_rows]
            for name, (_, uncertainty, _) in product_names.items()},
        bounds={name: (0.0, None) for name in product_names},
        provenance=_provenance(
            "Products.csv", ar_ion_energy_eV=215.0, incidence_angle_deg=0.0,
            substrate_temperature_k=298.0,
            source_figure="Figures 13-14 and Supplementary Figures S1-S2",
            note=(
                "The separate Figure 12 ALE morphology sequence uses 80 eV; "
                "the Products.csv points are the 215 eV product-yield sequence.")))

    completed_cycles = np.asarray([1.0, 2.0, 3.0, 4.0])
    ale_energies = np.asarray([40.0, 60.0, 80.0, 100.0])
    keyed_rows = {}
    expected_source_rows = {1: 3254, 2: 6509, 3: 9764, 4: 13019}
    for row in ale_cycle_rows:
        cycle = int(row["Completed cycle"])
        energy_node = float(row["Ion energy (eV)"])
        key = (cycle, energy_node)
        if key in keyed_rows:
            raise ValueError(f"duplicate ALE cycle endpoint: {key}")
        if (cycle not in expected_source_rows
                or int(row["Source row index"]) != expected_source_rows[cycle]):
            raise ValueError(f"unexpected ALE source endpoint: {key}")
        keyed_rows[key] = row
    expected_keys = {
        (int(cycle), float(energy_node))
        for cycle in completed_cycles
        for energy_node in ale_energies
    }
    if set(keyed_rows) != expected_keys:
        raise ValueError("ALE cycle endpoint grid is incomplete or contains an unknown node")

    def ale_grid(column):
        return np.asarray([
            [
                float(keyed_rows[(int(cycle), float(energy_node))][column])
                for energy_node in ale_energies
            ]
            for cycle in completed_cycles
        ])

    ale_cycles = SurfaceInteractionTable(
        material="Si(100)",
        incident_species=("Ar+", "Cl2"),
        axes=(
            InteractionAxis("completed_cycle", completed_cycles, "cycle"),
            InteractionAxis("ion_energy", ale_energies, "eV"),
        ),
        outputs={
            "cumulative_si_etched_material_ml": ale_grid(
                "Cumulative Si etched (material ML)"),
            "post_bombardment_cl_material_ml": ale_grid(
                "Post-bombardment Cl (material ML)"),
        },
        output_units={
            "cumulative_si_etched_material_ml": (
                "72 Si atoms per 32.58 angstrom square cell"),
            "post_bombardment_cl_material_ml": (
                "72 Cl atoms per 32.58 angstrom square cell"),
        },
        bounds={
            "cumulative_si_etched_material_ml": (0.0, None),
            "post_bombardment_cl_material_ml": (0.0, None),
        },
        provenance=_provenance(
            "ALE_cycle_endpoints.csv",
            parent_source_table="DeepMDData/ALE/ALE.csv",
            parent_source_table_sha256=KOUNIS_MELAS_2024_ALE_TRAJECTORY_SHA256,
            reduction=(
                "For each completed cycle, exact final archived row with "
                "Cycle coordinate below the next integer"),
            source_row_indices=[3254, 6509, 9764, 13019],
            silicon_atoms_per_material_ml=72,
            chlorine_atoms_per_material_ml=72,
            simulation_cell_x_angstrom=32.58,
            simulation_cell_y_angstrom=32.58,
            cl2_impacts_per_cycle=2255,
            ar_impacts_per_cycle=1000,
            incidence_angle_deg=0.0,
            substrate_temperature_k=298.0,
            note=(
                "Material ML is an atom count in the simulation cell. It is "
                "not the source's separate 1e15 impacts/cm2 fluence shorthand.")))
    return KounisMelas2024Tables(
        sputtering, reactive_ion_etch, ale_products, ale_cycles)
