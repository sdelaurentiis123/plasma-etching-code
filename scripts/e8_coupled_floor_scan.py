"""E8 floor composition measured THROUGH the coupled radiosity solve.

The first E8 pass (RESULTS_E8_THERMALIZED_RETURN_2026-08-05.md) measured the
thermalized weight at the face where the cascade dropped it -- the raw gather
diagnostic -- and concluded the weight was "stranded on the sidewalls".  That
measurement could not see re-emission: the gather deposit is the *source* term,
not the delivered flux.  petch already solves multi-bounce diffuse re-emission
for plasma-sourced neutrals (``solve_diffuse_neutral_radiosity_3d``, H = D +
B(1-s)H), and the gather writes the E8 weight into the same per-species neutral
ledger that becomes ``D``.  So once ``thermalized_radical_return`` is plumbed to
the gather, the thermalized radicals diffuse at their own published sticking
exactly as Huang describes:

    "After losing energy through several collisions with the sidewalls and etch
    front, these energetic species become thermal CFx and CxFy radicals, which
    can passivate the oxide surface or deposit as polymer... As the AR increases
    to greater than 10, the neutralized and thermalized CFx+ and CxFy+ ions
    become the main source (> 95%) of radicals reaching the etch front."
    -- huang_thesis.txt L5714-5727 (research_sources/thesis_extracts/)

The solve is linear in its source term, so the E8-delivered floor flux is
measured exactly by differencing two solves on the identical operator (same
frozen geometry, same reaction probabilities): E8 off, then E8 on.

The reactive fraction remains a DECLARED input.  Krueger publishes an aggregate
positive-ion flux with a combined IEAD, so the CFx+/Ar+ split does not exist for
his reactor; it is swept over its physical band [0, 1], never fitted.
"""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from cascade_funnelling_scan import (  # noqa: E402
    make_straight_trench_geometry_3d, floor_face_mask, OPEN_WIDTH_UM,
    MASK_THICKNESS_UM, CELL_WIDTH_UM, CELL_LENGTH_UM, HEADSPACE_UM,
)
from mouth_equilibrium_probe import _boundary_and_mechanism  # noqa: E402
from petch.feature_step_3d import advance_feature_step_3d  # noqa: E402

# Fluorocarbon radicals in Krueger Table 6.1; the thermalized ion partners are
# returned as CF2, the dominant reactive fluorocarbon radical in that table.
# Huang sec. 6.4.3 excludes non-fluorocarbon partners explicitly ("non-reactive
# species... diffuse out of the feature with no surface reactions").
E8_TARGET_SPECIES = "CF2"


def coupled_gather(geometry, *, thermalized_radical_return=None,
                   transport_device="cpu", ion_azimuthal_order=16,
                   face_quadrature_points=3):
    """Pilot-configuration step: ballistic gather + deterministic radiosity."""
    boundary, mechanism, role, realized, source_z = _boundary_and_mechanism(
        geometry, ion_azimuthal_order=ion_azimuthal_order)
    radiosity = {
        "form_factor_backend": "deterministic_extruded_2d",
        "periodic_lateral": True,
        "domain_size": tuple(float(v) for v in realized),
        "relative_tolerance": 1e-4,
        "maximum_iterations": 200,
        "deterministic_extruded_options": {
            "exchange_method": "analytic_occlusion",
            "exchange_relative_tolerance": 1e-4,
            "maximum_refinement_level": 24,
        },
    }
    return advance_feature_step_3d(
        geometry, boundary, role, mechanism, etchable_material_ids=(1, 2),
        duration_s=0.0, source_bounds=(0.0, realized[0], 0.0, realized[1]),
        source_z=source_z, n_position=1, seed=1, cfl_number=0.25,
        reinitialize=True, reinitialization_method="cr2",
        profile_periodic_lateral=True, transport_device=str(transport_device),
        neutral_radiosity_options=radiosity,
        ballistic_transport="face_gather", grazing_ion_reflection={},
        thermalized_radical_return=thermalized_radical_return,
        ballistic_face_quadrature_points=int(face_quadrature_points),
        topology_change_policy="continue_gas_cavity",
        surface_state_remap_backend="common_refinement")


def _floor_flux(result, mask, species):
    """Area-weighted incident flux of one species on the floor band."""
    active = np.asarray(result.active_face_index, dtype=int)
    areas = np.asarray(result.active_face_area, dtype=float)
    flux = np.asarray(
        result.transport.surface_fluxes.neutral_flux_m2_s[species], dtype=float)
    weight = areas[mask]
    return float((flux[active][mask] * weight).sum() / max(weight.sum(), 1e-300))


def scan(*, etched_depth, dx, mask_thickness, fractions, label,
         transport_device="cpu"):
    geometry, floor_z = make_straight_trench_geometry_3d(
        etched_depth=etched_depth, dx=dx, mask_thickness=mask_thickness)
    base = coupled_gather(geometry, transport_device=transport_device)
    mask = floor_face_mask(base, floor_z_um=floor_z)
    plasma = _floor_flux(base, mask, E8_TARGET_SPECIES)
    records = []
    for fraction in fractions:
        if fraction == 0.0:
            delivered, total = 0.0, plasma
        else:
            withe8 = coupled_gather(
                geometry, thermalized_radical_return={
                    E8_TARGET_SPECIES: float(fraction)},
                transport_device=transport_device)
            total = _floor_flux(withe8, mask, E8_TARGET_SPECIES)
            delivered = total - plasma
        records.append(dict(
            label=label, oxide_ar=float(etched_depth / OPEN_WIDTH_UM),
            mask_um=float(mask_thickness), fc_fraction=float(fraction),
            plasma_floor_flux=plasma, e8_delivered_floor_flux=delivered,
            total_floor_flux=total,
            e8_share=delivered / max(total, 1e-300),
            floor_faces=int(mask.sum())))
        print(json.dumps(records[-1]), flush=True)
    return records


if __name__ == "__main__":
    out = ROOT / "results/curated/e8_thermalized_return"
    out.mkdir(parents=True, exist_ok=True)
    dx = float(sys.argv[1]) if len(sys.argv) > 1 else 0.01
    fractions = (0.0, 0.3, 1.0)
    records = []
    # Krueger cell: 0.85 um mask over a 0.09 um opening (mask-dominated).
    records += scan(etched_depth=12.0 * OPEN_WIDTH_UM, dx=dx,
                    mask_thickness=MASK_THICKNESS_UM, fractions=fractions,
                    label="krueger_cell_ar12")
    # Huang-like: thin mask over the same opening, so the oxide aspect ratio --
    # not the mask stack -- sets neutral transmission (his >95% statement is for
    # a feature where plasma-neutral delivery to the front is far lower).
    records += scan(etched_depth=20.0 * OPEN_WIDTH_UM, dx=dx,
                    mask_thickness=0.05, fractions=fractions,
                    label="huang_like_ar20_thinmask")
    json.dump(records, open(out / "coupled_floor_composition.json", "w"),
              indent=1)
    print("DONE")
