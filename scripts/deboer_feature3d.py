#!/usr/bin/env python3
"""de Boer SF6/O2 Si ARDE through the COMMON feature-3d engine with COUPLED ion+neutral chemistry.

Replaces the additive two-channel diagnostic (scripts/deboer_two_channel.py) with the real coupled
surface chemistry the engine already supports. The reduced coupled network of
`ReducedSiO2FluorocarbonMechanism` is generic (a coverage that neutrals BUILD and ions REMOVE, gated
by a film): reparameterized for Si-F cryo it is
  - complex_fraction  = SiF_x fluorinated-layer coverage : F radical builds it, ion removes it (high
    ion-assisted yield vs bare Si), so the etch rate is MULTIPLICATIVE in the F-built coverage x ion
    flux -- genuine synergy, not an additive sum.
  - polymer inventory = O passivation (SiO_xF_y) : O builds it, it blocks F access
    (access = exp(-passiv/monolayer)), ion sputter removes it.
The ARDE then EMERGES: F flux falls with AR (radiosity), ion flux falls with AR (angular transport),
and the ion-assisted removal is multiplicative in the F-built coverage, so the floor rate collapses
faster than either flux alone -- exactly the coupling the additive model faked.

Honest scope: parameters are declared calibration inputs (labeled nonpredictive); the material-exchange
product labels are inherited from the SiO2 mechanism (SiO2_formula_unit / fluorocarbon_film_unit) --
physically these are Si_atom / SiOxFy_passivation; velocity is label-independent and correct. A
dedicated relabeled Si class is a follow-up cleanup. Run: python scripts/deboer_feature3d.py
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np
from scipy.stats import qmc

sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from petch.boundary_state import (
    IonEnergyTransverseMaxwellianDensity, MaxwellianFluxVelocityDensity,
    PlasmaBoundaryState, SpeciesBoundaryState,
)
from petch.feature_step_3d import make_rectangular_trench_geometry_3d, solve_feature_3d
from petch.surface_kinetics import (
    EnergeticYield, ParameterEvidence, ReducedSiO2FluorocarbonMechanism,
    ReducedSiO2FluorocarbonParameters,
)


def _ev(note, supports=False):
    return ParameterEvidence(source="de Boer SF6/O2 cryo calibration closure", evidence_type="calibration_closure",
                             note=note, supports_prediction_within_declared_domain=supports)


def build_deboer_si_mechanism(*, s_F=0.06, o_dep=0.02, bare_yield=0.10, sifx_yield=1.2,
                              passiv_sputter_yield=0.3, ion_threshold_eV=15.0):
    """Coupled Si-F mechanism = the validated reduced kernel with Si physics parameters.

    F builds SiF_x coverage; ion removes SiF_x (high yield) and bare Si (low); O builds passivation
    that blocks F; ion sputters passivation. All params are declared calibration inputs.
    """
    evidence = {k: _ev(k) for k in (
        "site_density_m2", "polymer_monolayer_density_m2", "complex_formation_probability",
        "polymer_deposition_probability_on_substrate", "polymer_deposition_probability_on_polymer",
        "oxygen_polymer_etch_probability", "bare_sio2_yield", "complex_sio2_yield",
        "polymer_sputter_yield")}
    evidence["bulk_formula_density_m3"] = ParameterEvidence(
        source="Si atomic density 5.0e28 m^-3", evidence_type="material_constant_derived",
        note="crystalline Si", supports_prediction_within_declared_domain=True)
    return ReducedSiO2FluorocarbonMechanism(ReducedSiO2FluorocarbonParameters(
        site_density_m2=6.8e18,                    # Si(100) surface site density
        bulk_formula_density_m3=5.0e28,            # Si atoms / m^3
        polymer_monolayer_density_m2=6.8e18,       # passivation monolayer
        complex_formation_probability={"F": s_F},          # F builds SiF_x coverage
        polymer_deposition_probability_on_substrate={"O": o_dep},   # O builds passivation on bare Si
        polymer_deposition_probability_on_polymer={"O": o_dep},     # ... and on existing passivation
        oxygen_species="O", oxygen_polymer_etch_probability=0.0,     # O does not etch its own film
        bare_sio2_yield=EnergeticYield(bare_yield, ion_threshold_eV, 100.0, energy_exponent=0.5,
                                       angular_model="none"),
        complex_sio2_yield=EnergeticYield(sifx_yield, ion_threshold_eV, 100.0, energy_exponent=0.5,
                                          angular_model="none"),   # ion-assisted SiF_x removal (>> bare)
        polymer_sputter_yield=EnergeticYield(passiv_sputter_yield, ion_threshold_eV, 100.0,
                                             energy_exponent=0.5, angular_model="none"),
        evidence=evidence))


def thermal_neutral(name, mass_amu, flux, ref_m, *, T=0.05, log2=15, seed=0):
    density = MaxwellianFluxVelocityDensity(T)
    u = qmc.Sobol(3, scramble=True, seed=seed).random_base2(log2)
    vel = density.sample_flux_velocity(u)
    return SpeciesBoundaryState(name=name, charge_number=0, mass_amu=mass_amu, flux_m2_s=flux,
                                velocity_sqrt_eV=vel, weight=np.full(vel.shape[0], 1.0 / vel.shape[0]),
                                density_model=density, provenance={"model": "thermal_flux_neutral"})


def ion_species(flux, ref_m, *, energy_eV=100.0, iad_sigma_deg=1.0, mass_amu=131.0, log2=15, seed=1):
    # tangential temperature giving the requested cross-slot angular sigma at this energy:
    # sigma ~ sqrt(T_tang / (2 E))  ->  T_tang = 2 E sigma^2
    sigma = np.deg2rad(iad_sigma_deg)
    T_tang = max(2.0 * energy_eV * sigma * sigma, 1e-4)
    density = IonEnergyTransverseMaxwellianDensity(
        normal_energy_edges_eV=[energy_eV - 1.0, energy_eV + 1.0], probability_mass=[1.0],
        tangential_temperature_eV=T_tang)
    u = qmc.Sobol(3, scramble=True, seed=seed).random_base2(log2)
    vel = density.sample_flux_velocity(u)
    return SpeciesBoundaryState(name="ion", charge_number=1, mass_amu=mass_amu, flux_m2_s=flux,
                                velocity_sqrt_eV=vel, weight=np.full(vel.shape[0], 1.0 / vel.shape[0]),
                                density_model=density, provenance={"model": "narrow_iad_ion"})


def floor_rate(aspect_ratio, mechanism, *, opening_um=0.10, dx_um=0.01, mask_um=0.05,
               ion_flux=2e19, f_flux=2e20, o_flux=4e19, iad_sigma_deg=1.0, duration_s=1.0, seed=0):
    etched = aspect_ratio * opening_um
    substrate_top = etched + max(4.0 * dx_um, 0.05)
    geometry = make_rectangular_trench_geometry_3d(
        cell_width=2.0 * opening_um, cell_length=max(6.0 * dx_um, 0.06),
        domain_height=substrate_top + mask_um + max(6.0 * dx_um, 0.06), dx=dx_um,
        opening_width=opening_um, mask_thickness=mask_um, substrate_top=substrate_top,
        etched_depth=etched)
    domain = (np.asarray(geometry.phi.shape) - 1) * geometry.dx
    source_z = float(domain[2])
    ref_m = source_z * geometry.mesh_length_unit_m
    boundary = PlasmaBoundaryState(species=(
        ion_species(ion_flux, ref_m, iad_sigma_deg=iad_sigma_deg, seed=seed + 1),
        thermal_neutral("F", 19.0, f_flux, ref_m, seed=seed + 2),
        thermal_neutral("O", 16.0, o_flux, ref_m, seed=seed + 3)),
        reference_plane_m=ref_m)
    result = solve_feature_3d(
        geometry, boundary,
        {"ion": "energetic_bombardment", "F": "neutral_reactant", "O": "neutral_reactant"},
        mechanism, etchable_material_ids=(1,), duration_s=duration_s, n_steps=1,
        source_bounds=(0.0, float(domain[0]), 0.0, float(domain[1])), source_z=source_z,
        n_position=64, seed=seed, cfl_number=0.3, reinitialize=False, transport_device="cpu",
        neutral_radiosity_options={
            "rays_per_face": 64, "seed": seed + 100, "periodic_lateral": True,
            "domain_size": domain,
            "nonetchable_reaction_probability_by_material": {2: {"F": 1e-3, "O": 1e-3}}},
        ballistic_transport="forward", ballistic_face_quadrature_points=1,
        reinitialization_method="cr2")
    step = result.steps[-1]
    vel = step.face_velocity_mesh_units_s[step.active_face_index] * geometry.mesh_length_unit_m
    cz = step.active_face_centroid[:, 2]
    floor = cz <= float(np.min(cz)) + geometry.dx
    return float(np.mean(vel[floor]))


def deboer_calibrate_predict(*, dx_um=0.01, iad_sigma_deg=1.0, sF_grid=(0.08, 0.14, 0.22)):
    """Calibrate F sticking on the de Boer knee (AR10,20); predict AR40 held-out.

    Verified result (dx=0.01, iad=1deg): best knee s_F=0.08 gives NR10/NR20 = 0.476/0.289 vs exp
    0.43/0.29 (RMSE 0.033); held-out AR40 predicted 0.166 vs exp 0.20 (error 0.034). The AR40
    under-prediction is the sub-degree-IAD / charging frontier (Phase 2). AR40 runs in ~10 s.
    """
    exp = {10.0: 0.43, 20.0: 0.29, 40.0: 0.20}
    print("de Boer through the COUPLED engine: calibrate knee (AR10,20), predict AR40 held-out")
    print(f"{'s_F':>5} {'NR10':>6} {'NR20':>6} {'knee_RMSE':>9}")
    best = None
    for sF in sF_grid:
        m = build_deboer_si_mechanism(s_F=sF)
        r0 = floor_rate(0.0, m, dx_um=dx_um, iad_sigma_deg=iad_sigma_deg)
        nr = {ar: floor_rate(ar, m, dx_um=dx_um, iad_sigma_deg=iad_sigma_deg) / r0 for ar in (10.0, 20.0)}
        err = float(np.sqrt(np.mean([(nr[a] - exp[a]) ** 2 for a in (10.0, 20.0)])))
        print(f"{sF:>5.2f} {nr[10.0]:>6.3f} {nr[20.0]:>6.3f} {err:>9.4f}", flush=True)
        if best is None or err < best[0]:
            best = (err, sF)
    err, sF = best
    m = build_deboer_si_mechanism(s_F=sF)
    r0 = floor_rate(0.0, m, dx_um=dx_um, iad_sigma_deg=iad_sigma_deg)
    nr40 = floor_rate(40.0, m, dx_um=dx_um, iad_sigma_deg=iad_sigma_deg) / r0
    print(f"\nbest knee s_F={sF} (RMSE {err:.4f}); HELD-OUT AR40 = {nr40:.3f} vs exp 0.20 "
          f"(error {abs(nr40 - 0.20):.3f})")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--deboer", action="store_true",
                        help="calibrate the de Boer knee and predict AR40 held-out, then exit")
    parser.add_argument("--aspect-ratios", type=float, nargs="+", default=[0.0, 2.0])
    parser.add_argument("--dx-um", type=float, default=0.01)
    parser.add_argument("--iad-sigma-deg", type=float, default=1.0)
    args = parser.parse_args()
    if args.deboer:
        deboer_calibrate_predict(dx_um=args.dx_um, iad_sigma_deg=args.iad_sigma_deg)
        return
    mech = build_deboer_si_mechanism()
    print("de Boer Si SF6/O2 ARDE through feature-3d (COUPLED ion+neutral chemistry)")
    print(f"{'AR':>5} {'floor_rate_m_s':>16} {'normalized':>12}")
    rates = {}
    for ar in args.aspect_ratios:
        r = floor_rate(ar, mech, dx_um=args.dx_um, iad_sigma_deg=args.iad_sigma_deg)
        rates[ar] = r
        print(f"{ar:>5.1f} {r:>16.6e}", flush=True)
    r0 = rates[min(rates)]
    if r0 > 0:
        print("\nnormalized (rate / rate@shallowest):")
        for ar in args.aspect_ratios:
            print(f"  AR{ar:>4.1f}  {rates[ar] / r0:.4f}")


if __name__ == "__main__":
    main()
