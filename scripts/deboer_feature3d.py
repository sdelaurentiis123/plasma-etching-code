#!/usr/bin/env python3
"""de Boer SF6/O2 Si ARDE adapters for the common feature-3d engine.

``build_deboer_si_mechanism`` preserves the historical generic-SiO2 adapter so archived runs remain
replayable.  It was useful diagnostically but omits standalone fluorine chemical removal and must not
support a silicon-chemistry claim.  Its reduced network is
  - complex_fraction  = SiF_x fluorinated-layer coverage : F radical builds it, ion removes it (high
    ion-assisted yield vs bare Si), so the etch rate is MULTIPLICATIVE in the F-built coverage x ion
    flux -- genuine synergy, not an additive sum.
  - polymer inventory = O passivation (SiO_xF_y) : O builds it, it blocks F access
    (access = exp(-passiv/monolayer)), ion sputter removes it.
``build_common_belen_si_mechanism`` is the repaired production-shaped path.  It uses the explicit
common-engine Belen silicon class, including coupled F/O coverages, direct F chemical removal,
physical sputtering, and ion-enhanced removal, with a closed target-material ledger and declared
sources/bounds.  Figure-9 data already seen by development are not validation evidence for it.
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
from petch.feature_step_3d import (
    advance_feature_step_3d, make_rectangular_trench_geometry_3d, solve_feature_3d,
)
from petch.surface_kinetics import (
    EnergeticYield, ParameterEvidence, ReducedSiO2FluorocarbonMechanism,
    ReducedSiO2FluorocarbonParameters, SteinbruchelYield,
)
from petch.silicon_sf6o2 import (
    BelenSiliconParameters, BelenSiliconSF6O2Mechanism,
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


def build_common_belen_si_mechanism(
        *, s_F=0.5, s_O=1.0, k_sigma_m2_s=3.0e21,
        beta_sigma_m2_s=4.0e17, a_sp=0.0337, a_ie=7.0, a_o=3.0):
    """Return the common-engine Belen SF6/O2 silicon law.

    This is the production-shaped repair for the generic SiO2 adapter above.  ``s_F`` is the
    neutral reaction probability controlling both Knudsen/radiosity loss and the coupled site
    balance.  The default 0.5 is the cryogenic-Si value reported by Blauw et al.; it remains marked
    nonpredictive because reaction probability depends on flux, temperature, and surface state.
    Belen/ViennaPS calibration values are converted from ``1e15 cm^-2 s^-1`` to SI.
    """
    predictive = {
        "site_density_m2": True,
        "bulk_si_atom_density_m3": True,
        "fluorine_atoms_per_removed_si": True,
        "ion_enhanced_fluorine_release_per_si": True,
    }
    source = {
        "site_density_m2": "Si(100) surface-site density; material constant",
        "bulk_si_atom_density_m3": "crystalline Si atomic density 5.0e28 m^-3",
        "fluorine_sticking_probability": (
            "Blauw et al., JVST B 18, 3453 (2000), DOI 10.1116/1.1313578; "
            "effective cryogenic atomic-F reaction probability about 0.5"),
        "oxygen_sticking_probability": (
            "Belen et al., JVST A 23, 99 (2005), DOI 10.1116/1.1830495; "
            "ViennaPS calibrated closure"),
        "spontaneous_fluorine_removal_rate_m2_s": (
            "Belen et al. 2005/ViennaPS k_sigma=300 in 1e15 cm^-2 s^-1 units"),
        "oxygen_desorption_rate_m2_s": (
            "Belen et al. 2005/ViennaPS beta_sigma=0.04 in 1e15 cm^-2 s^-1 units"),
        "physical_sputter_yield": (
            "Belen et al. 2005 calibration with Steinbruchel square-root energy law"),
        "ion_enhanced_yield": (
            "Belen et al. 2005 calibration with Steinbruchel square-root energy law"),
        "oxygen_sputter_yield": (
            "Belen et al. 2005 calibration with Steinbruchel square-root energy law"),
        "fluorine_atoms_per_removed_si": "SiF4 chemical-product stoichiometry",
        "ion_enhanced_fluorine_release_per_si": "Belen coupled-coverage equation coefficient",
    }
    evidence = {
        name: ParameterEvidence(
            source=text,
            evidence_type=("material_constant_or_stoichiometry" if predictive.get(name, False)
                           else "literature_calibrated_model_input"),
            note="de Boer transfer points are not used as evidence for this value",
            supports_prediction_within_declared_domain=predictive.get(name, False))
        for name, text in source.items()}
    bounds = {
        "site_density_m2": (5.0e18, 8.0e18),
        "bulk_si_atom_density_m3": (4.9e28, 5.1e28),
        "fluorine_sticking_probability": (0.03, 0.7),
        "oxygen_sticking_probability": (0.0, 1.0),
        "spontaneous_fluorine_removal_rate_m2_s": (1.0e20, 1.0e22),
        "oxygen_desorption_rate_m2_s": (1.0e16, 1.0e19),
        "physical_sputter_yield": {
            "prefactor_per_sqrt_eV": (0.0, 0.2),
            "threshold_energy_eV": (5.0, 50.0),
            "angular_parameter": (0.0, 15.0)},
        "ion_enhanced_yield": {
            "prefactor_per_sqrt_eV": (0.0, 15.0),
            "threshold_energy_eV": (5.0, 50.0)},
        "oxygen_sputter_yield": {
            "prefactor_per_sqrt_eV": (0.0, 8.0),
            "threshold_energy_eV": (5.0, 50.0)},
        "fluorine_atoms_per_removed_si": (4.0, 4.0),
        "ion_enhanced_fluorine_release_per_si": (2.0, 2.0),
    }
    return BelenSiliconSF6O2Mechanism(BelenSiliconParameters(
        material_name="crystalline_Si", material_inventory_name="Si_atom",
        fluorine_species="F", oxygen_species="O", projectile_species=("ion",),
        site_density_m2=6.8e18, bulk_si_atom_density_m3=5.0e28,
        fluorine_sticking_probability=s_F, oxygen_sticking_probability=s_O,
        spontaneous_fluorine_removal_rate_m2_s=k_sigma_m2_s,
        oxygen_desorption_rate_m2_s=beta_sigma_m2_s,
        physical_sputter_yield=SteinbruchelYield(
            a_sp, 20.0, angular_model="kress_1999", angular_parameter=9.3),
        ion_enhanced_yield=SteinbruchelYield(
            a_ie, 15.0, angular_model="chang_sawin_1997"),
        oxygen_sputter_yield=SteinbruchelYield(
            a_o, 10.0, angular_model="chang_sawin_1997"),
        fluorine_atoms_per_removed_si=4.0,
        ion_enhanced_fluorine_release_per_si=2.0,
        evidence=evidence, parameter_bounds=bounds))


def thermal_neutral(name, mass_amu, flux, ref_m, *, T=0.05, log2=15, seed=0):
    density = MaxwellianFluxVelocityDensity(T)
    u = qmc.Sobol(3, scramble=True, seed=seed).random_base2(log2)
    vel = density.sample_flux_velocity(u)
    return SpeciesBoundaryState(name=name, charge_number=0, mass_amu=mass_amu, flux_m2_s=flux,
                                velocity_sqrt_eV=vel, weight=np.full(vel.shape[0], 1.0 / vel.shape[0]),
                                density_model=density, provenance={"model": "thermal_flux_neutral"})


def ion_species(flux, ref_m, *, energy_eV=100.0, iad_sigma_deg=1.0, mass_amu=131.0,
                log2=15, seed=1):
    """Return a provenance-bearing narrow-IED, transverse-Maxwellian ion population.

    ``iad_sigma_deg`` is the small-angle standard deviation of either transverse angular
    component, not an ambiguously defined polar-angle FWHM.  The corresponding polar magnitude is
    approximately Rayleigh distributed; its FWHM is about ``1.60 * iad_sigma_deg``.  Keeping that
    distinction explicit matters for the de Boer low-bias RIE-lag experiment, whose source paper
    discusses a roughly five-degree polar IAD width.
    """
    if (not np.isfinite(energy_eV) or energy_eV <= 1.0
            or not np.isfinite(iad_sigma_deg) or iad_sigma_deg <= 0.0
            or not np.isfinite(mass_amu) or mass_amu <= 0.0):
        raise ValueError("invalid ion energy, angular spread, or mass")
    # Tangential temperature giving the requested component angular sigma at this energy:
    # sigma ~ sqrt(T_tang / (2 E))  ->  T_tang = 2 E sigma^2.
    sigma = np.deg2rad(iad_sigma_deg)
    T_tang = max(2.0 * energy_eV * sigma * sigma, 1e-4)
    density = IonEnergyTransverseMaxwellianDensity(
        normal_energy_edges_eV=[energy_eV - 1.0, energy_eV + 1.0], probability_mass=[1.0],
        tangential_temperature_eV=T_tang)
    u = qmc.Sobol(3, scramble=True, seed=seed).random_base2(log2)
    vel = density.sample_flux_velocity(u)
    return SpeciesBoundaryState(
        name="ion", charge_number=1, mass_amu=mass_amu, flux_m2_s=flux,
        velocity_sqrt_eV=vel, weight=np.full(vel.shape[0], 1.0 / vel.shape[0]),
        density_model=density, provenance={
            "model": "narrow_ied_transverse_maxwellian_iad",
            "representative_normal_energy_eV": float(energy_eV),
            "normal_energy_half_width_eV": 1.0,
            "iad_component_sigma_deg": float(iad_sigma_deg),
            "approximate_polar_iad_fwhm_deg": float(1.6025 * iad_sigma_deg),
            "tangential_temperature_eV": float(T_tang),
        })


def floor_rate(aspect_ratio, mechanism, *, opening_um=0.10, dx_um=0.01, mask_um=0.05,
               ion_flux=2e19, f_flux=2e20, o_flux=4e19, ion_energy_eV=100.0,
               iad_sigma_deg=1.0, duration_s=1.0,
               seed=0, floor_average="face", surface_equilibration_steps=1,
               surface_fixed_point_tolerance=None, surface_fixed_point_max_iterations=20,
               charged_surface_response=None, charged_surface_response_options=None,
               return_diagnostics=False):
    """Return the trench-floor normal velocity from one common-engine step.

    ``floor_average='face'`` preserves the historical diagnostic, which gives every marching-cubes
    triangle equal weight.  ``'area'`` is the physical surface average and is the authoritative
    choice for new direct validation.  The explicit switch prevents an evidence-breaking silent
    reinterpretation of archived runs.
    """
    if floor_average not in {"face", "area"}:
        raise ValueError("floor_average must be 'face' or 'area'")
    if (int(surface_equilibration_steps) != surface_equilibration_steps
            or surface_equilibration_steps <= 0):
        raise ValueError("surface_equilibration_steps must be a positive integer")
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
        ion_species(
            ion_flux, ref_m, energy_eV=ion_energy_eV,
            iad_sigma_deg=iad_sigma_deg, seed=seed + 1),
        thermal_neutral("F", 19.0, f_flux, ref_m, seed=seed + 2),
        thermal_neutral("O", 16.0, o_flux, ref_m, seed=seed + 3)),
        reference_plane_m=ref_m)
    roles = {"ion": "energetic_bombardment", "F": "neutral_reactant", "O": "neutral_reactant"}
    common = dict(
        etchable_material_ids=(1,), duration_s=duration_s,
        source_bounds=(0.0, float(domain[0]), 0.0, float(domain[1])), source_z=source_z,
        n_position=64, seed=seed, cfl_number=0.3, reinitialize=False, transport_device="cpu",
        neutral_radiosity_options={
            "rays_per_face": 64, "seed": seed + 100, "periodic_lateral": True,
            "domain_size": domain,
            "nonetchable_reaction_probability_by_material": {2: {"F": 1e-3, "O": 1e-3}}},
        neutral_surface_fixed_point_tolerance=surface_fixed_point_tolerance,
        neutral_surface_fixed_point_max_iterations=surface_fixed_point_max_iterations,
        charged_surface_response=charged_surface_response,
        charged_surface_response_options=charged_surface_response_options,
        ballistic_transport="forward", ballistic_face_quadrature_points=1,
        reinitialization_method="cr2")
    state = None; fingerprint = None; step = None
    # Hold geometry fixed while the fast surface inventory relaxes.  Using ``step.surface.state``
    # (the state on the original mesh) rather than ``next_surface_state`` (remapped to the advected
    # mesh) is essential.  This repeatedly executes the same common transport/chemistry operator;
    # it is not a separate benchmark-specific rate law.
    for iteration in range(int(surface_equilibration_steps)):
        # A quasi-steady mechanism can update its transport-coupled coverage at zero profile time.
        # Only the final fixed-point evaluation advances its authoritative removal inventory.
        local_duration = duration_s
        if (getattr(mechanism, "quasi_steady_surface_state", False)
                and iteration + 1 < int(surface_equilibration_steps)):
            local_duration = 0.0
        step = advance_feature_step_3d(
            geometry, boundary, roles, mechanism, surface_state=state,
            surface_state_mesh_fingerprint=fingerprint,
            **dict(common, duration_s=local_duration))
        state = step.surface.state
        fingerprint = step.surface_state_mesh_fingerprint
    vel = step.face_velocity_mesh_units_s[step.active_face_index] * geometry.mesh_length_unit_m
    cz = step.active_face_centroid[:, 2]
    area = np.asarray(step.active_face_area, dtype=float)
    material = np.asarray(step.face_material_id)[step.active_face_index]
    floor = (cz <= float(np.min(cz)) + geometry.dx) & (material == 1)
    face_mean = float(np.mean(vel[floor]))
    area_mean = float(np.average(vel[floor], weights=area[floor]))
    value = face_mean if floor_average == "face" else area_mean
    if not return_diagnostics:
        return value
    floor_level = substrate_top - etched
    diagnostics = {
        "face_mean_m_s": face_mean,
        "area_mean_m_s": area_mean,
        "relative_area_vs_face": area_mean / face_mean - 1.0,
        "active_face_count": int(len(cz)),
        "floor_face_count": int(np.sum(floor)),
        "floor_area_mesh2": float(np.sum(area[floor])),
        "floor_centroid_z_min": float(np.min(cz[floor])),
        "floor_centroid_z_max": float(np.max(cz[floor])),
        "declared_floor_z": float(floor_level),
        "floor_z_offset_dx_quantiles": np.quantile(
            (cz[floor] - floor_level) / geometry.dx, [0.0, 0.25, 0.5, 0.75, 1.0]).tolist(),
        "surface_equilibration_steps": int(surface_equilibration_steps),
        "surface_fixed_point_iterations": step.diagnostics[
            "neutral_surface_fixed_point_iterations"],
        "surface_fixed_point_residual": step.diagnostics[
            "neutral_surface_fixed_point_residual"],
        "ion_energy_eV": float(ion_energy_eV),
        "iad_component_sigma_deg": float(iad_sigma_deg),
        "approximate_polar_iad_fwhm_deg": float(1.6025 * iad_sigma_deg),
        "charged_surface_response_applied": step.diagnostics[
            "charged_surface_response_applied"],
        "charged_surface_response_field": step.diagnostics[
            "charged_surface_response_field"],
        "charged_surface_response_bounces": step.diagnostics[
            "charged_surface_response_bounces"],
        "charged_surface_response_reimpact_events": step.diagnostics[
            "charged_surface_response_reimpact_events"],
        "charged_surface_response_relative_charge_error": step.diagnostics[
            "charged_surface_response_relative_charge_error"],
        "charged_surface_response_maximum_energy_error": step.diagnostics[
            "charged_surface_response_maximum_energy_error"],
        "charged_surface_response_tail_l1_error_bound": step.diagnostics[
            "charged_surface_response_tail_l1_error_bound"],
    }
    if hasattr(state, "complex_fraction"):
        diagnostics.update(
            floor_complex_fraction_mean=float(np.average(
                np.asarray(state.complex_fraction)[floor], weights=area[floor])),
            floor_polymer_units_m2_mean=float(np.average(
                np.asarray(state.polymer_units_m2)[floor], weights=area[floor])))
    surface = step.surface
    if hasattr(surface, "fluorine_coverage"):
        for label, field in (
                ("fluorine_coverage", surface.fluorine_coverage),
                ("oxygen_coverage", surface.oxygen_coverage),
                ("available_site_fraction", surface.available_site_fraction),
                ("chemical_removal_rate_m2_s", surface.chemical_removal_rate_m2_s),
                ("physical_sputter_rate_m2_s", surface.physical_sputter_rate_m2_s),
                ("ion_enhanced_removal_rate_m2_s", surface.ion_enhanced_removal_rate_m2_s),
                ("transport_fixed_point_change", np.abs(surface.transport_fixed_point_change))):
            diagnostics[f"floor_{label}_mean"] = float(np.average(
                np.asarray(field)[floor], weights=area[floor]))
    return value, diagnostics


def deboer_calibrate_predict(*, dx_um=0.01, iad_sigma_deg=1.0, sF_grid=(0.08, 0.14, 0.22)):
    """Legacy comparison to the Blauw/Clausing fitted curve—not direct experimental pixels.

    The 1/.43/.29/.20 values were historically mislabeled as measured de Boer points. They are a
    transport-model curve evaluated with a fitted sticking probability. Direct Figure-9 validation
    lives in ``scripts/deboer_2002_direct_validation.py`` and must be used for experimental claims.
    """
    reference = {10.0: 0.43, 20.0: 0.29, 40.0: 0.20}
    print("legacy Blauw/Clausing curve cross-check (not direct experimental evidence)")
    print(f"{'s_F':>5} {'NR10':>6} {'NR20':>6} {'knee_RMSE':>9}")
    best = None
    for sF in sF_grid:
        m = build_deboer_si_mechanism(s_F=sF)
        r0 = floor_rate(0.0, m, dx_um=dx_um, iad_sigma_deg=iad_sigma_deg)
        nr = {ar: floor_rate(ar, m, dx_um=dx_um, iad_sigma_deg=iad_sigma_deg) / r0 for ar in (10.0, 20.0)}
        err = float(np.sqrt(np.mean([(nr[a] - reference[a]) ** 2 for a in (10.0, 20.0)])))
        print(f"{sF:>5.2f} {nr[10.0]:>6.3f} {nr[20.0]:>6.3f} {err:>9.4f}", flush=True)
        if best is None or err < best[0]:
            best = (err, sF)
    err, sF = best
    m = build_deboer_si_mechanism(s_F=sF)
    r0 = floor_rate(0.0, m, dx_um=dx_um, iad_sigma_deg=iad_sigma_deg)
    nr40 = floor_rate(40.0, m, dx_um=dx_um, iad_sigma_deg=iad_sigma_deg) / r0
    print(f"\nbest knee s_F={sF} (RMSE {err:.4f}); AR40 model cross-check = {nr40:.3f} vs 0.20 "
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
