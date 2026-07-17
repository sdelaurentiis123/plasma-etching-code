from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
import petch.feature_step_3d as feature_step_module

from petch.boundary_state import (
    IonEnergyTransverseMaxwellianDensity, PlasmaBoundaryState, SpeciesBoundaryState,
    maxwellian_electron_boundary_state, mixture_boundary_proposal, qmc_boundary_proposal,
)
from petch.charging_poisson_3d import NodalPoissonSystem3D
from petch.charged_surface_response_3d import GrazingSpecularIonReflection3D
from petch.feature_step_3d import (
    FeatureGeometry3D,
    SurfaceTopologyChangeError,
    _apply_subcell_cleanup_to_material_levelsets,
    _apply_subcell_gas_fill_to_material_levelsets,
    _face_material_ids,
    _periodic_physical_volume_topology_signature,
    _physical_volume_topology_signature,
    _remove_unresolved_subcell_solid_components,
    _surface_gas_normals,
    _unresolved_subcell_gas_cavity_mask,
    advance_feature_step_3d,
    conservative_remap_surface_state,
    make_rectangular_trench_geometry_3d,
    solve_feature_3d,
)
from petch.boundary_transport_3d import (
    gather_boundary_state_field_adjoint_3d, trace_boundary_state_field_3d,
)
from petch.interaction_data import load_kounis_melas_2024_tables
from petch.neutral_radiosity_3d import DiffuseNeutralNoSinkError
from petch.physical_api import COMMON_FEATURE_ENGINE, PhysicalProcess
from petch.physical_sputtering import (
    PhysicalSputterMechanism, PhysicalSputterParameters,
)
from petch.surface_kinetics import (
    EnergeticYield,
    MechanismValidity,
    ParameterEvidence,
    ReducedSiO2FluorocarbonMechanism,
    ReducedSiO2FluorocarbonParameters,
    SiO2SurfaceState,
)
from petch.surface_exchange import SurfaceMaterialExchange
from petch.tabulated_chemistry import TabulatedSiClArMechanism, TabulatedSiSurfaceState
from petch.threed import advect_3d, extract_mesh_3d, reinit_narrow


@pytest.mark.parametrize(
    ("old", "new", "expected"),
    [
        ((1, 0, False, ((1, 1),)), (1, 1, False, ((1, 1),)),
         "gas_cavity_enclosed"),
        ((1, 0, False, ((1, 1),)), (1, 0, True, ((1, 1),)),
         "domain_gas_breakthrough"),
        ((1, 0, False, ((1, 1),)), (2, 0, False, ((1, 1),)),
         "solid_component_change"),
    ],
)
def test_surface_topology_change_error_exposes_geometry_event(old, new, expected):
    error = SurfaceTopologyChangeError(
        "manufactured topology event", method="manufactured",
        old_topology=old, new_topology=new,
        old_mesh_topology=(1, 1), new_mesh_topology=(1, 1),
        changed_slice_topology={})

    assert isinstance(error, ValueError)
    assert error.event_kind == expected


def test_periodic_component_size_diagnostic_merges_wrapped_faces():
    occupied = np.zeros((5, 3, 2), dtype=bool)
    occupied[0, 1, 0] = True
    occupied[-1, 1, 0] = True
    occupied[2, 1, 1] = True

    assert feature_step_module._periodic_component_sizes(occupied) == (2, 1)


def test_default_feature_policy_refuses_resolved_gas_cavity_enclosure(monkeypatch):
    opened = _periodic_pinchoff_geometry(sealed=False)
    sealed = _periodic_pinchoff_geometry(sealed=True)

    monkeypatch.setattr(
        feature_step_module, "_advect_exposed_material_levelsets",
        lambda *args, **kwargs: {1: np.asarray(sealed.phi).copy()})

    with pytest.raises(SurfaceTopologyChangeError) as info:
        advance_feature_step_3d(
            opened, _pinchoff_boundary(1.6e-6),
            {"Ar+": "energetic_bombardment", "CF2": "neutral_reactant"},
            _mechanism(), etchable_material_ids=(1,), duration_s=0.01,
            source_bounds=(-0.01, 1.21, -0.01, 0.81), source_z=1.6,
            n_position=4, seed=11, reinitialize=False,
            profile_periodic_lateral=True, transport_device="cpu",
            ballistic_transport="face_gather",
            ballistic_face_quadrature_points=1)

    assert info.value.event_kind == "gas_cavity_enclosed"


def test_explicit_cavity_policy_continues_pinchoff_and_reopening_with_closed_flux_gate(
        monkeypatch):
    opened = _periodic_pinchoff_geometry(sealed=False)
    sealed = _periodic_pinchoff_geometry(sealed=True)
    target = iter((sealed.phi, opened.phi))

    def manufactured_transition(*args, **kwargs):
        return {1: np.asarray(next(target)).copy()}

    monkeypatch.setattr(
        feature_step_module, "_advect_exposed_material_levelsets",
        manufactured_transition)
    result = solve_feature_3d(
        opened, _pinchoff_boundary(1.6e-6),
        {"Ar+": "energetic_bombardment", "CF2": "neutral_reactant"},
        _mechanism(), etchable_material_ids=(1,), duration_s=0.02, n_steps=2,
        source_bounds=(-0.01, 1.21, -0.01, 0.81), source_z=1.6,
        n_position=4, seed=11, reinitialize=False,
        profile_periodic_lateral=True, transport_device="cpu",
        ballistic_transport="face_gather", ballistic_face_quadrature_points=1,
        topology_change_policy="continue_gas_cavity")

    assert len(result.steps) == 2
    assert result.duration_s == 0.02
    assert [step.diagnostics["topology_event"]["kind"] for step in result.steps] == [
        "gas_cavity_enclosed", "gas_cavity_opened"]
    assert all(
        step.state_remap_diagnostics["materials"][1][
            "max_relative_conservation_residual"] <= 1e-14
        for step in result.steps)

    # Step two evaluates transport on the closed geometry.  The lower central
    # faces bound the sealed gas component and must receive exactly no external
    # first-hit source before the manufactured reopening is applied.
    closed_step = result.steps[1]
    centroid = closed_step.active_face_centroid
    cavity_face = (
        (centroid[:, 0] > 0.5) & (centroid[:, 0] < 0.9)
        & (centroid[:, 2] < 1.25))
    assert np.any(cavity_face)
    ion = next(
        item for item in closed_step.transport.surface_fluxes.energetic_fluxes
        if item.name == "Ar+")
    assert np.all(ion.flux_m2_s[cavity_face] == 0.0)
    assert np.array_equal(result.geometry.phi, opened.phi)


def test_cavity_continue_policy_still_refuses_solid_component_change(monkeypatch):
    opened = _periodic_pinchoff_geometry(sealed=False)
    split = np.asarray(opened.phi).copy()
    split[:, :, 2] = -1.0

    monkeypatch.setattr(
        feature_step_module, "_advect_exposed_material_levelsets",
        lambda *args, **kwargs: {1: split.copy()})

    with pytest.raises(SurfaceTopologyChangeError) as info:
        solve_feature_3d(
            opened, _pinchoff_boundary(1.6e-6),
            {"Ar+": "energetic_bombardment", "CF2": "neutral_reactant"},
            _mechanism(), etchable_material_ids=(1,), duration_s=0.01, n_steps=1,
            source_bounds=(-0.01, 1.21, -0.01, 0.81), source_z=1.6,
            n_position=4, seed=11, reinitialize=False,
            profile_periodic_lateral=True, transport_device="cpu",
            ballistic_transport="face_gather",
            ballistic_face_quadrature_points=1,
            topology_change_policy="continue_gas_cavity")

    assert info.value.event_kind in {
        "solid_component_change", "material_component_change"}


def test_public_engine_physically_closes_continues_and_reopens_keyhole():
    geometry = _public_keyhole_geometry()
    mechanism = _ManufacturedReversibleMotion()
    state = None
    fingerprint = None
    physical_time = 0.0
    trajectory = []

    def take_step(mode, duration_s):
        nonlocal geometry, state, fingerprint, physical_time
        boundary = _public_keyhole_boundary(mode, 1.0e-6)
        role = {
            species.name: (
                "neutral_reactant" if species.charge_number == 0
                else "energetic_bombardment")
            for species in boundary.species}
        result = advance_feature_step_3d(
            geometry, boundary, role, mechanism,
            etchable_material_ids=(1,), duration_s=duration_s,
            source_bounds=(-0.005, 0.605, -0.005, 0.105), source_z=1.0,
            surface_state=state,
            surface_state_mesh_fingerprint=fingerprint,
            n_position=4, seed=29, cfl_number=0.25,
            reinitialize=True, reinitialization_method="cr2",
            profile_periodic_lateral=True, transport_device="cpu",
            ballistic_transport="face_gather",
            ballistic_face_quadrature_points=1,
            topology_change_policy="continue_gas_cavity")
        geometry = result.geometry
        state = result.next_surface_state
        fingerprint = result.next_surface_state_mesh_fingerprint
        physical_time += duration_s
        trajectory.append(result)
        return result

    closure = None
    for _ in range(20):
        result = take_step("coat", 0.25)
        event = result.diagnostics["topology_event"]
        if event is not None:
            assert event["kind"] == "gas_cavity_enclosed"
            closure = result
            break
    assert closure is not None
    closure_time = physical_time

    # This is a real positive-duration step on the sealed geometry, not a
    # geometry injection or a zero-time audit: accepted physical time advances.
    sealed_step = take_step("coat", 0.10)
    assert physical_time > closure_time
    assert sealed_step.diagnostics["topology_event"] is None
    centroid = sealed_step.active_face_centroid
    cavity_face = (
        (np.abs(centroid[:, 0] - 0.3) < 0.12)
        & (centroid[:, 2] < 0.68))
    assert np.any(cavity_face)
    probe = next(
        item for item in sealed_step.transport.surface_fluxes.energetic_fluxes
        if item.name == "probe+")
    assert np.all(probe.flux_m2_s[cavity_face] == 0.0)
    assert np.all(
        sealed_step.transport.surface_fluxes.neutral_flux_m2_s[
            "coat"][cavity_face] == 0.0)

    reopened = None
    for _ in range(12):
        result = take_step("etch", 0.25)
        event = result.diagnostics["topology_event"]
        if event is not None:
            assert event["kind"] == "gas_cavity_opened"
            reopened = result
            break
    assert reopened is not None

    access = take_step("etch", 0.0)
    centroid = access.active_face_centroid
    floor_face = (
        (np.abs(centroid[:, 0] - 0.3) < 0.10)
        & (centroid[:, 2] < 0.30))
    assert np.any(floor_face)
    etchant = next(
        item for item in access.transport.surface_fluxes.energetic_fluxes
        if item.name == "etch+")
    assert np.any(etchant.flux_m2_s[floor_face] > 0.0)
    assert all(
        max(
            material["max_relative_conservation_residual"]
            for material in step.state_remap_diagnostics["materials"].values())
        <= 1e-14
        for step in trajectory)
    assert all(
        np.all(step.surface.material_exchange.residual_units_m2("solid_unit") == 0.0)
        for step in trajectory)


def test_subcell_material_cleanup_selects_only_new_unresolved_components():
    phi = np.ones((5, 5, 5))
    previous = np.ones(phi.shape, dtype=int)
    candidate = previous.copy()
    candidate[2, 2, 2] = 2

    repair, count = (
        feature_step_module._new_unresolved_subcell_material_component_mask(
            phi, candidate, previous, (1, 2), periodic_lateral=False))

    assert count == 1
    assert np.array_equal(np.argwhere(repair), [[2, 2, 2]])

    # An already-existing fragment is a physical split and remains a refusal.
    previous[2, 2, 2] = 2
    repair, count = (
        feature_step_module._new_unresolved_subcell_material_component_mask(
            phi, candidate, previous, (1, 2), periodic_lateral=False))
    assert count == 0
    assert not np.any(repair)

    # Eight nodes support one resolved hexahedral volume cell and are never cleaned.
    previous.fill(1)
    candidate.fill(1)
    candidate[1:3, 1:3, 1:3] = 2
    repair, count = (
        feature_step_module._new_unresolved_subcell_material_component_mask(
            phi, candidate, previous, (1, 2), periodic_lateral=False))
    assert count == 0
    assert not np.any(repair)


INTERACTION_DATA = (
    Path(__file__).parents[1] / "data" / "surface_interactions" / "kounis_melas_2024")
SI_ATOM_DENSITY_M3 = 8.0 / (5.43e-10) ** 3


def _evidence():
    names = {
        "site_density_m2", "bulk_formula_density_m3", "polymer_monolayer_density_m2",
        "complex_formation_probability", "polymer_deposition_probability_on_substrate",
        "polymer_deposition_probability_on_polymer", "oxygen_polymer_etch_probability",
        "bare_sio2_yield", "complex_sio2_yield", "polymer_sputter_yield",
    }
    return {name: ParameterEvidence("manufactured moving-plane gate", "analytic") for name in names}


def _mechanism():
    yield_law = EnergeticYield(0.2, 20.0, 100.0)
    return ReducedSiO2FluorocarbonMechanism(ReducedSiO2FluorocarbonParameters(
        site_density_m2=5e18, bulk_formula_density_m3=2.2e28,
        polymer_monolayer_density_m2=4e18,
        complex_formation_probability={"CF2": 0.0},
        polymer_deposition_probability_on_substrate={},
        polymer_deposition_probability_on_polymer={}, oxygen_species="O",
        oxygen_polymer_etch_probability=0.0,
        bare_sio2_yield=yield_law, complex_sio2_yield=yield_law,
        polymer_sputter_yield=yield_law, evidence=_evidence()))


def _plane_geometry():
    dx = 0.25; shape = (4, 4, 8); top = 0.95
    z = np.arange(shape[2]) * dx
    phi = np.broadcast_to(top - z, shape).copy()
    material = np.where(phi > 0.0, 1, 0)
    return FeatureGeometry3D(phi, material, dx, 1e-6), top


def _periodic_pinchoff_geometry(*, sealed):
    """Resolved translational cavity used to certify enclosure and reopening."""
    dx = 0.2
    phi = np.ones((7, 5, 9), dtype=float)
    phi[:, :, 7:] = -1.0
    phi[3:5, :, 3:7] = -1.0
    if sealed:
        phi[3:5, :, 6] = 1.0
    material = np.where(phi > 0.0, 1, 0)
    return FeatureGeometry3D(
        phi, material, dx, 1e-6, material_levelsets={1: phi})


def _pinchoff_boundary(reference_plane_m):
    ion = SpeciesBoundaryState(
        "Ar+", 1, 40.0, 2.2e21, [[0.0, 0.0, 10.0]], [1.0])
    neutral = SpeciesBoundaryState(
        "CF2", 0, 50.0, 0.0, [[0.0, 0.0, 1.0]], [1.0])
    return PlasmaBoundaryState(
        (ion, neutral), reference_plane_m=reference_plane_m)


def _public_keyhole_geometry():
    """One resolved periodic trench for a real level-set close/reopen cycle."""
    dx = 0.05
    shape = (13, 3, 21)
    x, y, z = (np.arange(size) * dx for size in shape)
    X, _, Z = np.meshgrid(x, y, z, indexing="ij")
    center = 0.3
    half_width = 0.075
    floor = 0.15
    top = 0.75
    wall_slab = np.minimum(Z - floor, top - Z)
    wall = np.minimum(wall_slab, np.abs(X - center) - half_width)
    phi = reinit_narrow(np.maximum(floor - Z, wall), dx, 2.0)
    material = np.where(phi > 0.0, 1, 0)
    return FeatureGeometry3D(
        phi, material, dx, 1e-6, material_levelsets={1: phi})


def _public_keyhole_boundary(mode, reference_plane_m):
    if mode == "coat":
        species = (
            SpeciesBoundaryState(
                "coat", 0, 50.0, 1.0e20, [[0.0, 0.0, 1.0]], [1.0]),
            SpeciesBoundaryState(
                "probe+", 1, 40.0, 1.0e20, [[0.0, 0.0, 10.0]], [1.0]),
        )
    elif mode == "etch":
        species = (SpeciesBoundaryState(
            "etch+", 1, 40.0, 1.0e20, [[0.0, 0.0, 10.0]], [1.0]),)
    else:
        raise ValueError(mode)
    return PlasmaBoundaryState(species, reference_plane_m=reference_plane_m)


class _ManufacturedReversibleMotion:
    """Test-only conformal coat plus flux-limited directional strip.

    The mechanism does not replace a numerical operator. It supplies two declared
    nonnegative normal-velocity laws to the public engine: an ideal conformal film
    closes the keyhole; a directional energetic flux removes only externally
    visible material and therefore strips the exposed cap until access reopens.
    """

    density_m3 = 1.0e28
    coat_velocity_m_s = 2.5e-8
    etch_velocity_m_s = 5.0e-8

    @staticmethod
    def initial_state(shape=()):
        return SiO2SurfaceState.bare(shape)

    @staticmethod
    def _validity():
        return MechanismValidity(
            within_declared_scope=True, reasons=(), unsupported_neutral_species=(),
            known_model_form_omissions=(
                "manufactured conformal-coat/directional-strip certification law",),
            parameter_evidence_supports_prediction=False,
            nonpredictive_parameters=("manufactured_normal_velocity",))

    def advance(self, state, fluxes, duration_s):
        shape = state.polymer_units_m2.shape
        if "coat" in fluxes.neutral_flux_m2_s:
            growth = np.full(shape, self.coat_velocity_m_s)
            etch = np.zeros(shape)
        else:
            population = next(
                item for item in fluxes.energetic_fluxes
                if item.name == "etch+")
            incident = np.asarray(population.flux_m2_s, dtype=float)
            maximum = float(np.max(incident)) if incident.size else 0.0
            etch = (
                np.zeros(shape) if maximum == 0.0 else
                self.etch_velocity_m_s * incident / maximum)
            growth = np.zeros(shape)
        removed = etch * self.density_m3 * float(duration_s)
        deposited = growth * self.density_m3 * float(duration_s)
        exchange = SurfaceMaterialExchange(
            removed_units_m2={"solid_unit": removed},
            outgoing_units_m2={},
            unresolved_units_m2={"solid_unit": removed},
            deposited_units_m2={"solid_unit": deposited},
            known_limitations=("manufactured topology-continuation gate",))
        return SimpleNamespace(
            state=state,
            etch_velocity_m_s=etch,
            normal_growth_velocity_m_s=growth,
            material_exchange=exchange,
            product_populations=(),
            validity=self._validity())


def _boundary():
    # Y=0.2 at 100 eV and Gamma=2.2e21 m^-2 s^-1 gives V=2e-8 m/s = 0.02 um/s.
    ion = SpeciesBoundaryState(
        "Ar+", 1, 40.0, 2.2e21, [[0.0, 0.0, 10.0]], [1.0])
    neutral = SpeciesBoundaryState(
        "CF2", 0, 50.0, 0.0, [[0.0, 0.0, 1.0]], [1.0])
    return PlasmaBoundaryState((ion, neutral), reference_plane_m=1.75e-6)


def _charging_boundary():
    ion = SpeciesBoundaryState(
        "Ar+", 1, 40.0, 2.2e21, [[0.0, 0.0, 10.0]], [1.0])
    electron = SpeciesBoundaryState(
        "electron", -1, 5.4858e-4, 2.2e22,
        [[0.0, 0.0, 1.0], [0.0, 0.0, np.sqrt(20.0)]], [0.9, 0.1])
    return PlasmaBoundaryState((ion, electron), reference_plane_m=1.75e-6)


def _continuous_trench_charging_boundary(reference_plane_m):
    ion_flux = 2.2e21
    ion = SpeciesBoundaryState(
        "Ar+", 1, 40.0, ion_flux, [[0.0, 0.0, 10.0]], [1.0],
        density_model=IonEnergyTransverseMaxwellianDensity(
            np.array([99.0, 101.0]), np.array([1.0]), 0.05))
    electron = maxwellian_electron_boundary_state(
        4.0, 10.0 * ion_flux, n_transverse=5, n_normal=8,
        reference_plane_m=reference_plane_m).species[0]
    return PlasmaBoundaryState((ion, electron), reference_plane_m=reference_plane_m)


def _trench_adjoint_proposals(boundary):
    physical_ion = boundary.get("Ar+")
    broad_ion = SpeciesBoundaryState(
        "Ar+", 1, physical_ion.mass_amu, physical_ion.flux_m2_s,
        [[0.0, 0.0, 10.0]], [1.0],
        density_model=IonEnergyTransverseMaxwellianDensity(
            np.array([1.0, 201.0]), np.array([1.0]), 2.0))
    ion_mixture = mixture_boundary_proposal(
        (physical_ion, broad_ion), (0.8, 0.2), name="Ar+")
    ion_proposal = qmc_boundary_proposal(ion_mixture, 10, seed=79)
    broad_electron = maxwellian_electron_boundary_state(
        20.0, boundary.get("electron").flux_m2_s,
        n_transverse=5, n_normal=8, electron_name="electron",
        reference_plane_m=boundary.reference_plane_m).species[0]
    electron_mixture = mixture_boundary_proposal(
        (boundary.get("electron"), broad_electron), (0.8, 0.2), name="electron")
    electron_proposal = qmc_boundary_proposal(electron_mixture, 10, seed=83)
    return {"Ar+": ion_proposal, "electron": electron_proposal}


def _si_cl_ar_boundary():
    ion = SpeciesBoundaryState(
        "Ar+", 1, 40.0, 2e21, [[0.0, 0.0, 10.0]], [1.0])
    chlorine = SpeciesBoundaryState(
        "Cl2", 0, 70.906, 2e22, [[0.0, 0.0, 1.0]], [1.0])
    return PlasmaBoundaryState((ion, chlorine), reference_plane_m=1.75e-6)


def _si_cl_ar_mechanism():
    table = load_kounis_melas_2024_tables(INTERACTION_DATA).reactive_ion_etch
    return TabulatedSiClArMechanism(
        table, SI_ATOM_DENSITY_M3,
        ParameterEvidence(
            "Kounis-Melas OSTI 2589032 RIE in.lammps: diamond-Si lattice a=5.43 angstrom",
            "source_derived", supports_prediction_within_declared_domain=True))


def _plane_poisson_system(geometry):
    fixed = np.zeros(geometry.phi.shape, dtype=bool); fixed[:, :, -1] = True
    # Q1 value at the cell centre; this manufactured planar gate has one dielectric material.
    phi_center = sum(
        geometry.phi[i:i + geometry.phi.shape[0] - 1,
                     j:j + geometry.phi.shape[1] - 1,
                     k:k + geometry.phi.shape[2] - 1]
        for i in (0, 1) for j in (0, 1) for k in (0, 1)) / 8.0
    epsilon_r = np.where(phi_center > 0.0, 3.9, 1.0)
    return NodalPoissonSystem3D(
        epsilon_r,
        geometry.dx * geometry.mesh_length_unit_m, fixed)


def _area_weighted_height(phi, dx):
    _, _, centroids, areas = extract_mesh_3d(phi, dx)
    return float(np.dot(centroids[:, 2], areas) / areas.sum())


@pytest.mark.parametrize("dx", (0.02, 0.01, 0.005))
@pytest.mark.parametrize("opening_width", (0.06, 0.08, 0.10, 0.15, 0.18, 0.20))
def test_rectangular_trench_is_one_connected_substrate_at_jeon_widths(dx, opening_width):
    geometry = make_rectangular_trench_geometry_3d(
        cell_width=0.5, cell_length=3.0 * dx, domain_height=2.35, dx=dx,
        opening_width=opening_width, mask_thickness=0.7,
        substrate_top=1.4, etched_depth=3.0 * dx)

    assert _physical_volume_topology_signature(geometry, (1,)) == (1, 1)


def test_rectangular_trench_nodal_grid_includes_requested_physical_endpoints():
    geometry = make_rectangular_trench_geometry_3d(
        cell_width=0.5, cell_length=0.1, domain_height=2.35, dx=0.01,
        opening_width=0.08, mask_thickness=0.7,
        substrate_top=1.4, etched_depth=0.03)

    assert np.array_equal(
        (np.asarray(geometry.phi.shape) - 1) * geometry.dx,
        np.array([0.5, 0.1, 2.35]))


def test_periodic_topology_does_not_invent_event_as_trench_becomes_resolved():
    common = dict(
        cell_width=0.13, cell_length=0.02, domain_height=2.0, dx=0.01,
        opening_width=0.09, mask_thickness=0.85, substrate_top=1.0)
    unetched = make_rectangular_trench_geometry_3d(
        **common, etched_depth=0.0)
    resolved = make_rectangular_trench_geometry_3d(
        **common, etched_depth=0.06)

    assert _periodic_physical_volume_topology_signature(
        unetched, (1, 2)) == _periodic_physical_volume_topology_signature(
            resolved, (1, 2))


def test_unetched_opening_surface_is_owned_by_substrate_not_adjacent_mask():
    geometry = make_rectangular_trench_geometry_3d(
        cell_width=0.5, cell_length=0.1, domain_height=2.35, dx=0.02,
        opening_width=0.2, mask_thickness=0.7, substrate_top=1.4, etched_depth=0.0)
    _, faces, centroids, _ = extract_mesh_3d(geometry.phi, geometry.dx)
    material = _face_material_ids(centroids, geometry)
    opening = ((np.abs(centroids[:, 0] - 0.25) < 0.09)
               & (np.abs(centroids[:, 2] - 1.4) < 0.03))

    assert np.any(opening)
    assert np.all(material[opening] == 1)


def test_one_physical_3d_step_moves_a_uniform_sio2_plane_by_flux_yield_over_density():
    geometry, initial_height = _plane_geometry()
    result = advance_feature_step_3d(
        geometry, _boundary(),
        {"Ar+": "energetic_bombardment", "CF2": "neutral_reactant"},
        _mechanism(), etchable_material_ids=(1,), duration_s=1.0,
        source_bounds=(0.0, 0.75, 0.0, 0.75), source_z=1.75,
        n_position=16384, seed=3, cfl_number=0.3, reinitialize=False,
        transport_device="cpu")

    final_height = _area_weighted_height(result.geometry.phi, geometry.dx)
    assert np.isclose(initial_height - final_height, 0.02, atol=0.002)
    assert np.isclose(result.diagnostics["max_velocity_m_s"], 2e-8, rtol=0.08)
    assert result.diagnostics["cfl_substeps"] == 1
    assert result.validity.within_declared_scope
    assert not result.validity.parameter_evidence_supports_prediction
    assert "bare_sio2_yield" in result.validity.nonpredictive_parameters
    assert "conservative surface-state remap" in " ".join(result.validity.known_limitations)
    assert "product identities and branching are unresolved" in " ".join(
        result.validity.known_limitations)
    assert result.diagnostics["product_routing_complete"] is False
    assert result.state_remap_diagnostics["old_topology"] == (1, 1)
    assert result.state_remap_diagnostics["new_topology"] == (1, 1)


def test_adaptive_outer_step_replays_large_coupling_move_without_changing_operator():
    geometry, initial_height = _plane_geometry()
    result = solve_feature_3d(
        geometry, _boundary(),
        {"Ar+": "energetic_bombardment", "CF2": "neutral_reactant"},
        _mechanism(), etchable_material_ids=(1,), duration_s=20.0, n_steps=1,
        source_bounds=(0.0, 0.75, 0.0, 0.75), source_z=1.75,
        n_position=16, seed=3, cfl_number=0.3, reinitialize=False,
        transport_device="cpu", ballistic_transport="face_gather",
        adaptive_timestep_options={
            "initial_step_duration_s": 20.0,
            "minimum_step_duration_s": 0.1,
            "maximum_step_duration_s": 20.0,
            "target_displacement_cells": 0.35,
            "maximum_displacement_cells": 0.75,
        })

    accepted_dt = np.asarray([
        step.diagnostics["accepted_step_duration_s"] for step in result.steps])
    assert len(result.steps) > 1
    assert np.isclose(accepted_dt.sum(), 20.0, rtol=0.0, atol=1e-12)
    assert result.steps[0].diagnostics["adaptive_retry_count"] >= 1
    assert result.steps[0].diagnostics["adaptive_rejected_trials"][0][
        "classification"] == "inline_recovery_retry"
    assert all(
        step.diagnostics["max_displacement_mesh_units"]
        <= 0.75 * geometry.dx * (1.0 + 1e-12)
        for step in result.steps)
    final_height = _area_weighted_height(result.geometry.phi, geometry.dx)
    assert np.isclose(initial_height - final_height, 0.4, atol=0.01)


def test_adaptive_outer_step_attempts_declared_minimum_before_refusing(monkeypatch):
    geometry, _ = _plane_geometry()
    original = feature_step_module.advance_feature_step_3d
    attempted_duration = []

    def manufactured_topology_threshold(*args, **kwargs):
        duration = float(kwargs["duration_s"])
        attempted_duration.append(duration)
        if duration > 0.1:
            raise ValueError(
                "surface topology changed under manufactured minimum-step gate")
        return original(*args, **kwargs)

    monkeypatch.setattr(
        feature_step_module, "advance_feature_step_3d",
        manufactured_topology_threshold)
    result = solve_feature_3d(
        geometry, _boundary(),
        {"Ar+": "energetic_bombardment", "CF2": "neutral_reactant"},
        _mechanism(), etchable_material_ids=(1,), duration_s=0.2, n_steps=1,
        source_bounds=(0.0, 0.75, 0.0, 0.75), source_z=1.75,
        n_position=4, seed=3, cfl_number=0.3, reinitialize=False,
        transport_device="cpu", ballistic_transport="face_gather",
        adaptive_timestep_options={
            "initial_step_duration_s": 0.2,
            "minimum_step_duration_s": 0.1,
            "maximum_step_duration_s": 0.2,
            "target_displacement_cells": 0.35,
            "maximum_displacement_cells": 0.75,
            "shrink_factor": 0.25,
        })

    assert attempted_duration[:2] == [0.2, 0.1]
    assert np.isclose(sum(
        step.diagnostics["accepted_step_duration_s"]
        for step in result.steps), 0.2, rtol=0.0, atol=1e-15)
    assert result.steps[0].diagnostics["adaptive_retry_count"] == 1


def test_cfl_and_outer_displacement_use_grid_resolved_velocity_not_raw_face_max(
        monkeypatch):
    geometry, _ = _plane_geometry()

    def manufactured_grid_projection(values, centroids, extension_geometry, band):
        return np.full_like(extension_geometry["phi"], 0.001)

    monkeypatch.setattr(
        feature_step_module, "extend_velocity_3d", manufactured_grid_projection)
    result = advance_feature_step_3d(
        geometry, _boundary(),
        {"Ar+": "energetic_bombardment", "CF2": "neutral_reactant"},
        _mechanism(), etchable_material_ids=(1,), duration_s=1.0,
        source_bounds=(0.0, 0.75, 0.0, 0.75), source_z=1.75,
        n_position=4, seed=3, cfl_number=0.3, reinitialize=False,
        transport_device="cpu", ballistic_transport="face_gather")

    assert np.isclose(result.diagnostics["raw_maximum_face_velocity_m_s"], 2e-8)
    assert np.isclose(result.diagnostics["max_velocity_m_s"], 1e-9)
    assert np.isclose(result.diagnostics["max_displacement_mesh_units"], 0.001)
    assert result.diagnostics["cfl_substeps"] == 1


def test_periodic_subcell_gas_cleanup_selects_only_unresolved_cavities():
    phi = np.ones((5, 5, 6))
    phi[:, :, -1] = -1.0
    # One physical node stored on both copies of the periodic y endpoint.
    phi[2, 0, 2] = -1e-4
    phi[2, -1, 2] = -1e-4
    mask, count = _unresolved_subcell_gas_cavity_mask(
        phi, periodic_lateral=True)

    assert count == 1
    assert mask[2, 0, 2]
    assert mask[2, -1, 2]
    assert np.count_nonzero(mask) == 2

    layers, healed, owner = _apply_subcell_gas_fill_to_material_levelsets(
        {1: phi}, mask, (1,), 1.0, "cr2", True)
    healed_mask, healed_count = _unresolved_subcell_gas_cavity_mask(
        healed, periodic_lateral=True)
    assert healed_count == 0
    assert not np.any(healed_mask)
    assert np.all(owner[mask] == 1)
    assert np.array_equal(layers[1][0, :, :], layers[1][-1, :, :])
    assert np.array_equal(layers[1][:, 0, :], layers[1][:, -1, :])

    resolved = np.ones((6, 6, 7))
    resolved[:, :, -1] = -1.0
    resolved[2:4, 2:4, 2:4] = -1.0
    resolved_mask, resolved_count = _unresolved_subcell_gas_cavity_mask(
        resolved, periodic_lateral=True)
    assert resolved_count == 0
    assert not np.any(resolved_mask)


def test_ordinary_feature_step_reuses_certified_charged_response_lineage():
    geometry, _ = _plane_geometry()
    response = GrazingSpecularIonReflection3D.literature_bounded_sensitivity(
        1, ion_species_name="Ar+")
    result = advance_feature_step_3d(
        geometry, _boundary(),
        {"Ar+": "energetic_bombardment", "CF2": "neutral_reactant"},
        _mechanism(), etchable_material_ids=(1,), duration_s=0.0,
        source_bounds=(0.0, 0.75, 0.0, 0.75), source_z=1.75,
        n_position=64, seed=3, reinitialize=False, transport_device="cpu",
        charged_surface_response=response,
        charged_surface_response_options={
            "fixed_dt": 0.01, "max_steps": 256, "max_bounces": 16})

    cascade = result.charged_surface_cascade
    primary_ion = next(
        population for population in result.transport.surface_fluxes.energetic_fluxes
        if population.name == "Ar+")
    assert cascade is not None and cascade.completed
    assert np.array_equal(result.geometry.phi, geometry.phi)
    assert np.array_equal(result.geometry.material_id, geometry.material_id)
    assert primary_ion.event_position is not None
    assert result.diagnostics["charged_surface_response_applied"] is True
    assert result.diagnostics["charged_surface_response_field"] == "explicit_zero_field"
    assert result.diagnostics["charged_surface_response_bounces"] == 1
    assert cascade.relative_charge_balance_error < 5e-15
    assert max(
        transfer.relative_kinetic_energy_balance_error
        for transfer in cascade.transfers) < 5e-15
    assert "charged_surface_reimpact_cascade" in result.transport.transport_model
    assert "no surface reflection or neutral re-emission" not in result.transport.known_limitations


def test_periodic_godunov_advection_wraps_unique_nodal_core_and_closes_seam():
    shape = (7, 6, 4)
    x = 2.0 * np.pi * np.arange(shape[0] - 1) / (shape[0] - 1)
    y = 2.0 * np.pi * np.arange(shape[1] - 1) / (shape[1] - 1)
    core = (np.sin(x)[:, None, None] + 0.3 * np.cos(y)[None, :, None]
            + 0.1 * np.arange(shape[2])[None, None, :])
    phi = np.empty(shape)
    phi[:-1, :-1] = core
    phi[-1, :-1] = core[0]
    phi[:-1, -1] = core[:, 0]
    phi[-1, -1] = core[0, 0]
    speed = np.ones(shape)

    periodic = advect_3d(phi, speed, 1.0, 0.05, periodic_axes=(0, 1))
    open_boundary = advect_3d(phi, speed, 1.0, 0.05)

    assert np.array_equal(periodic[0], periodic[-1])
    assert np.array_equal(periodic[:, 0], periodic[:, -1])
    assert not np.array_equal(periodic[0, 2], open_boundary[0, 2])


def test_multistep_periodic_profile_keeps_duplicate_endpoint_planes_identical():
    geometry, _ = _plane_geometry()
    response = GrazingSpecularIonReflection3D.literature_bounded_sensitivity(
        1, ion_species_name="Ar+")
    result = solve_feature_3d(
        geometry, _boundary(),
        {"Ar+": "energetic_bombardment", "CF2": "neutral_reactant"},
        _mechanism(), etchable_material_ids=(1,), duration_s=1.0, n_steps=2,
        source_bounds=(0.0, 0.75, 0.0, 0.75), source_z=1.75,
        n_position=64, seed=3, reinitialize=True, reinitialization_method="cr2",
        transport_device="cpu", charged_surface_response=response,
        charged_surface_response_options={
            "fixed_dt": 0.01, "max_steps": 256, "max_bounces": 16,
            "periodic_lateral": True})

    assert np.array_equal(result.geometry.phi[0], result.geometry.phi[-1])
    assert np.array_equal(result.geometry.phi[:, 0], result.geometry.phi[:, -1])
    for step in result.steps:
        assert step.diagnostics["profile_periodic_lateral"] is True
        assert (step.diagnostics["periodic_seam_projection_max_mesh_units"]
                < 0.25 * geometry.dx)
        assert step.charged_surface_cascade.completed


def test_zero_duration_periodic_audit_does_not_compare_velocity_to_cell_length(
        monkeypatch):
    geometry, _ = _plane_geometry()

    def asymmetric_extension(_velocity, _centroid, extension_geometry, _distance):
        output = np.zeros_like(extension_geometry["phi"])
        output[-1] = 2.0
        return output

    monkeypatch.setattr(
        "petch.feature_step_3d.extend_velocity_3d", asymmetric_extension)
    result = advance_feature_step_3d(
        geometry, _boundary(),
        {"Ar+": "energetic_bombardment", "CF2": "neutral_reactant"},
        _mechanism(), etchable_material_ids=(1,), duration_s=0.0,
        source_bounds=(0.0, 0.75, 0.0, 0.75), source_z=1.75,
        n_position=16, seed=3, reinitialize=False,
        profile_periodic_lateral=True, transport_device="cpu")

    assert result.diagnostics["periodic_seam_projection_max_mesh_units"] == 0.0
    assert (
        result.diagnostics[
            "periodic_seam_velocity_projection_max_mesh_units_s"] == 1.0)
    assert (
        result.diagnostics[
            "periodic_seam_velocity_projection_displacement_mesh_units"] == 0.0)
    assert np.array_equal(result.geometry.phi, geometry.phi)


def test_public_physical_process_uses_common_engine_and_exposes_validity():
    geometry, _ = _plane_geometry()
    process = PhysicalProcess(
        geometry=geometry, boundary=_boundary(),
        species_role={"Ar+": "energetic_bombardment", "CF2": "neutral_reactant"},
        mechanism=_mechanism(), etchable_material_ids=(1,), duration_s=1.0, n_steps=1,
        source_bounds=(0.0, 0.75, 0.0, 0.75), source_z=1.75,
        solver_options={
            "n_position": 16384, "seed": 3, "cfl_number": 0.3,
            "reinitialize": False, "transport_device": "cpu",
        })

    result = process.run()

    assert process.engine == COMMON_FEATURE_ENGINE
    assert result.engine == COMMON_FEATURE_ENGINE
    assert result.run_manifest["engine"] == COMMON_FEATURE_ENGINE
    assert result.run_manifest["initial_geometry"]["phi"]["sha256"]
    assert result.run_manifest["surface_mechanism"]["type"]
    assert result.validity.within_declared_scope
    assert not result.validity.parameter_evidence_supports_prediction
    assert len(result.steps) == 1
    assert result.wall_time_s >= 0.0


def test_physical_sputter_mechanism_uses_same_feature_engine_and_reports_product_readiness():
    evidence = {
        name: ParameterEvidence(
            "manufactured common-engine sputter gate", "analytic",
            supports_prediction_within_declared_domain=True)
        for name in (
            "bulk_material_unit_density_m3", "sputter_yield",
            "emitted_product_mass_amu", "emission_angular_model", "emission_energy_model")}
    mechanism = PhysicalSputterMechanism(PhysicalSputterParameters(
        material_name="SiO2", material_inventory_name="SiO2_formula_unit",
        projectile_species=("Ar+",), bulk_material_unit_density_m3=2.2e28,
        sputter_yield=EnergeticYield(0.2, 20.0, 100.0),
        emitted_product_name="sputtered_SiO2_unit", emitted_product_mass_amu=60.084,
        emitted_material_units_per_particle=1.0, emission_angular_model="diffuse_cosine",
        emission_energy_model="thompson", emission_energy_parameters={
            "surface_binding_energy_eV": 4.7, "maximum_energy_eV": 100.0},
        evidence=evidence))
    geometry, _ = _plane_geometry()
    ion_boundary = PlasmaBoundaryState(
        (_boundary().species[0],), reference_plane_m=1.75e-6)
    result = advance_feature_step_3d(
        geometry, ion_boundary, {"Ar+": "energetic_bombardment"},
        mechanism, etchable_material_ids=(1,), duration_s=1.0,
        source_bounds=(0.0, 0.75, 0.0, 0.75), source_z=1.75,
        n_position=16384, seed=3, cfl_number=0.3, reinitialize=False,
        transport_device="cpu")

    assert result.surface.material_exchange.product_routing_complete
    assert result.diagnostics["product_population_count"] == 1
    assert result.diagnostics["product_transport_ready"] is True


def test_deterministic_face_gather_feature_step_is_independent_of_forward_particle_budget():
    geometry, _ = _plane_geometry()
    common = dict(
        geometry=geometry, boundary=_boundary(),
        species_role={"Ar+": "energetic_bombardment", "CF2": "neutral_reactant"},
        mechanism=_mechanism(), etchable_material_ids=(1,), duration_s=1.0,
        source_bounds=(0.0, 0.75, 0.0, 0.75), source_z=1.75,
        ballistic_transport="face_gather", ballistic_face_quadrature_points=3,
        cfl_number=0.3, reinitialize=False, transport_device="cpu")

    one = advance_feature_step_3d(**common, n_position=1, seed=3)
    many = advance_feature_step_3d(**common, n_position=1024, seed=99)

    assert one.transport.transport_model == "collisionless_deterministic_face_gather_3d"
    assert np.array_equal(one.face_velocity_mesh_units_s, many.face_velocity_mesh_units_s)
    assert np.array_equal(one.geometry.phi, many.geometry.phi)


def test_deterministic_face_gather_moves_plane_by_same_dimensional_flux_law():
    geometry, initial_height = _plane_geometry()
    result = advance_feature_step_3d(
        geometry, _boundary(),
        {"Ar+": "energetic_bombardment", "CF2": "neutral_reactant"},
        _mechanism(), etchable_material_ids=(1,), duration_s=1.0,
        source_bounds=(0.0, 0.75, 0.0, 0.75), source_z=1.75,
        ballistic_transport="face_gather", ballistic_face_quadrature_points=3,
        seed=3, cfl_number=0.3, reinitialize=False, transport_device="cpu")

    final_height = _area_weighted_height(result.geometry.phi, geometry.dx)
    assert np.isclose(initial_height - final_height, 0.02, atol=0.002)
    assert np.isclose(result.diagnostics["max_velocity_m_s"], 2e-8, rtol=1e-12)
    assert result.transport.transport_model == "collisionless_deterministic_face_gather_3d"


def test_feature_step_diffusely_reemits_unreacted_neutrals_with_global_balance():
    geometry, _ = _plane_geometry()
    ion = SpeciesBoundaryState(
        "Ar+", 1, 40.0, 2.2e21, [[0.0, 0.0, 10.0]], [1.0])
    neutral = SpeciesBoundaryState(
        "CF2", 0, 50.0, 3e20, [[0.0, 0.0, 1.0]], [1.0])
    boundary = PlasmaBoundaryState((ion, neutral), reference_plane_m=1.75e-6)
    result = advance_feature_step_3d(
        geometry, boundary,
        {"Ar+": "energetic_bombardment", "CF2": "neutral_reactant"},
        _mechanism(), etchable_material_ids=(1,), duration_s=0.0,
        source_bounds=(0.0, 0.75, 0.0, 0.75), source_z=1.75,
        n_position=64, seed=3, reinitialize=False, transport_device="cpu",
        neutral_radiosity_options={"rays_per_face": 16, "seed": 5})

    audit = result.diagnostics["neutral_radiosity"]["CF2"]
    assert audit["source_rate_s"] > 0.0
    assert audit["reacted_rate_s"] == 0.0
    assert np.isclose(audit["source_rate_s"], audit["escaped_rate_s"], rtol=1e-12)
    assert audit["relative_balance_error"] < 1e-12
    assert audit["solver_method"] == "analytic_zero_reaction_elision"
    assert audit["repeated_incident_flux_elided"] is True
    direct = result.transport.surface_fluxes.neutral_flux_m2_s["CF2"]
    assert np.all(np.isfinite(direct))
    assert "flux_conservative_diffuse_radiosity" in result.transport.transport_model


def test_periodic_ballistic_first_hit_is_independent_of_radiosity(monkeypatch):
    geometry, _ = _plane_geometry()
    original = feature_step_module.trace_boundary_state_first_hit_3d
    calls = []

    def capture_periodicity(*args, **kwargs):
        calls.append((kwargs.get("periodic_lateral"), kwargs.get("domain_size")))
        return original(*args, **kwargs)

    monkeypatch.setattr(
        feature_step_module, "trace_boundary_state_first_hit_3d", capture_periodicity)
    result = advance_feature_step_3d(
        geometry, _boundary(),
        {"Ar+": "energetic_bombardment", "CF2": "neutral_reactant"},
        _mechanism(), etchable_material_ids=(1,), duration_s=0.0,
        source_bounds=(0.0, 0.75, 0.0, 0.75), source_z=1.75,
        n_position=16, seed=3, reinitialize=False, transport_device="cpu",
        ballistic_periodic_lateral=True, profile_periodic_lateral=True)

    assert len(calls) == 1
    assert calls[0][0] is True
    assert np.array_equal(calls[0][1], np.array([0.75, 0.75, 1.75]))
    assert result.diagnostics["ballistic_periodic_lateral"] is True
    assert result.diagnostics["neutral_radiosity"] == {}
    assert "periodic_cell" in result.transport.transport_model


def test_periodic_radiosity_cannot_disable_periodic_ballistic_first_hit():
    geometry, _ = _plane_geometry()
    with pytest.raises(
            ValueError,
            match="periodic neutral radiosity requires periodic ballistic first-hit"):
        advance_feature_step_3d(
            geometry, _boundary(),
            {"Ar+": "energetic_bombardment", "CF2": "neutral_reactant"},
            _mechanism(), etchable_material_ids=(1,), duration_s=0.0,
            source_bounds=(0.0, 0.75, 0.0, 0.75), source_z=1.75,
            n_position=4, seed=3, reinitialize=False, transport_device="cpu",
            ballistic_periodic_lateral=False,
            neutral_radiosity_options={
                "rays_per_face": 8, "seed": 5, "periodic_lateral": True})


def test_ballistic_periodic_control_refuses_field_path_where_it_would_be_ignored():
    geometry, _ = _plane_geometry()
    with pytest.raises(ValueError, match="controls only field-free first-hit"):
        advance_feature_step_3d(
            geometry, _boundary(),
            {"Ar+": "energetic_bombardment", "CF2": "neutral_reactant"},
            _mechanism(), etchable_material_ids=(1,), duration_s=0.0,
            source_bounds=(0.0, 0.75, 0.0, 0.75), source_z=1.75,
            n_position=4, seed=3, reinitialize=False, transport_device="cpu",
            ballistic_periodic_lateral=True,
            nodal_potential_v=np.zeros(geometry.phi.shape),
            potential_origin=(0.0, 0.0, 0.0), potential_spacing=geometry.dx,
            trajectory_fixed_dt=0.01)


def test_feature_step_refines_sampled_no_sink_radiosity_class(monkeypatch):
    geometry, _ = _plane_geometry()
    original = feature_step_module.solve_diffuse_neutral_radiosity_3d
    call_count = [0]

    def one_sampled_graph_failure(*args, **kwargs):
        call_count[0] += 1
        if call_count[0] == 1:
            raise DiffuseNeutralNoSinkError(2)
        return original(*args, **kwargs)

    monkeypatch.setattr(
        feature_step_module, "solve_diffuse_neutral_radiosity_3d",
        one_sampled_graph_failure)
    reactive = ReducedSiO2FluorocarbonMechanism(replace(
        _mechanism().parameters,
        complex_formation_probability={"CF2": 0.2}))
    result = advance_feature_step_3d(
        geometry, _boundary(),
        {"Ar+": "energetic_bombardment", "CF2": "neutral_reactant"},
        reactive, etchable_material_ids=(1,), duration_s=0.0,
        source_bounds=(0.0, 0.75, 0.0, 0.75), source_z=1.75,
        n_position=4, seed=3, reinitialize=False, transport_device="cpu",
        neutral_radiosity_options={
            "rays_per_face": 8, "maximum_rays_per_face": 16, "seed": 5})

    audit = result.diagnostics["neutral_radiosity"]["CF2"]
    assert call_count[0] == 2
    assert audit["form_factor_rays_per_face"] == 16
    assert audit["form_factor_refinement_count"] == 1
    assert audit["form_factor_refinement"][0]["classification"] == (
        "nested_form_factor_refinement")


def test_feature_step_radiosity_requires_explicit_probability_for_every_pinned_material():
    geometry, _ = _plane_geometry()
    material = np.array(geometry.material_id, copy=True)
    material[2:, :, :] = np.where(material[2:, :, :] > 0, 2, 0)
    mixed = FeatureGeometry3D(
        geometry.phi, material, geometry.dx, geometry.mesh_length_unit_m)
    with pytest.raises(ValueError, match="missing neutral reaction probability for material 2"):
        advance_feature_step_3d(
            mixed, _boundary(),
            {"Ar+": "energetic_bombardment", "CF2": "neutral_reactant"},
            _mechanism(), etchable_material_ids=(1,), duration_s=0.0,
            source_bounds=(0.0, 0.75, 0.0, 0.75), source_z=1.75,
            n_position=16, seed=3, reinitialize=False, transport_device="cpu",
            neutral_radiosity_options={"rays_per_face": 8, "seed": 5})


def _static_trench_floor_neutral_flux(opening_width, *, rays_per_face=32):
    geometry = make_rectangular_trench_geometry_3d(
        cell_width=0.8, cell_length=0.15, domain_height=1.4, dx=0.05,
        opening_width=opening_width, mask_thickness=0.2,
        substrate_top=0.9, etched_depth=0.5)
    source_z = geometry.phi.shape[2] * geometry.dx - geometry.dx
    quadrature = maxwellian_electron_boundary_state(
        0.026, 3e20, n_transverse=3, n_normal=4,
        reference_plane_m=source_z * geometry.mesh_length_unit_m).species[0]
    neutral = SpeciesBoundaryState(
        "CF2", 0, 50.0, quadrature.flux_m2_s,
        quadrature.velocity_sqrt_eV, quadrature.weight)
    ion = SpeciesBoundaryState(
        "Ar+", 1, 40.0, 1e19, [[0.0, 0.0, 30.0]], [1.0])
    boundary = PlasmaBoundaryState(
        (ion, neutral), reference_plane_m=source_z * geometry.mesh_length_unit_m)
    parameters = replace(
        _mechanism().parameters, complex_formation_probability={"CF2": 0.2})
    mechanism = ReducedSiO2FluorocarbonMechanism(parameters)
    result = advance_feature_step_3d(
        geometry, boundary,
        {"Ar+": "energetic_bombardment", "CF2": "neutral_reactant"},
        mechanism, etchable_material_ids=(1,), duration_s=0.0,
        source_bounds=(0.0, 0.8, 0.0, 0.15), source_z=source_z,
        n_position=32, seed=7, reinitialize=False, transport_device="cpu",
        neutral_radiosity_options={
            "rays_per_face": rays_per_face, "seed": 11, "periodic_lateral": True,
            "domain_size": (np.asarray(geometry.phi.shape) - 1) * geometry.dx,
            "nonetchable_reaction_probability_by_material": {2: {"CF2": 0.2}},
        })
    floor = result.active_face_centroid[:, 2] < 0.45
    flux = result.transport.surface_fluxes.neutral_flux_m2_s["CF2"][
        result.active_face_index[floor]]
    return float(np.average(flux, weights=result.active_face_area[floor]))


def test_diffuse_neutral_transport_widens_from_local_plane_to_trench_width_ring():
    narrow = _static_trench_floor_neutral_flux(0.2)
    wide = _static_trench_floor_neutral_flux(0.4)

    assert 0.0 < narrow < wide < 3e20


def test_unetched_rectangular_trench_mesh_drops_only_zero_measure_csg_faces():
    geometry = make_rectangular_trench_geometry_3d(
        cell_width=0.5, cell_length=0.06, domain_height=2.35, dx=0.02,
        opening_width=0.2, mask_thickness=0.7, substrate_top=1.4, etched_depth=0.0)
    _, _, _, areas = extract_mesh_3d(geometry.phi, geometry.dx)

    assert areas.size > 0
    assert np.all(areas > 0.0)


def test_subcell_solid_component_is_removed_but_one_cell_support_is_preserved():
    phi = -np.ones((5, 5, 5))
    material = np.ones_like(phi, dtype=int)
    phi[0:2, 0:2, 0:2] = 1.0
    phi[4, 4, 4] = 0.1

    cleaned, removed, removal_mask = _remove_unresolved_subcell_solid_components(
        phi, material, (1,), 1.0)

    assert removed == 1
    assert np.count_nonzero(removal_mask) == 1
    assert removal_mask[4, 4, 4]
    assert not np.any(removal_mask[0:2, 0:2, 0:2])
    assert np.all(cleaned[0:2, 0:2, 0:2] > 0.0)
    assert cleaned[4, 4, 4] < 0.0


def test_subcell_solid_cleanup_uses_physical_periodic_connectivity():
    phi = -np.ones((5, 5, 5))
    material = np.zeros(phi.shape, dtype=int)
    # y=4 duplicates y=0. Each bounded strip is subcell-sized, while the
    # unique y=0 and y=3 strips form one resolved component through wrapping.
    for y_index in (0, 3, 4):
        phi[1:3, y_index, 1:3] = 1.0
        material[1:3, y_index, 1:3] = 1

    _, bounded_removed, _ = _remove_unresolved_subcell_solid_components(
        phi, material, (1,), 0.25)
    cleaned, periodic_removed, periodic_mask = (
        _remove_unresolved_subcell_solid_components(
            phi, material, (1,), 0.25, periodic_lateral=True))

    assert bounded_removed == 4
    assert periodic_removed == 0
    assert not np.any(periodic_mask)
    assert np.array_equal(cleaned, phi)


def test_subcell_cleanup_updates_the_authoritative_material_levelset():
    oxide = -np.ones((7, 7, 7))
    oxide[1:3, 1:3, 1:3] = 1.0
    oxide[5, 5, 5] = 0.1
    mask = -2.0 * np.ones_like(oxide)
    combined = np.maximum(oxide, mask)
    owner = np.where(combined >= 0.0, 1, 0)
    _, removed, removal_mask = _remove_unresolved_subcell_solid_components(
        combined, owner, (1,), 1.0)

    layers, cleaned, cleaned_owner = _apply_subcell_cleanup_to_material_levelsets(
        {1: oxide, 2: mask}, removal_mask, owner, (1,), 1.0, "cr2", False)

    assert removed == 1
    assert layers[1][5, 5, 5] < 0.0
    assert cleaned[5, 5, 5] < 0.0
    assert cleaned_owner[5, 5, 5] == 0
    assert np.all(cleaned_owner[1:3, 1:3, 1:3] == 1)


def test_diffuse_neutral_trench_floor_flux_converges_with_form_factor_rule():
    coarse = _static_trench_floor_neutral_flux(0.2, rays_per_face=16)
    fine = _static_trench_floor_neutral_flux(0.2, rays_per_face=32)

    assert abs(coarse - fine) / fine < 0.15


def test_feature_step_refuses_surface_history_without_matching_mesh_fingerprint():
    geometry, _ = _plane_geometry()
    from petch.surface_kinetics import SiO2SurfaceState
    with pytest.raises(ValueError, match="fingerprint mismatch"):
        advance_feature_step_3d(
            geometry, _boundary(),
            {"Ar+": "energetic_bombardment", "CF2": "neutral_reactant"},
            _mechanism(), etchable_material_ids=(1,), duration_s=1.0,
            source_bounds=(0.0, 0.75, 0.0, 0.75), source_z=1.75,
            surface_state=SiO2SurfaceState.bare((18,)), n_position=8,
            surface_state_mesh_fingerprint="not-the-current-mesh",
            reinitialize=False, transport_device="cpu")


def test_feature_step_accepts_surface_history_with_exact_current_mesh_fingerprint():
    geometry, _ = _plane_geometry()
    common = dict(
        geometry=geometry, boundary=_boundary(),
        species_role={"Ar+": "energetic_bombardment", "CF2": "neutral_reactant"},
        mechanism=_mechanism(), etchable_material_ids=(1,), duration_s=0.0,
        source_bounds=(0.0, 0.75, 0.0, 0.75), source_z=1.75,
        n_position=8, reinitialize=False, transport_device="cpu")
    first = advance_feature_step_3d(**common)
    replay = advance_feature_step_3d(
        **common, surface_state=first.surface.state,
        surface_state_mesh_fingerprint=first.surface_state_mesh_fingerprint)
    assert replay.surface_state_mesh_fingerprint == first.surface_state_mesh_fingerprint
    assert np.array_equal(replay.surface.state.complex_fraction, first.surface.state.complex_fraction)


def test_feature_step_routes_old_and_new_surface_extraction_through_uniform_backend(
        monkeypatch):
    geometry, _ = _plane_geometry()
    calls = []
    backend = feature_step_module.UniformFeatureGeometryBackend3D
    original = backend.extract_surface

    def capture(self):
        calls.append(self.geometry)
        return original(self)

    monkeypatch.setattr(backend, "extract_surface", capture)
    result = advance_feature_step_3d(
        geometry, _boundary(),
        {"Ar+": "energetic_bombardment", "CF2": "neutral_reactant"},
        _mechanism(), etchable_material_ids=(1,), duration_s=0.0,
        source_bounds=(0.0, 0.75, 0.0, 0.75), source_z=1.75,
        n_position=8, reinitialize=False, transport_device="cpu")

    assert len(calls) == 2
    assert calls[0] is geometry
    assert calls[1] is result.geometry
    assert np.array_equal(calls[0].phi, calls[1].phi)


def test_surface_remap_preserves_material_integrals_and_coverage_bounds():
    from petch.surface_kinetics import SiO2SurfaceState
    old_centroid = np.array([
        [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [1.0, 1.0, 0.0]])
    new_centroid = old_centroid + [0.0, 0.0, -0.1]
    old_area = np.array([1.0, 2.0, 1.5, 0.5])
    new_area = np.array([0.8, 2.2, 1.4, 0.6])
    material = np.array([1, 1, 2, 2])
    state = SiO2SurfaceState(
        [0.1, 0.8, 0.3, 0.6], [1e18, 2e18, 3e18, 4e18],
        [2e17, 4e17, 6e17, 8e17], [0.05, 0.4, 0.1, 0.3],
        [0.2, 0.6, 0.4, 0.8])
    remapped, diagnostics = conservative_remap_surface_state(
        state, old_centroid, old_area, material, new_centroid, new_area, material,
        dx=1.0, mesh_length_unit_m=1e-6)

    for material_id in (1, 2):
        old = material == material_id; new = material == material_id
        for before, after in (
                (state.complex_fraction, remapped.complex_fraction),
                (state.polymer_units_m2, remapped.polymer_units_m2),
                (state.removed_formula_units_m2, remapped.removed_formula_units_m2),
                (state.activated_complex_fraction, remapped.activated_complex_fraction),
                (state.activated_polymer_fraction, remapped.activated_polymer_fraction)):
            assert np.isclose(np.dot(before[old], old_area[old]),
                              np.dot(after[new], new_area[new]), rtol=2e-13)
    assert np.all((remapped.complex_fraction >= 0.0) & (remapped.complex_fraction <= 1.0))
    assert diagnostics["maximum_nearest_distance"] <= 0.1 + 1e-12


def test_surface_remap_certifies_interface_distance_not_retriangulated_centroids():
    from petch.surface_kinetics import SiO2SurfaceState

    old_triangle = np.array([[[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 2.0, 0.0]]])
    old_centroid = np.mean(old_triangle, axis=1)
    new_centroid = np.array([[1.8, 0.1, 0.0]])
    area = np.array([2.0])
    material = np.array([1])
    state = SiO2SurfaceState([0.2], [1e18], [2e17], [0.1], [0.4])

    with pytest.raises(ValueError, match="surface remap distance"):
        conservative_remap_surface_state(
            state, old_centroid, area, material,
            new_centroid, area, material,
            dx=1.0, mesh_length_unit_m=1e-6,
            maximum_distance=0.01)
    with pytest.raises(ValueError, match="point-to-triangle"):
        conservative_remap_surface_state(
            state, old_centroid, area, material,
            np.array([[1.8, 0.3, 0.0]]), area, material,
            dx=1.0, mesh_length_unit_m=1e-6,
            maximum_distance=0.01, old_triangles=old_triangle)

    remapped, diagnostics = conservative_remap_surface_state(
        state, old_centroid, area, material,
        new_centroid, area, material,
        dx=1.0, mesh_length_unit_m=1e-6,
        maximum_distance=0.01, old_triangles=old_triangle)

    assert np.array_equal(remapped.complex_fraction, state.complex_fraction)
    assert np.array_equal(remapped.polymer_units_m2, state.polymer_units_m2)
    assert diagnostics["distance_metric"] == "point_to_material_triangle"
    assert diagnostics["maximum_nearest_distance"] <= 1e-14
    assert diagnostics["maximum_nearest_centroid_distance"] > 0.01


def test_surface_remap_uses_wrapped_distance_across_periodic_seam():
    from petch.surface_kinetics import SiO2SurfaceState

    old_triangle = np.array([[[0.0, 0.98, 0.0], [0.1, 0.98, 0.0], [0.0, 1.0, 0.0]]])
    old_centroid = np.mean(old_triangle, axis=1)
    new_centroid = np.array([[0.02, 0.005, 0.0]])
    area = np.array([0.001])
    material = np.array([1])
    state = SiO2SurfaceState([0.2], [1e18], [2e17], [0.1], [0.4])

    with pytest.raises(ValueError, match="surface remap distance"):
        conservative_remap_surface_state(
            state, old_centroid, area, material,
            new_centroid, area, material,
            dx=0.01, mesh_length_unit_m=1e-6,
            maximum_distance=0.01, old_triangles=old_triangle)

    remapped, diagnostics = conservative_remap_surface_state(
        state, old_centroid, area, material,
        new_centroid, area, material,
        dx=0.01, mesh_length_unit_m=1e-6,
        maximum_distance=0.01, old_triangles=old_triangle,
        periodic_lengths=(None, 1.0, None))

    assert np.array_equal(remapped.polymer_units_m2, state.polymer_units_m2)
    assert diagnostics["periodic_lengths"] == (None, 1.0, None)
    assert diagnostics["maximum_nearest_distance"] <= 0.009


def test_surface_remap_is_exact_identity_on_unchanged_heterogeneous_mesh():
    from petch.surface_kinetics import SiO2SurfaceState
    centroid = np.array([
        [0.0, 0.0, 0.0], [0.5, 0.0, 0.0],
        [0.0, 0.5, 0.0], [0.5, 0.5, 0.0]])
    area = np.array([0.7, 1.1, 0.9, 1.3])
    material = np.ones(4, dtype=int)
    state = SiO2SurfaceState(
        [0.0, 0.2, 0.7, 1.0],
        [0.0, 2e18, 9e18, 4e19],
        [1e16, 3e17, 2e18, 8e18],
        [0.0, 0.1, 0.35, 0.7],
        [0.0, 0.2, 0.6, 1.0])

    remapped, diagnostics = conservative_remap_surface_state(
        state, centroid, area, material, centroid.copy(), area.copy(), material.copy(),
        dx=0.1, mesh_length_unit_m=1e-6)

    assert np.array_equal(remapped.complex_fraction, state.complex_fraction)
    assert np.array_equal(remapped.polymer_units_m2, state.polymer_units_m2)
    assert np.array_equal(remapped.removed_formula_units_m2, state.removed_formula_units_m2)
    assert np.array_equal(
        remapped.activated_complex_fraction, state.activated_complex_fraction)
    assert np.array_equal(
        remapped.activated_polymer_fraction, state.activated_polymer_fraction)
    assert diagnostics["maximum_nearest_distance"] == 0.0


def test_multistep_solver_carries_remapped_state_and_matches_planar_total_motion():
    geometry, initial_height = _plane_geometry()
    result = solve_feature_3d(
        geometry, _boundary(),
        {"Ar+": "energetic_bombardment", "CF2": "neutral_reactant"},
        _mechanism(), etchable_material_ids=(1,), duration_s=2.0, n_steps=2,
        source_bounds=(0.0, 0.75, 0.0, 0.75), source_z=1.75,
        n_position=16384, seed=13, cfl_number=0.3, reinitialize=False,
        transport_device="cpu")

    final_height = _area_weighted_height(result.geometry.phi, geometry.dx)
    assert np.isclose(initial_height - final_height, 0.04, atol=0.004)
    assert len(result.steps) == 2
    assert result.steps[0].next_surface_state_mesh_fingerprint == (
        result.steps[1].surface_state_mesh_fingerprint)
    assert result.surface_state_mesh_fingerprint == (
        result.steps[-1].next_surface_state_mesh_fingerprint)
    assert result.validity.within_declared_scope


def test_feature_step_uses_supplied_3d_potential_for_ion_energy_and_surface_velocity():
    geometry, initial_height = _plane_geometry()
    z = np.arange(geometry.phi.shape[2]) * geometry.dx
    potential = np.broadcast_to(10.0 * z / 1.75, geometry.phi.shape).copy()
    result = advance_feature_step_3d(
        geometry, _boundary(),
        {"Ar+": "energetic_bombardment", "CF2": "neutral_reactant"},
        _mechanism(), etchable_material_ids=(1,), duration_s=1.0,
        source_bounds=(0.0, 0.75, 0.0, 0.75), source_z=1.75,
        nodal_potential_v=potential, potential_origin=(0.0, 0.0, 0.0),
        potential_spacing=geometry.dx, trajectory_fixed_dt=0.005,
        trajectory_max_steps=1000, n_position=16384, seed=37,
        cfl_number=0.3, reinitialize=False, transport_device="cpu")

    expected_energy = 100.0 + 10.0 * (1.0 - initial_height / 1.75)
    expected_yield = 0.2 * (expected_energy - 20.0) / (100.0 - 20.0)
    expected_velocity = 2.2e21 * expected_yield / 2.2e28
    average_velocity = np.dot(
        result.surface.etch_velocity_m_s, result.active_face_area) / result.active_face_area.sum()
    assert result.transport.transport_model == "collisionless_fixed_step_nodal_field_3d"
    assert np.isclose(average_velocity, expected_velocity, rtol=3e-4)
    assert _area_weighted_height(result.geometry.phi, geometry.dx) < initial_height - 0.02


def test_feature_step_consumes_a_precomputed_exact_transport_without_retracing():
    geometry, initial_height = _plane_geometry()
    boundary = _boundary()
    role = {"Ar+": "energetic_bombardment", "CF2": "neutral_reactant"}
    verts, faces, centroids, areas = extract_mesh_3d(geometry.phi, geometry.dx)
    normals = _surface_gas_normals(verts, faces, centroids, geometry)
    potential = np.zeros(geometry.phi.shape)
    transport = trace_boundary_state_field_3d(
        boundary, role, verts, faces, areas,
        source_bounds=(0.0, 0.75, 0.0, 0.75), source_z=1.75,
        nodal_potential_v=potential, potential_origin=(0.0, 0.0, 0.0),
        potential_spacing=geometry.dx, mesh_length_unit_m=geometry.mesh_length_unit_m,
        n_position=1024, seed=37, fixed_dt=0.005, max_steps=1000,
        face_gas_normals=normals, device="cpu")
    mechanism = _mechanism()
    result = advance_feature_step_3d(
        geometry, boundary, role, mechanism,
        etchable_material_ids=(1,), duration_s=1.0,
        source_bounds=(0.0, 0.75, 0.0, 0.75), source_z=1.75,
        precomputed_transport=transport, n_position=1, seed=999,
        cfl_number=0.3, reinitialize=False, transport_device="cpu")

    assert result.transport is transport
    assert result.charging is None
    ion = next(item for item in transport.surface_fluxes.energetic_fluxes
               if item.name == "Ar+")
    expected_velocity = (
        ion.yield_rate_m2_s(mechanism.parameters.bare_sio2_yield)
        [result.active_face_index] / mechanism.parameters.bulk_formula_density_m3)
    assert np.allclose(result.surface.etch_velocity_m_s, expected_velocity, rtol=2e-13)
    assert _area_weighted_height(result.geometry.phi, geometry.dx) < initial_height


def test_feature_step_solves_charge_reuses_ion_events_and_excludes_electron_from_chemistry():
    geometry, initial_height = _plane_geometry()
    mechanism = _mechanism()
    result = advance_feature_step_3d(
        geometry, _charging_boundary(),
        {"Ar+": "energetic_bombardment", "electron": "charge_carrier"},
        mechanism, etchable_material_ids=(1,), duration_s=1.0,
        source_bounds=(0.0, 0.75, 0.0, 0.75), source_z=1.75,
        potential_origin=(0.0, 0.0, 0.0), potential_spacing=geometry.dx,
        trajectory_fixed_dt=0.005, trajectory_max_steps=2000,
        charging_poisson_system=_plane_poisson_system(geometry),
        charging_options=dict(
            max_iter=30, min_iter=2, current_balance_tol=1e-12,
            beta=0.5, response_energy_eV=4.0),
        n_position=64, seed=61, cfl_number=0.3, reinitialize=False,
        transport_device="cpu")

    assert result.charging is not None and result.charging.converged
    assert result.transport is result.charging.transport
    support = (result.charging.positive_current_node_a
               + result.charging.negative_current_node_a) > 0.0
    assert np.allclose(
        result.charging.positive_current_node_a[support],
        result.charging.negative_current_node_a[support], rtol=1e-14)
    surface_voltage = result.charging.potential_v[:, :, 4]
    assert np.all((-20.0 < surface_voltage) & (surface_voltage < -1.0))

    populations = {item.name: item for item in result.transport.surface_fluxes.energetic_fluxes}
    assert set(populations) == {"Ar+", "electron"}
    ion_only_velocity = (
        populations["Ar+"].yield_rate_m2_s(mechanism.parameters.bare_sio2_yield)
        [result.active_face_index] / mechanism.parameters.bulk_formula_density_m3)
    assert np.allclose(result.surface.etch_velocity_m_s, ion_only_velocity, rtol=2e-13)
    assert result.diagnostics["self_consistent_charging"]
    assert result.diagnostics["charging_converged"]
    assert _area_weighted_height(result.geometry.phi, geometry.dx) < initial_height - 0.02


def test_periodic_trench_forward_and_adjoint_electron_currents_close_by_region():
    geometry = make_rectangular_trench_geometry_3d(
        cell_width=1.0, cell_length=0.5, domain_height=2.0, dx=0.25,
        opening_width=0.5, mask_thickness=0.25,
        substrate_top=1.25, etched_depth=0.75)
    source_z = 2.0
    boundary = _continuous_trench_charging_boundary(
        source_z * geometry.mesh_length_unit_m)
    verts, faces, centroids, areas = extract_mesh_3d(geometry.phi, geometry.dx)
    normals = _surface_gas_normals(verts, faces, centroids, geometry)
    electron_boundary = PlasmaBoundaryState(
        (boundary.get("electron"),), boundary.reference_plane_m)
    common = dict(
        boundary=electron_boundary, species_role={"electron": "charge_carrier"},
        verts=verts, faces=faces, areas=areas,
        source_bounds=(0.0, 1.0, 0.0, 0.5), source_z=source_z,
        nodal_potential_v=np.zeros(geometry.phi.shape),
        potential_origin=(0.0, 0.0, 0.0), potential_spacing=geometry.dx,
        mesh_length_unit_m=geometry.mesh_length_unit_m,
        fixed_dt=0.005, max_steps=50000, periodic_lateral=True, device="cpu")
    forward = trace_boundary_state_field_3d(
        **common, phase_space_log2_samples=14, face_gas_normals=normals)
    adjoint = gather_boundary_state_field_adjoint_3d(
        **common, centroids=centroids, gas_normals=normals,
        face_quadrature_points=3, ray_offset=1e-4,
        proposal_by_species={"electron": _trench_adjoint_proposals(boundary)["electron"]})
    physical_flux = boundary.get("electron").flux_m2_s
    forward_flux = forward.surface_fluxes.energetic_fluxes[0].flux_m2_s / physical_flux
    adjoint_flux = adjoint.surface_fluxes.energetic_fluxes[0].flux_m2_s / physical_flux
    regions = (
        centroids[:, 2] < 0.75,
        (centroids[:, 2] >= 0.75) & (centroids[:, 2] < 1.25),
        (centroids[:, 2] >= 1.25) & (centroids[:, 2] < 1.45),
        centroids[:, 2] >= 1.45,
    )

    assert np.isclose(forward.hit_probability["electron"], 1.0, atol=1.0 / 2 ** 14)
    assert np.isclose(adjoint.hit_probability["electron"], 1.0, rtol=2e-3)
    for region in regions:
        forward_measure = np.dot(forward_flux[region], areas[region]) / 0.5
        adjoint_measure = np.dot(adjoint_flux[region], areas[region]) / 0.5
        assert np.isclose(adjoint_measure, forward_measure, atol=0.012)


def test_periodic_trench_source_aligned_adjoint_resolves_directional_ion_current():
    geometry = make_rectangular_trench_geometry_3d(
        cell_width=1.0, cell_length=0.5, domain_height=2.0, dx=0.25,
        opening_width=0.5, mask_thickness=0.25,
        substrate_top=1.25, etched_depth=0.75)
    source_z = 2.0
    boundary = _continuous_trench_charging_boundary(
        source_z * geometry.mesh_length_unit_m)
    verts, faces, centroids, areas = extract_mesh_3d(geometry.phi, geometry.dx)
    normals = _surface_gas_normals(verts, faces, centroids, geometry)
    ion_boundary = PlasmaBoundaryState((boundary.get("Ar+"),), boundary.reference_plane_m)
    common = dict(
        boundary=ion_boundary, species_role={"Ar+": "energetic_bombardment"},
        verts=verts, faces=faces, areas=areas,
        source_bounds=(0.0, 1.0, 0.0, 0.5), source_z=source_z,
        nodal_potential_v=np.zeros(geometry.phi.shape),
        potential_origin=(0.0, 0.0, 0.0), potential_spacing=geometry.dx,
        mesh_length_unit_m=geometry.mesh_length_unit_m,
        fixed_dt=0.005, max_steps=50000, periodic_lateral=True, device="cpu")
    forward = trace_boundary_state_field_3d(
        **common, phase_space_log2_samples=14, face_gas_normals=normals)
    proposal = {"Ar+": qmc_boundary_proposal(boundary.get("Ar+"), 10, seed=79)}
    coarse_adjoint = gather_boundary_state_field_adjoint_3d(
        **common, centroids=centroids, gas_normals=normals,
        face_quadrature_points=3, ray_offset=1e-4,
        proposal_by_species=proposal,
        proposal_frame_by_species={"Ar+": "source_aligned"})
    adjoint = gather_boundary_state_field_adjoint_3d(
        **common, centroids=centroids, gas_normals=normals,
        face_quadrature_points=7, ray_offset=1e-4,
        proposal_by_species=proposal,
        proposal_frame_by_species={"Ar+": "source_aligned"})
    physical_flux = boundary.get("Ar+").flux_m2_s
    forward_flux = forward.surface_fluxes.energetic_fluxes[0].flux_m2_s / physical_flux
    adjoint_flux = adjoint.surface_fluxes.energetic_fluxes[0].flux_m2_s / physical_flux
    regions = (
        centroids[:, 2] < 0.75,
        (centroids[:, 2] >= 0.75) & (centroids[:, 2] < 1.25),
        (centroids[:, 2] >= 1.25) & (centroids[:, 2] < 1.45),
        centroids[:, 2] >= 1.45,
    )

    assert np.isclose(forward.hit_probability["Ar+"], 1.0, atol=1.0 / 2 ** 14)
    assert abs(adjoint.hit_probability["Ar+"] - 1.0) < abs(
        coarse_adjoint.hit_probability["Ar+"] - 1.0)
    assert np.isclose(adjoint.hit_probability["Ar+"], 1.0, rtol=6e-3)
    for region in regions:
        forward_measure = np.dot(forward_flux[region], areas[region]) / 0.5
        adjoint_measure = np.dot(adjoint_flux[region], areas[region]) / 0.5
        assert np.isclose(adjoint_measure, forward_measure, atol=0.012)


def test_multistep_charged_profile_refuses_a_fixed_geometry_poisson_operator():
    geometry, _ = _plane_geometry()
    with pytest.raises(ValueError, match="geometry-dependent Poisson builder"):
        solve_feature_3d(
            geometry, _charging_boundary(),
            {"Ar+": "energetic_bombardment", "electron": "charge_carrier"},
            _mechanism(), etchable_material_ids=(1,), duration_s=2.0, n_steps=2,
            source_bounds=(0.0, 0.75, 0.0, 0.75), source_z=1.75,
            potential_origin=(0.0, 0.0, 0.0), potential_spacing=geometry.dx,
            trajectory_fixed_dt=0.005,
            charging_poisson_system=_plane_poisson_system(geometry))


def test_multistep_quasistatic_charging_rebuilds_material_operator_and_reconverges():
    geometry, initial_height = _plane_geometry(); systems = []

    def build(current_geometry):
        system = _plane_poisson_system(current_geometry); systems.append(system)
        return system

    result = solve_feature_3d(
        geometry, _charging_boundary(),
        {"Ar+": "energetic_bombardment", "electron": "charge_carrier"},
        _mechanism(), etchable_material_ids=(1,), duration_s=10.0, n_steps=2,
        source_bounds=(0.0, 0.75, 0.0, 0.75), source_z=1.75,
        potential_origin=(0.0, 0.0, 0.0), potential_spacing=geometry.dx,
        trajectory_fixed_dt=0.005, trajectory_max_steps=2000,
        charging_system_builder=build,
        charging_options=dict(
            max_iter=30, min_iter=2, current_balance_tol=1e-12,
            beta=0.5, response_energy_eV=4.0),
        n_position=64, seed=67, cfl_number=0.3, reinitialize=False,
        transport_device="cpu")

    assert len(systems) == 2 and len(result.steps) == 2
    assert all(step.charging is not None and step.charging.converged for step in result.steps)
    assert np.count_nonzero(systems[1].epsilon_r == 3.9) < np.count_nonzero(
        systems[0].epsilon_r == 3.9)
    assert _area_weighted_height(result.geometry.phi, geometry.dx) < initial_height - 0.15
    assert "quasi-static charging" in " ".join(result.validity.known_limitations)


def test_feature_step_refuses_nonconverged_charging_for_profile_motion():
    geometry, _ = _plane_geometry()
    with pytest.raises(ValueError, match="requires a converged"):
        advance_feature_step_3d(
            geometry, _charging_boundary(),
            {"Ar+": "energetic_bombardment", "electron": "charge_carrier"},
            _mechanism(), etchable_material_ids=(1,), duration_s=1.0,
            source_bounds=(0.0, 0.75, 0.0, 0.75), source_z=1.75,
            potential_origin=(0.0, 0.0, 0.0), potential_spacing=geometry.dx,
            trajectory_fixed_dt=0.005,
            charging_poisson_system=_plane_poisson_system(geometry),
            charging_options={"require_converged": False})


def test_second_chemistry_runs_through_unchanged_transport_remap_and_interface_engine():
    geometry, initial_height = _plane_geometry(); mechanism = _si_cl_ar_mechanism()
    result = solve_feature_3d(
        geometry, _si_cl_ar_boundary(),
        {"Ar+": "energetic_bombardment", "Cl2": "neutral_reactant"},
        mechanism, etchable_material_ids=(1,), duration_s=4.0, n_steps=2,
        source_bounds=(0.0, 0.75, 0.0, 0.75), source_z=1.75,
        n_position=4096, seed=73, cfl_number=0.3, reinitialize=False,
        transport_device="cpu")

    expected_velocity_m_s = 2e21 * 0.24182079610957588 / SI_ATOM_DENSITY_M3
    mean_velocity = np.mean(result.steps[0].surface.etch_velocity_m_s)
    assert isinstance(result.surface_state, TabulatedSiSurfaceState)
    assert np.isclose(mean_velocity, expected_velocity_m_s, rtol=0.01)
    assert all(step.surface.table_fingerprint == mechanism.table.fingerprint
               for step in result.steps)
    assert result.validity.parameter_evidence_supports_prediction
    assert result.validity.nonpredictive_parameters == ()
    assert result.surface_state.removed_atoms_m2.size == result.steps[-1].next_active_face_area.size
    assert _area_weighted_height(result.geometry.phi, geometry.dx) < initial_height - 0.03
