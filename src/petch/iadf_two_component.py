"""Two-component (thermal core + collisional tail) ion angular/energy distribution.

Every measured sub-degree ion angular distribution at high-aspect-ratio etch
conditions is **bi-Gaussian**: a narrow thermal core carrying the sheath-edge
transverse temperature, plus a much wider non-thermal tail produced by
ion--neutral collisions *inside* the sheath.  A collisionless (single-Gaussian)
boundary cannot express the tail at any quadrature order, and at AR >= 100 the
tail carries essentially all of the sidewall flux --- see
``RESEARCH_IADF_SUBDEGREE_AND_REACTOR_2026-07-29.md`` sections A.3 and A.7.

Model
-----
At a fixed normal (impact) energy ``E`` the transverse velocity components are
two independent zero-mean Gaussians per component, so with ``v_z = sqrt(E)`` in
the engine's sqrt-eV convention::

    tan(theta_x), tan(theta_y)  ~  N(0, s^2)  independently,
    s(E) = sqrt(T_perp / (2 E))                                          (1)

``s`` is the **1-D projected (planar) angular width** in radians --- the width a
signed-angle IEAD plot shows.  It is mass independent (M cancels between the
transverse thermal speed and the sheath-accelerated normal speed), which is why
Ar+ and Kr+ are measured with the same main-component width.

Two identities follow exactly and are gated in the tests:

* **Polar vs planar.**  ``E[tan^2 theta] = 2 s^2`` while ``E[tan^2 theta_x] =
  s^2``: the polar rms exceeds the planar rms by exactly ``sqrt(2)``, for any
  axisymmetric measure.  This is the same identity :mod:`petch.angular_lift`
  inverts when it lifts a *published* planar marginal; here the distribution is
  generated analytically, so the lift is unnecessary and the factor is built in.
* **Characteristic (2-D radial) spread.**  ``theta_th = sqrt(T_perp / E) =
  s * sqrt(2)`` --- the quantity Khrabrov & Kaganovich quote as
  "``theta_th = sqrt(T_perp/E_b)`` rad or approximately 0.4 degrees".

Acceptance
----------
The AR-200 observable is the fraction of the incident population inside the
geometric acceptance of a straight feature, ``alpha = arctan(1/AR)``.  Under (1)
both conventions integrate in closed form, with **no quadrature noise at any
aspect ratio**:

* round hole (axisymmetric cone, the physical acceptance of a via)::

      P(theta <= alpha) = sum_i f_i * (1 - exp(-tan^2(alpha) / (2 s_i^2)))   (2)

* planar/trench projection (the convention of the section A.7 table, i.e. the
  signed-angle marginal)::

      P(|theta_x| <= alpha) = sum_i f_i * erf(tan(alpha) / (s_i * sqrt(2)))   (3)

(2) is *not* (3); a round hole rejects more flux than the planar projection
suggests because both transverse components must be small simultaneously.  Both
are exposed, both are gated against numerical quadrature, and callers must pick
the one matching their geometry.

Provenance
----------
The three parameters are **declared, never fitted**.  The reference set is the
Nagoya/Kioxia measurement (better than 0.1 degree angular resolution):

    core T_perp = 0.044 eV, tail T_perp = 0.57 eV
    Kim, Kawamura, Naito, Iino, Fukumizu, Kurihara, Suzuki, Toyoda,
    Jpn. J. Appl. Phys. 64, 05SP15 (2025), doi:10.35848/1347-4065/adce84

The measured *tail fraction* is not printed in the accessible text and is
therefore carried as an explicit ``[VERIFY]`` in the provenance of the reference
set: ``KIM_2025_TAIL_FRACTION_VERIFY``.  The physically-derived alternative is
the sheath collided fraction ``1 - exp(-s/lambda)``, which for the Krueger base
case (s = 14 mm, lambda = 13.8 mm) is 0.64 and independently matches the 0.65
recovered by fitting the digitized Krueger Figure-4 marginal.  Deriving that
fraction from a collision operator instead of declaring it is stage S2 of the
research plan; this module is stage S1 and only *represents* it.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from types import MappingProxyType
from typing import Mapping

import numpy as np
from scipy.special import erf

__all__ = [
    "AngularComponent",
    "TwoComponentIADF",
    "KIM_2025_CORE_TEMPERATURE_EV",
    "KIM_2025_TAIL_TEMPERATURE_EV",
    "KIM_2025_TAIL_FRACTION_VERIFY",
    "KIM_2025_PROVENANCE",
    "KRUEGER_2024_FIGURE4_BIGAUSSIAN",
    "acceptance_half_angle_deg",
    "build_two_component_boundary",
]


# --- Declared parameter sets -------------------------------------------------

#: Measured transverse temperature of the thermal core (Kim 2025).  Equals the
#: neutral gas temperature (0.044 eV = 511 K) and is reported as
#: condition-independent, which is the signature of a thermodynamic property
#: rather than a fitted width.
KIM_2025_CORE_TEMPERATURE_EV = 0.044

#: Measured transverse temperature of the collisional tail (Kim 2025): 13x
#: hotter than the core, hence sqrt(13) = 3.6x wider.
KIM_2025_TAIL_TEMPERATURE_EV = 0.57

#: [VERIFY] The tail *fraction* is not printed in the accessible text.  0.65 is
#: the value recovered by fitting the repo's digitized Krueger Figure-4 marginal
#: and independently equals the sheath collided fraction 1 - exp(-s/lambda)
#: = 0.64 at s = 14 mm, lambda = 13.8 mm.  Declared, not measured.
KIM_2025_TAIL_FRACTION_VERIFY = 0.65

KIM_2025_PROVENANCE = MappingProxyType({
    "core_temperature_eV": "Kim et al., Jpn. J. Appl. Phys. 64, 05SP15 (2025), "
                           "doi:10.35848/1347-4065/adce84 -- effective ion "
                           "temperature, main component",
    "tail_temperature_eV": "Kim et al. 2025 (same) -- tail component; mechanism "
                           "identified as Ar+ collisions with Ar in the sheath, "
                           "Jpn. J. Appl. Phys. 64, 096002 (2025), "
                           "doi:10.35848/1347-4065/ae0105",
    "tail_fraction": "[VERIFY] not printed in accessible text; declared from the "
                     "digitized Krueger Figure-4 bi-Gaussian fit (0.65) which "
                     "matches the collided fraction 1-exp(-s/lambda)=0.64",
    "measurement_conditions": "dual-frequency CCP, Ar, 2.4 Pa, V_pp 2.7 kV, "
                              "V_dc 950 V, ion energy 1.4-2.0 keV, MCP imaging "
                              "at better than 0.1 degree resolution",
    "evidence_kind": "measured",
})

#: Widths (1-D planar, degrees) of the bi-Gaussian fitted to the repo's
#: digitized Krueger Figure-4 marginal, with its declared limitation: the
#: digitization grid is 0.25 degrees, so the fitted core width is an upper
#: bound (the true core cannot be resolved on that grid).
KRUEGER_2024_FIGURE4_BIGAUSSIAN = MappingProxyType({
    "core_sigma_planar_deg": 0.600,
    "tail_sigma_planar_deg": 0.946,
    "tail_fraction": 0.65,
    "marginal_sigma_band_deg": (0.822, 0.860),
    "evidence_kind": "digitized",
    "limitation": "0.25 degree digitization grid; fitted core width is an upper "
                  "bound and the marginal sigma is the only gradeable statistic",
})


def acceptance_half_angle_deg(aspect_ratio):
    """Geometric acceptance half-angle ``arctan(1/AR)`` of a straight feature."""
    aspect = np.asarray(aspect_ratio, dtype=float)
    if np.any(aspect <= 0.0) or np.any(~np.isfinite(aspect)):
        raise ValueError("aspect ratio must be positive and finite")
    return np.rad2deg(np.arctan(1.0 / aspect))


# --- Components --------------------------------------------------------------

@dataclass(frozen=True)
class AngularComponent:
    """One axisymmetric Gaussian component of the angular distribution.

    Exactly one width source must be given: ``temperature_eV`` (physical, width
    then varies as ``E**-0.5``) or ``sigma_planar_deg`` (a fixed digitized width,
    energy independent).  ``fraction`` is the flux share.
    """

    fraction: float
    temperature_eV: float | None = None
    sigma_planar_deg: float | None = None
    label: str = ""

    def __post_init__(self):
        fraction = float(self.fraction)
        if not np.isfinite(fraction) or fraction < 0.0:
            raise ValueError("component fraction must be finite and nonnegative")
        has_temperature = self.temperature_eV is not None
        has_sigma = self.sigma_planar_deg is not None
        if has_temperature == has_sigma:
            raise ValueError(
                "give exactly one of temperature_eV or sigma_planar_deg")
        if has_temperature:
            temperature = float(self.temperature_eV)
            if not np.isfinite(temperature) or temperature <= 0.0:
                raise ValueError("temperature must be positive and finite")
            object.__setattr__(self, "temperature_eV", temperature)
        else:
            sigma = float(self.sigma_planar_deg)
            if not np.isfinite(sigma) or sigma <= 0.0:
                raise ValueError("sigma must be positive and finite")
            object.__setattr__(self, "sigma_planar_deg", sigma)
        object.__setattr__(self, "fraction", fraction)

    def sigma_planar_rad(self, energy_eV):
        """1-D projected width ``s`` of equation (1), in radians."""
        energy = np.asarray(energy_eV, dtype=float)
        if np.any(energy <= 0.0) or np.any(~np.isfinite(energy)):
            raise ValueError("energy must be positive and finite")
        if self.temperature_eV is not None:
            return np.sqrt(self.temperature_eV / (2.0 * energy))
        return np.full(np.shape(energy), np.deg2rad(self.sigma_planar_deg))


@dataclass(frozen=True)
class TwoComponentIADF:
    """Weighted sum of axisymmetric Gaussian angular components.

    The name records the physics (core + collisional tail); the implementation
    accepts any number of components so that the single-component limit and
    future multi-collision decompositions are expressible without a second type.
    """

    components: tuple[AngularComponent, ...]
    provenance: Mapping[str, object] = field(default_factory=dict)

    def __post_init__(self):
        components = tuple(self.components)
        if not components:
            raise ValueError("at least one angular component is required")
        total = float(sum(component.fraction for component in components))
        if not np.isfinite(total) or total <= 0.0:
            raise ValueError("component fractions must carry positive mass")
        if abs(total - 1.0) > 1e-12:
            raise ValueError(
                f"component fractions must sum to 1 (got {total!r}); the split "
                "is declared physics, never normalized silently")
        if not self.provenance:
            raise ValueError(
                "an IADF must carry provenance for its declared parameters")
        object.__setattr__(self, "components", components)
        object.__setattr__(self, "provenance", MappingProxyType(dict(self.provenance)))

    # -- widths --------------------------------------------------------------

    def sigma_planar_deg(self, energy_eV):
        """Per-component 1-D projected widths in degrees, shape ``(n_components,)``."""
        return np.array([
            np.rad2deg(component.sigma_planar_rad(energy_eV))
            for component in self.components])

    def theta_thermal_deg(self, energy_eV):
        """Per-component 2-D radial spread ``sqrt(T/E) = sqrt(2) * s``, degrees."""
        return np.sqrt(2.0) * self.sigma_planar_deg(energy_eV)

    def marginal_sigma_planar_deg(self, energy_eV):
        """rms of the *combined* planar (signed-angle) marginal, in degrees.

        This is the statistic a digitized signed-angle IEAD can be graded on.
        """
        sigma = self.sigma_planar_deg(energy_eV)
        fractions = np.array([component.fraction for component in self.components])
        weighted = np.tensordot(fractions, sigma ** 2, axes=(0, 0))
        return np.sqrt(weighted)

    def polar_rms_deg(self, energy_eV):
        """rms polar angle of the combined measure: ``sqrt(2)`` x planar rms."""
        return np.sqrt(2.0) * self.marginal_sigma_planar_deg(energy_eV)

    # -- acceptance ----------------------------------------------------------

    def acceptance_fraction_cone(self, half_angle_deg, energy_eV):
        """Flux fraction inside a cone of half-angle ``alpha`` -- equation (2).

        The physical acceptance of a round hole.  Exact under the model; no
        quadrature is involved at any aspect ratio.
        """
        tangent = self._acceptance_tangent(half_angle_deg)
        total = 0.0
        for component in self.components:
            sigma = component.sigma_planar_rad(energy_eV)
            total = total + component.fraction * (
                1.0 - np.exp(-0.5 * (tangent / sigma) ** 2))
        return total

    def acceptance_fraction_planar(self, half_angle_deg, energy_eV):
        """Fraction with ``|theta_x| <= alpha`` in one plane -- equation (3).

        The trench/signed-angle-marginal convention (and the convention of the
        research document's section A.7 acceptance table).
        """
        tangent = self._acceptance_tangent(half_angle_deg)
        total = 0.0
        for component in self.components:
            sigma = component.sigma_planar_rad(energy_eV)
            total = total + component.fraction * erf(
                tangent / (sigma * np.sqrt(2.0)))
        return total

    @staticmethod
    def _acceptance_tangent(half_angle_deg):
        angle = np.asarray(half_angle_deg, dtype=float)
        if np.any(angle <= 0.0) or np.any(angle >= 90.0) or np.any(~np.isfinite(angle)):
            raise ValueError("acceptance half-angle must lie in (0, 90) degrees")
        return np.tan(np.deg2rad(angle))

    # -- quadrature export ---------------------------------------------------

    def polar_quadrature(self, energy_eV, *, n_polar=64, max_sigma=6.0):
        """Polar nodes and weights of the combined measure at one energy.

        Nodes are midpoints of ``n_polar`` equal-probability-free bins spanning
        ``[0, max_sigma * widest sigma]``; each node carries the exact analytic
        cone mass of its shell, so the returned weights are the measure itself
        (not a sampled approximation) and sum to the mass inside the outer
        radius.  The residual beyond ``max_sigma`` is folded into the outermost
        node so the export conserves flux to machine precision.
        """
        if int(n_polar) < 1:
            raise ValueError("n_polar must be positive")
        sigma = np.array([
            float(np.asarray(component.sigma_planar_rad(energy_eV)).reshape(()))
            for component in self.components])
        outer = float(max_sigma) * float(sigma.max()) * np.sqrt(2.0)
        if not np.isfinite(outer) or outer <= 0.0:
            raise ValueError("quadrature support degenerated")
        edges = np.linspace(0.0, outer, int(n_polar) + 1)
        cumulative = self._cone_mass_at_tangent(np.tan(edges), sigma)
        mass = np.diff(cumulative)
        mass[-1] += 1.0 - cumulative[-1]
        centres = 0.5 * (edges[:-1] + edges[1:])
        return np.rad2deg(centres), mass

    def _cone_mass_at_tangent(self, tangent, sigma):
        total = np.zeros(np.shape(tangent), dtype=float)
        for component, width in zip(self.components, sigma):
            total = total + component.fraction * (
                1.0 - np.exp(-0.5 * (tangent / width) ** 2))
        return total

    def discrete_nodes(self, energy_eV, energy_weight=None, *, n_polar=64,
                       azimuthal_order=16, max_sigma=6.0):
        """Discrete ``(energy, polar_deg, weight)`` nodes for the boundary state.

        Azimuth is closed uniformly at ``azimuthal_order`` nodes exactly as the
        axisymmetric closure in :mod:`petch.reactor_boundary` does; the polar
        weight is divided evenly among them by the caller-facing
        :func:`build_two_component_boundary`.  Returned weights sum to 1.
        """
        energies = np.atleast_1d(np.asarray(energy_eV, dtype=float))
        if energy_weight is None:
            energy_weight = np.full(energies.shape, 1.0 / energies.size)
        energy_weight = np.asarray(energy_weight, dtype=float)
        if energy_weight.shape != energies.shape:
            raise ValueError("energy weights must match energies")
        if np.any(energy_weight < 0.0) or energy_weight.sum() <= 0.0:
            raise ValueError("energy weights must be nonnegative with positive mass")
        energy_weight = energy_weight / energy_weight.sum()
        if int(azimuthal_order) < 1:
            raise ValueError("azimuthal_order must be positive")
        out_energy, out_polar, out_weight = [], [], []
        for value, share in zip(energies, energy_weight):
            polar_deg, mass = self.polar_quadrature(
                value, n_polar=n_polar, max_sigma=max_sigma)
            out_energy.append(np.full(polar_deg.shape, value))
            out_polar.append(polar_deg)
            out_weight.append(share * mass)
        return (np.concatenate(out_energy), np.concatenate(out_polar),
                np.concatenate(out_weight))


def build_two_component_boundary(
        iadf, flux_m2_s, energy_eV, *, energy_weight=None, ion_mass_amu=40.0,
        name="ions", n_polar=64, azimuthal_order=16, max_sigma=6.0,
        reference_plane_m=0.0, extra_provenance=None):
    """Opt-in axisymmetric boundary state driven by a :class:`TwoComponentIADF`.

    Returns a :class:`~petch.boundary_state.PlasmaBoundaryState` whose single ion
    species carries the polar/azimuth quadrature of ``iadf`` and a matching
    :class:`~petch.boundary_state.DiscreteEnergyPolarAzimuthDensity3D`.  Nothing
    in the existing pipelines calls this; it is the S1 entry point for the
    axisymmetric hole path.
    """
    from .boundary_state import (
        DiscreteEnergyPolarAzimuthDensity3D,
        PlasmaBoundaryState,
        SpeciesBoundaryState,
    )

    energy, polar_deg, weight = iadf.discrete_nodes(
        energy_eV, energy_weight, n_polar=n_polar,
        azimuthal_order=azimuthal_order, max_sigma=max_sigma)
    order = int(azimuthal_order)
    azimuth = 2.0 * np.pi * (np.arange(order, dtype=float) + 0.5) / order
    speed = np.sqrt(energy)
    polar = np.deg2rad(polar_deg)
    transverse = (speed * np.sin(polar))[:, None]
    velocity = np.column_stack((
        (transverse * np.cos(azimuth)[None, :]).ravel(),
        (transverse * np.sin(azimuth)[None, :]).ravel(),
        np.repeat(speed * np.cos(polar), order),
    ))
    quadrature_weight = np.repeat(weight / order, order)
    provenance = {
        "provider": "two_component_iadf",
        "model": "bi_gaussian_core_plus_collisional_tail",
        "three_dimensional_azimuthal_closure":
            "axisymmetric uniform azimuth; the angular measure is generated "
            "analytically, so no planar-marginal inversion is required",
        "three_dimensional_azimuthal_order": order,
        "polar_node_count": int(n_polar),
        "iadf_provenance": dict(iadf.provenance),
        "component_fractions": tuple(
            component.fraction for component in iadf.components),
    }
    if extra_provenance:
        provenance.update(dict(extra_provenance))
    ion = SpeciesBoundaryState(
        name=name, charge_number=1, mass_amu=float(ion_mass_amu),
        flux_m2_s=float(flux_m2_s), velocity_sqrt_eV=velocity,
        weight=quadrature_weight,
        density_model=DiscreteEnergyPolarAzimuthDensity3D(
            energy, polar_deg, weight),
        provenance=provenance)
    return PlasmaBoundaryState(
        species=(ion,), reference_plane_m=float(reference_plane_m),
        provenance={"source": "TwoComponentIADF"})


def kim_2025_reference_iadf(tail_fraction=None):
    """The declared measured reference set (core 0.044 eV, tail 0.57 eV)."""
    fraction = (KIM_2025_TAIL_FRACTION_VERIFY if tail_fraction is None
                else float(tail_fraction))
    if not 0.0 <= fraction <= 1.0:
        raise ValueError("tail fraction must lie in [0, 1]")
    provenance = dict(KIM_2025_PROVENANCE)
    provenance["tail_fraction_value"] = fraction
    return TwoComponentIADF(
        components=(
            AngularComponent(fraction=1.0 - fraction,
                             temperature_eV=KIM_2025_CORE_TEMPERATURE_EV,
                             label="thermal_core"),
            AngularComponent(fraction=fraction,
                             temperature_eV=KIM_2025_TAIL_TEMPERATURE_EV,
                             label="collisional_tail"),
        ),
        provenance=provenance)


def krueger_2024_figure4_iadf():
    """The bi-Gaussian fitted to the repo's digitized Krueger Figure-4 marginal."""
    fraction = float(KRUEGER_2024_FIGURE4_BIGAUSSIAN["tail_fraction"])
    return TwoComponentIADF(
        components=(
            AngularComponent(
                fraction=1.0 - fraction,
                sigma_planar_deg=float(
                    KRUEGER_2024_FIGURE4_BIGAUSSIAN["core_sigma_planar_deg"]),
                label="fitted_core"),
            AngularComponent(
                fraction=fraction,
                sigma_planar_deg=float(
                    KRUEGER_2024_FIGURE4_BIGAUSSIAN["tail_sigma_planar_deg"]),
                label="fitted_tail"),
        ),
        provenance=dict(KRUEGER_2024_FIGURE4_BIGAUSSIAN))


__all__ += ["kim_2025_reference_iadf", "krueger_2024_figure4_iadf"]
