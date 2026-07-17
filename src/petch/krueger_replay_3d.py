"""Source-conditioned electrostatics for the Krüger 2024 HAR replay.

The published MCFPM configuration uses periodic/reflective lateral electrostatic boundaries, a
zero-normal-gradient top, and a grounded bottom. This module maps that boundary-value problem to
the common Q1 Poisson operator. The amorphous-carbon relative permittivity is intentionally a
required input because Krüger et al. do not publish it for the experimental hard mask.
"""
from __future__ import annotations

import numpy as np

from .charging_poisson_3d import NodalPoissonSystem3D
from .feature_step_3d import FeatureGeometry3D


KRUEGER_SIO2_MATERIAL_ID = 1
KRUEGER_AMORPHOUS_CARBON_MATERIAL_ID = 2


def _cell_center_average(field):
    value = np.asarray(field, dtype=float)
    return sum(
        value[i:i + value.shape[0] - 1,
              j:j + value.shape[1] - 1,
              k:k + value.shape[2] - 1]
        for i in (0, 1) for j in (0, 1) for k in (0, 1)) / 8.0


def make_krueger_2024_poisson_system_3d(
        geometry: FeatureGeometry3D, *, mask_relative_permittivity: float,
        sio2_relative_permittivity: float = 3.9):
    """Build the source-declared MCFPM electrostatic boundary-value problem.

    The bottom node plane is the grounded reservoir. Leaving the top node plane unconstrained is
    the natural zero-normal-gradient boundary of :class:`NodalPoissonSystem3D`; both lateral axes
    are periodic. No plasma volume charge or local analytic electron-current law is added.
    """
    if not isinstance(geometry, FeatureGeometry3D):
        raise TypeError("geometry must be a FeatureGeometry3D")
    values = np.asarray(
        [mask_relative_permittivity, sio2_relative_permittivity], dtype=float)
    if np.any(~np.isfinite(values)) or np.any(values <= 0.0):
        raise ValueError("material relative permittivities must be positive and finite")
    if (
        geometry.material_levelsets is None
        or set(geometry.material_levelsets) != {
            KRUEGER_SIO2_MATERIAL_ID,
            KRUEGER_AMORPHOUS_CARBON_MATERIAL_ID,
        }
    ):
        raise ValueError("geometry is not a Krueger SiO2/amorphous-carbon replay cell")

    material_ids = np.asarray(sorted(geometry.material_levelsets), dtype=int)
    center_fields = np.stack([
        _cell_center_average(geometry.material_levelsets[int(material_id)])
        for material_id in material_ids])
    owner = material_ids[np.argmax(center_fields, axis=0)]
    solid = np.max(center_fields, axis=0) >= 0.0
    epsilon_r = np.ones(owner.shape, dtype=float)
    epsilon_r[solid & (owner == KRUEGER_SIO2_MATERIAL_ID)] = float(
        sio2_relative_permittivity)
    epsilon_r[
        solid & (owner == KRUEGER_AMORPHOUS_CARBON_MATERIAL_ID)
    ] = float(mask_relative_permittivity)

    fixed = np.zeros(geometry.phi.shape, dtype=bool)
    fixed[:, :, 0] = True
    return NodalPoissonSystem3D(
        epsilon_r,
        geometry.dx * geometry.mesh_length_unit_m,
        fixed,
        periodic_axes=(0, 1),
    )
