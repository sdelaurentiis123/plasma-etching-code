"""Level-set advection, velocity extension, reinitialization. Ported verbatim from feature_etch.py.

`advect` is first-order Godunov upwind (PoC; contributor #5 swaps in WENO + TVD-RK later).
`reinit` restores the signed-distance property via fast marching (skfmm).
"""
import numpy as np
import skfmm
from scipy.spatial import cKDTree


def advect(phi, F, dx, dt):
    """phi_t + F|grad phi| = 0, F>=0 (etch shrinks solid). First-order upwind Godunov."""
    dxm = np.zeros_like(phi); dxp = np.zeros_like(phi)
    dym = np.zeros_like(phi); dyp = np.zeros_like(phi)
    dxm[1:, :] = (phi[1:, :] - phi[:-1, :]) / dx
    dxp[:-1, :] = (phi[1:, :] - phi[:-1, :]) / dx
    dym[:, 1:] = (phi[:, 1:] - phi[:, :-1]) / dx
    dyp[:, :-1] = (phi[:, 1:] - phi[:, :-1]) / dx
    grad = np.sqrt(np.maximum(dxm, 0)**2 + np.minimum(dxp, 0)**2 +
                   np.maximum(dym, 0)**2 + np.minimum(dyp, 0)**2)
    return phi - dt * F * grad


def extend_velocity(V, mid, phi, xs, ys, dx, band):
    """Extend surface velocities to grid narrow-band via nearest-segment lookup."""
    nx, ny = phi.shape
    F = np.zeros_like(phi)
    bandmask = np.abs(phi) < band
    ii, jj = np.where(bandmask)
    gx = xs[0] + ii * dx; gy = ys[0] + jj * dx
    if len(mid) == 0:
        return F
    tree = cKDTree(mid)
    _, idx = tree.query(np.column_stack([gx, gy]))
    F[ii, jj] = V[idx]
    return F


def reinit(phi, dx):
    """Reinitialize phi to a signed-distance function (fast marching)."""
    return skfmm.distance(phi, dx=dx)
