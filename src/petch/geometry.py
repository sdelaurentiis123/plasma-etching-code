"""Geometry / level-set surface extraction. Ported verbatim from feature_etch.py."""
import numpy as np
import skfmm
from skimage import measure


def make_trench(W, H, dx, trench_width, mask_thickness, sub_top):
    """Build phi (signed distance, >0 solid) and a static mask occupancy field."""
    nx, ny = int(round(W / dx)), int(round(H / dx))
    xs = (np.arange(nx) + 0.5) * dx
    ys = (np.arange(ny) + 0.5) * dx
    X, Y = np.meshgrid(xs, ys, indexing='ij')
    solid = Y < sub_top                                   # substrate fills below sub_top
    mask_band = (Y >= sub_top) & (Y < sub_top + mask_thickness)
    opening = np.abs(X - W / 2) < trench_width / 2
    mask = mask_band & (~opening)                         # mask everywhere except the opening
    solid = solid | mask
    phi0 = np.where(solid, 1.0, -1.0)
    phi = skfmm.distance(phi0, dx=dx)                     # signed distance, >0 solid
    return X, Y, xs, ys, phi, mask, nx, ny


def extract_surface(phi, xs, ys, dx):
    """Marching-squares contour of phi=0 -> ordered polyline segments + outward normals."""
    contours = measure.find_contours(phi, 0.0)
    segs = []
    for c in contours:
        px = xs[0] + c[:, 0] * dx
        py = ys[0] + c[:, 1] * dx
        pts = np.column_stack([px, py])
        for k in range(len(pts) - 1):
            segs.append((pts[k, 0], pts[k, 1], pts[k + 1, 0], pts[k + 1, 1]))
    segs = np.array(segs)  # (M,4): x0,y0,x1,y1
    if len(segs) == 0:
        return segs, segs, segs, segs
    mid = np.column_stack([(segs[:, 0] + segs[:, 2]) / 2, (segs[:, 1] + segs[:, 3]) / 2])
    tang = np.column_stack([segs[:, 2] - segs[:, 0], segs[:, 3] - segs[:, 1]])
    L = np.hypot(tang[:, 0], tang[:, 1]); L[L == 0] = 1e-12
    nrm = np.column_stack([-tang[:, 1], tang[:, 0]]) / L[:, None]   # one of the two normals
    return segs, mid, nrm, L


def orient_normals(mid, nrm, phi, xs, ys, dx):
    """Flip normals so they point into gas (phi<0)."""
    gx, gy = np.gradient(phi, dx)
    ix = np.clip(((mid[:, 0] - xs[0]) / dx).astype(int), 0, phi.shape[0] - 1)
    iy = np.clip(((mid[:, 1] - ys[0]) / dx).astype(int), 0, phi.shape[1] - 1)
    g = np.column_stack([gx[ix, iy], gy[ix, iy]])
    gl = np.hypot(g[:, 0], g[:, 1]); gl[gl == 0] = 1e-12
    outward = -g / gl[:, None]            # outward (into gas) = -grad phi
    dot = (nrm * outward).sum(1)
    nrm = np.where(dot[:, None] < 0, -nrm, nrm)
    return nrm


def seg_in_mask(mid, mask, xs, ys, dx):
    ix = np.clip(((mid[:, 0] - xs[0]) / dx).astype(int), 0, mask.shape[0] - 1)
    iy = np.clip(((mid[:, 1] - ys[0]) / dx).astype(int), 0, mask.shape[1] - 1)
    return mask[ix, iy]


def profile_bottom(phi, xs, ys, dx, W):
    """Lowest y of the phi=0 contour near the trench centre."""
    segs = extract_surface(phi, xs, ys, dx)[0]
    if len(segs) == 0:
        return ys[-1]
    cx = np.r_[segs[:, 0], segs[:, 2]]; cy = np.r_[segs[:, 1], segs[:, 3]]
    centre = np.abs(cx - W / 2) < W * 0.15
    return cy[centre].min() if centre.any() else cy.min()
