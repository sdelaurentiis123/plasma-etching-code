"""Bosch DRIE (time-multiplexed etch) cycle driver — SEM-gated scalloping.

Emulation tier of the standard cycle (Zhou'04 / Ertl'10 / ViennaPS, but GATED — see
BOSCH_BENCHMARK_SPEC.md): per cycle (1) conformal polymer passivation of the whole profile,
(2) directional ion punch-through clears the coat on ion-visible (up-facing) surfaces only,
(3) the etch step removes Si from CLEARED surfaces: an isotropic component (F neutrals — a
Huygens/disc front of radius r_iso, which is what carves the circular-arc scallop) plus a
directional component d_dir; coated sidewalls are protected. Repeat N cycles.

Scallops, undercut, pitch and total depth all EMERGE from the cycle mechanics on the grid;
nothing is drawn. Gates (published SEM data): Ayon 1999 (via the McVittie NNIN deck) and
Tillocher 2021 — see BOSCH_BENCHMARK_SPEC.md for the exact numbers and tolerances.

Optional aspect-ratio transport attenuation: r_iso and d_dir can be scaled per cycle by the
Ertl & Selberherr analytic bottom flux F/F0 = 1 - (2x/sqrt(1+4x^2))^(k+1) at the current
aspect ratio x (their Eq. 20, <1% vs their MC) — this produces ARDE and the scallop-shrink-
with-depth (ARDSA) behavior from transport, not fitting.
"""
import numpy as np
from scipy import ndimage


def ertl_flux_factor(ar, kappa=100.0):
    """Ertl & Selberherr Eq. 20: normalized bottom-center flux vs aspect ratio (kappa = source
    exponent; large kappa = collimated ions, small = diffuse neutrals)."""
    x = max(float(ar), 0.0)
    return 1.0 - (2.0 * x / np.sqrt(1.0 + 4.0 * x * x)) ** (kappa + 1.0)


def run_bosch(width_um, n_cycles, r_iso_um, d_dir_um, dx_um=0.02, mask_h_um=0.5,
              domain_depth_um=None, margin_um=1.0, transport_kappa=None, progress=0):
    """Run N Bosch cycles on a masked trench. Returns dict with the phi grid, profile metrics
    (depth, pitch, scallop depth via the Park tangent protocol, undercut) and the wall profile.

    width_um: mask opening; r_iso_um: isotropic bite per cycle; d_dir_um: directional advance
    per cycle; transport_kappa: if set, attenuate the per-cycle etch by ertl_flux_factor(AR)
    with this exponent (neutral-limited ~ 1-3, ion-limited ~ 50-200)."""
    dx = float(dx_um)
    if domain_depth_um is None:
        # per-cycle floor advance = d_dir + r_iso (swept-disc Minkowski sum), plus headroom
        domain_depth_um = n_cycles * (d_dir_um + r_iso_um) * 1.15 + 2.0
    nx = int(round((width_um + 2 * margin_um) / dx))
    nz = int(round((domain_depth_um + mask_h_um) / dx))
    solid = np.ones((nx, nz), bool)              # True = etchable Si
    mask = np.zeros((nx, nz), bool)              # PR mask (not etchable)
    zm = int(round(mask_h_um / dx))
    x0 = int(round(margin_um / dx)); x1 = nx - x0
    mask[:, :zm] = True
    mask[x0:x1, :zm] = False                     # mask opening
    solid[:, :zm] = False                        # mask rows are not Si
    gas = np.zeros((nx, nz), bool)
    gas[:, :zm] = ~mask[:, :zm]                  # open region above the wafer + in the opening
    r_cells = r_iso_um / dx
    d_cells = d_dir_um / dx
    depth_hist = []
    for c in range(n_cycles):
        # (2) ion punch-through clears UP-FACING surfaces (gas directly above) not shadowed by the
        # MASK. Only the mask casts a hard shadow: real ions have a few degrees of divergence, more
        # than enough to clear the ~100 nm scallop-crest slivers, so crest self-shadowing of the
        # near-wall floor is NOT applied (a zero-divergence ray rule tapers the trench unphysically).
        mask_shadow = np.logical_or.accumulate(mask, axis=1)
        surf = solid & ndimage.binary_dilation(gas)         # surface Si cells
        above_gas = np.zeros_like(gas); above_gas[:, 1:] = gas[:, :-1]
        # punch-through clears only NEAR-HORIZONTAL surface (ion flux ~ cos(incidence): the flat
        # floor clears; the curved scallop feet/crests -- cells with a lateral gas neighbour -- keep
        # their polymer). This is what preserves the per-cycle arc as a scallop instead of shaving it.
        lat_solid = np.ones_like(gas)
        lat_solid[1:, :] &= ~gas[:-1, :]
        lat_solid[:-1, :] &= ~gas[1:, :]
        cleared = surf & above_gas & lat_solid & ~mask_shadow
        if not cleared.any():
            break
        # transport attenuation at the current aspect ratio
        fac = 1.0
        if transport_kappa is not None:
            depth_now = (np.argmin(gas[nx // 2, ::-1]) if gas[nx // 2].any() else 0)
            depth_now = float(np.max(np.where(gas[nx // 2])[0]) - zm + 1) * dx if gas[nx // 2, zm:].any() else 0.0
            fac = ertl_flux_factor(depth_now / width_um, transport_kappa)
        # (3) etch, SEQUENTIAL semantics (this is what the published scallop geometry demands --
        # the simultaneous/swept-disc model provably gives s = r - sqrt(r^2 - ((p-d)/2)^2) ~ 32 nm,
        # 4x smaller than Ayon's 140 nm): the ion-driven directional punch advances the floor d_dir
        # FIRST (a straight column, no lateral spread), THEN the isotropic F-neutral front expands
        # radius r_iso from the NEW floor. Disc centered at the punch bottom gives, analytically:
        # advance = d + r = 434, crest at p/2 above center -> s = r - sqrt(r^2-(p/2)^2) = 140,
        # max lateral reach r at the disc-center level -> undercut ~ r = 238. All three match Ayon.
        dsteps = int(round(d_cells * fac))
        punched = cleared.copy()
        bottom = cleared.copy()
        if dsteps > 0:
            acc = cleared.copy()
            sh = cleared.copy()
            for _ in range(dsteps):
                nxt = np.zeros_like(sh); nxt[:, 1:] = sh[:, :-1]
                acc |= nxt
                sh = nxt
            punched = acc          # the straight punched column
            bottom = sh            # the new floor after the punch
        dist = ndimage.distance_transform_edt(~bottom)      # iso front from the NEW floor
        removed = ((dist <= (r_cells * fac)) | punched) & solid & ~mask
        solid[removed] = False
        gas |= removed
        # gas connectivity: only regions connected to the top are really gas (no tunnels)
        lab, _ = ndimage.label(gas)
        top_ids = np.unique(lab[:, 0]); top_ids = top_ids[top_ids > 0]
        gas = np.isin(lab, top_ids)
        if progress and (c + 1) % progress == 0:
            d_now = (np.max(np.where(gas[nx // 2])[0]) - zm + 1) * dx if gas[nx // 2, zm:].any() else 0
            print(f"  cycle {c+1}/{n_cycles}: depth={d_now:.2f} um", flush=True)
        depth_hist.append(gas[nx // 2, zm:].sum() * dx)
    # ---- metrics ----
    cc = nx // 2
    depth = float(gas[cc, zm:].sum()) * dx
    # wall profile: leftmost gas x per row inside the trench (left sidewall)
    rows = np.arange(zm, zm + int(depth / dx) - 2)
    wall = np.array([np.argmax(gas[:, r]) for r in rows], float) * dx   # first gas cell from left
    # undercut: how far the wall goes left of the mask opening edge
    undercut = float(max(0.0, (x0 * dx) - wall.min()))
    # scallop metrics on the mid-band of the wall (Park protocol: crest line minus valley)
    band = wall[len(wall) // 6: len(wall) // 2]
    z_band = rows[len(wall) // 6: len(wall) // 2] * dx
    # crests = local minima of wall x (peaks toward Si); valleys = local maxima
    from scipy.signal import find_peaks
    crest_i, _ = find_peaks(-band, distance=max(2, int(0.5 * d_dir_um / dx)))
    valley_i, _ = find_peaks(band, distance=max(2, int(0.5 * d_dir_um / dx)))
    pitch = float(np.mean(np.diff(z_band[crest_i]))) if len(crest_i) > 2 else float("nan")
    scallop = float(np.mean(band[valley_i]) - np.mean(band[crest_i])) if (len(crest_i) > 1 and len(valley_i) > 1) else float("nan")
    return dict(depth=depth, pitch=pitch, scallop=scallop, undercut=undercut,
                wall=wall, rows=rows, solid=solid, gas=gas, mask=mask, dx=dx, zm=zm,
                depth_hist=np.array(depth_hist))
