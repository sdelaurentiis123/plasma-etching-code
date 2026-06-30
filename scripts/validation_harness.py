"""Shared validation scorecard: petch (MC / Knudsen / radiosity) AND Craig Xu Chen's plasma_sim
(DDA / Knudsen) scored on the SAME experimental axes — the de Boer/Blauw trench ARDE (RMSE vs the
measured wafer) and the Gomez/Belen open-field rate — using the fast static-geometry method from
Craig's bench/industrial_validation.py: carve a feature at a fixed (width, depth), compute the
velocity field ONCE, take the mean floor rate. No etch evolution.

Both engines run with MATCHED Belen SF6/O2 chemistry (only the transport method differs). The
normalized ARDE RMSE is rate_scale-independent (it compares the rolloff SHAPE to the wafer); the
open-field rate is reported as the rate_scale needed to hit 1.3 um/min.

Gates (experimental_validation.json): de Boer trench ARDE RMSE <= 0.05; open-rate rel. error <= 0.10.
Expected honest result: both ballistic engines agree with EACH OTHER but sit above the real wafer at
high AR (the structural gas-conductance/charging gap), so default (ViennaPS-regime) params miss the
wafer gate -- which is the point of scoring them side by side.

petch runs on CPU (Warp-CPU); plasma_sim runs on the Apple Metal GPU (mlx).
Writes viz/validation_scorecard.png + validation_scorecard.npz and prints the scorecard table.
"""
import os, sys, time
import numpy as np
import skfmm
os.environ.setdefault("PETCH_DEVICE", "cpu")
sys.path.insert(0, "/Users/stanislavdelaurentiis/chip-etch")          # plasma_sim
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

import petch
from petch import threed as t3
import plasma_sim as ps
from plasma_sim.solver3d import Config3D, Solver3D

# ---- experimental targets + gates (bench/data/experimental_validation.json) ----
DEBOER_AR = np.array([0.0, 10.0, 20.0, 40.0])
DEBOER_NR = np.array([1.0, 0.43, 0.29, 0.20])
GOMEZ_OPEN_UM_MIN = 1.3
RMSE_GATE = 0.05
OPEN_ERR_GATE = 0.10

W = 2.0            # de Boer trench width (um)
DX = 0.25
MASK = 0.5
ARS = [0.0, 10.0, 20.0, 40.0]

# ---- matched Belen SF6/O2 chemistry ----
def petch_par():
    p = dict(petch.PAR)
    p.update(ied_mode='gauss', rate_scale=1.0, radiosity_solver='gmres',
             ionFlux=12.0, Fflux=1800.0, Oflux=100.0, k_sigma=300.0, beta_sigma=0.04, B_sp=9.3,
             betaE=0.7, betaO=1.0, A_ie=7.0, A_sp=0.0337, A_p=3.0, Eth_ie=15.0, Eth_sp=20.0,
             Eth_p=10.0, rho=5.02, Emean=100.0, Esig=10.0)
    return p

def craig_model():
    return ps.SF6O2.default(ion_flux=12, f_flux=1800, o_flux=100, cal_f=1.0, ied_mode='gauss',
                            mean_energy=100, sigma_energy=10, k_sigma=300, beta_sigma=0.04, b_sp=9.3,
                            beta_f=0.7, beta_o=1.0, a_ie=7, a_sp=0.0337, a_p=3,
                            eth_ie=15, eth_sp=20, eth_p=10, rho=5.02, rate_scale=1.0)

# ---- petch static-geometry floor rate (um/s) at a given aspect ratio ----
def petch_static_rate(ar, shape, transport, ion_refl):
    depth = ar * W
    sub_top = 6.0                                  # substrate-top z; feature etched into z < sub_top
    Lz = sub_top + MASK + 0.5
    z_min_feature = sub_top - depth
    Lz = max(Lz, sub_top + MASK + 0.5)
    sub = max(depth + 6 * DX, sub_top)             # ensure substrate deep enough below the floor
    sub_top = sub
    Lz = sub_top + MASK + 0.5
    Lx = Ly = max(4.0 * W, 6.0)
    nx = int(round(Lx / DX)); ny = int(round(Ly / DX)); nz = int(round(Lz / DX))
    xs = np.arange(nx) * DX; ys = np.arange(ny) * DX; zs = np.arange(nz) * DX
    X, Y, Z = np.meshgrid(xs, ys, zs, indexing='ij')
    cx, cy = xs.mean(), ys.mean()
    if shape == "hole":
        inside = np.sqrt((X - cx) ** 2 + (Y - cy) ** 2) < W / 2.0
    else:
        inside = np.abs(X - cx) < W / 2.0
    feature_void = inside & (Z > sub_top - depth) & (Z < sub_top + MASK + DX)   # carved column (+ open mouth)
    mask_slab = (Z >= sub_top) & (Z < sub_top + MASK) & (~inside)
    solid = ((Z < sub_top) & (~feature_void)) | mask_slab
    phi = skfmm.distance(np.where(solid, 1.0, -1.0).astype(np.float64), dx=DX)   # petch: phi>0 solid
    geo = dict(xs=xs, ys=ys, zs=zs, dx=DX, phi=phi, mask=mask_slab,
               Lx=Lx, Ly=Ly, Lz=Lz, sub_top=sub_top, trench_width=W, hole=(shape == "hole"))
    par = petch_par()
    fl = petch.Flags(chemistry="belen", yield_angular="viennaps", coverage_sticking=True,
                     ion_reflection=ion_refl, sampling="sobol", neutral_transport=transport)
    verts, faces, centroids, areas = t3.extract_mesh_3d(phi, DX)
    centroids = centroids + np.array([0, 0, 0.0])
    mesh = t3.wp.Mesh(points=t3.wp.array(verts.astype(np.float32), dtype=t3.wp.vec3, device=t3.DEVICE),
                      indices=t3.wp.array(faces.flatten().astype(np.int32), dtype=int, device=t3.DEVICE))
    _nt = transport
    if _nt == "radiosity":
        m_i, m_F, m_O, cos_i = t3.mc_flux_3d_radiosity(mesh, verts, faces, centroids, areas, geo, par,
                                                       n_ion=8000, seed=1, flags=fl)
    elif _nt == "knudsen":
        m_i, m_F, m_O, cos_i = t3.mc_flux_3d_knudsen(mesh, verts, faces, centroids, areas, geo, par,
                                                     n_ion=8000, seed=1, flags=fl)
    else:
        m_i, m_F, m_O, cos_i, _ = t3.mc_flux_3d_coupled(mesh, verts, faces, areas, geo, par,
                                                        n_ion=8000, n_neu=8000, seed=1,
                                                        sampling="sobol", flags=fl)
    is_mask = t3.faces_in_mask(centroids, geo, MASK, W, hole=(shape == "hole"))
    V = t3.surface_rate(m_i, m_F, m_O, cos_i, is_mask, par, flags=fl)
    V = np.nan_to_num(V, nan=0.0, posinf=0.0, neginf=0.0)
    # floor faces: near the carved floor plane, near center, facing up, not mask
    cz, crr = centroids[:, 2], np.sqrt((centroids[:, 0] - cx) ** 2 + (centroids[:, 1] - cy) ** 2)
    if shape == "trench":
        crr = np.abs(centroids[:, 0] - cx)
    floor = (np.abs(cz - (sub_top - depth)) < 1.5 * DX) & (crr < 0.30 * W) & (~is_mask)
    if depth < DX:                                   # AR~0: open field = the whole top surface
        floor = (np.abs(cz - sub_top) < 1.5 * DX) & (~is_mask)
    return float(np.mean(V[floor])) if floor.any() else 0.0

# ---- Craig static-geometry floor rate (um/s) ----
def craig_static_rate(ar, shape, solver_name):
    depth = ar * W
    half = max(2.0 * W, 3.0)
    cfg = Config3D(x_half=half, y_half=half, z_min=-max(depth + 2 * DX, 4 * DX), z_max=1.0, dx=DX,
                   via_radius=W / 2.0, feature_shape=("via" if shape == "hole" else "trench"),
                   mask_height=MASK, sf6o2=craig_model().params, coverage_sticking=True,
                   coverage_iters=3, neutral_solver=solver_name, use_gpu=True)
    solver = Solver3D(cfg)
    if shape == "trench":
        inside = np.abs(solver.X) < W / 2.0
    else:
        inside = np.sqrt(solver.X ** 2 + solver.Y ** 2) < W / 2.0
    gas = (solver.Z > 0.0) | (inside & (solver.Z > -depth))
    solver.phi = skfmm.distance(np.where(gas, 1.0, -1.0), dx=DX)   # plasma_sim: phi>0 gas
    V, _ = solver.velocity_field()
    band = np.abs(solver.phi) < 1.25 * DX
    if shape == "trench":
        cen = (np.abs(solver.X) < 0.30 * W) & (np.abs(solver.Y) < max(0.75 * solver.cfg.y_half, DX))
    else:
        cen = np.sqrt(solver.X ** 2 + solver.Y ** 2) < max(0.30 * W, DX)
    zsel = (np.abs(solver.Z + depth) < 1.25 * DX) if depth >= DX else (np.abs(solver.Z) < 1.25 * DX)
    floor = band & cen & zsel & (~solver.mask)
    return float(np.mean(V[floor])) if floor.any() else 0.0

def score(rates):
    rates = np.asarray(rates, float)
    nr = rates / max(rates[0], 1e-30)
    at_db = np.interp(DEBOER_AR, ARS, nr)
    rmse = float(np.sqrt(np.mean((at_db - DEBOER_NR) ** 2)))
    open_um_min = 60.0 * rates[0]
    rate_scale_for_gomez = GOMEZ_OPEN_UM_MIN / max(open_um_min, 1e-30)
    return nr, rmse, open_um_min, rate_scale_for_gomez

CASES = [
    ("petch mc",        "petch", "mc",        False),
    ("petch knudsen",   "petch", "knudsen",   False),
    ("petch radiosity", "petch", "radiosity", False),
    ("petch mc+ionrefl","petch", "mc",        True),
    ("Craig DDA",       "craig", "gmres",     None),
    ("Craig knudsen",   "craig", "knudsen",   None),
]

def main():
    for shape in ("trench", "hole"):
        print(f"\n================ {shape}  (W={W} um, dx={DX}, ARs={ARS}) ================", flush=True)
        print(f"{'config':18s} {'nr@AR':>28s}   {'RMSE':>6s} {'gate':>5s}  {'open um/min':>11s}", flush=True)
        print(f"{'de Boer wafer':18s} {str(list(DEBOER_NR)):>28s}   {'--':>6s} {'--':>5s}", flush=True)
        rows = {}
        for label, engine, mode, ion_refl in CASES:
            try:
                t0 = time.time()
                rates = [(petch_static_rate(ar, shape, mode, ion_refl) if engine == "petch"
                          else craig_static_rate(ar, shape, mode)) for ar in ARS]
                nr, rmse, openr, rsg = score(rates)
                wall = time.time() - t0
                rows[label] = dict(rates=rates, nr=nr.tolist(), rmse=rmse, open_um_min=openr,
                                   rate_scale_for_gomez=rsg, wall=wall)
                gate = "PASS" if rmse <= RMSE_GATE else "fail"
                print(f"{label:18s} {np.array2string(nr,precision=2,floatmode='fixed'):>28s}   "
                      f"{rmse:6.3f} {gate:>5s}  {openr:11.3g}   ({wall:.0f}s)", flush=True)
            except Exception as e:
                print(f"{label:18s} FAILED: {type(e).__name__}: {e}", flush=True)
        np.savez(os.path.join(HERE, f"validation_scorecard_{shape}.npz"),
                 ars=ARS, deboer_ar=DEBOER_AR, deboer_nr=DEBOER_NR, rows=rows, allow_pickle=True)

if __name__ == "__main__":
    main()
