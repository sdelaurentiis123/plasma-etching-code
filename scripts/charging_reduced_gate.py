#!/usr/bin/env python3
"""No-write reduced Hwang-Giapis charging/notching gate.

This is intentionally separate from charging_gate.py and notching_gate.py because those scripts write
the canonical .npz files used by the docs figures. This runner is for fast local diagnosis: same HG
digitized targets, explicit config printout, and no filesystem output unless this file is edited.
"""
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from petch.charging import floor_balance
from petch.charging2d import solve_edge_array_charging, solve_trench_charging


HG_AR = np.array([1.0, 1.2, 1.6, 2.0, 2.6, 3.0, 3.6, 4.0])
HG_FLUX = np.array([0.59, 0.55, 0.47, 0.40, 0.34, 0.30, 0.26, 0.22])
HG_FOOT_E = np.array([15.0, 16.5, 17.5, 20.0, 23.0, 25.0, 26.5, 27.5])
HG_VPOLY = np.array([6.0, 9.0, 15.0, 20.0, 27.0, 31.0, 36.0, 39.0])


def getenv_float(name, default):
    return float(os.environ.get(name, str(default)))


def getenv_int(name, default):
    return int(os.environ.get(name, str(default)))


def main():
    cfg = {
        "see_model": os.environ.get("PETCH_SEE_MODEL", "none"),
        "see_generations": getenv_int("PETCH_SEE_GENERATIONS", 1),
        "source_model": os.environ.get("PETCH_SOURCE_MODEL", "analytic"),
        "geometry": os.environ.get("PETCH_CHARGING_GEOMETRY", "edge_array"),
        "poly_mode": os.environ.get("PETCH_POLY_MODE", "tied"),
        "poly_bias_V": getenv_float("PETCH_POLY_BIAS_V", 0.0),
        "edge_open_model": os.environ.get("PETCH_EDGE_OPEN_MODEL", "line_of_sight"),
        "edge_open_electron_flux": os.environ.get("PETCH_EDGE_OPEN_ELECTRON_FLUX"),
        "W": getenv_int("PETCH_CHARGING_W", 16),
        "mouth": getenv_int("PETCH_CHARGING_MOUTH", 80),
        "n_per_iter": getenv_int("PETCH_N_PER_ITER", 1200),
        "n_iter": getenv_int("PETCH_N_ITER", 100),
    }
    cfg["edge_open_electron_flux"] = (
        None if cfg["edge_open_electron_flux"] is None else float(cfg["edge_open_electron_flux"])
    )

    print("=== REDUCED HG CHARGING + NOTCHING GATE (NO WRITE) ===", flush=True)
    print("config: " + " ".join(f"{k}={v}" for k, v in cfg.items()), flush=True)

    def run_solver(ar, seed, **kw):
        common = dict(
            W=cfg["W"],
            mouth=cfg["mouth"],
            n_per_iter=cfg["n_per_iter"],
            n_iter=cfg["n_iter"],
            seed=seed,
            see_model=cfg["see_model"],
            see_generations=cfg["see_generations"],
            source_model=cfg["source_model"],
            edge_open_model=cfg["edge_open_model"],
            edge_open_electron_flux=cfg["edge_open_electron_flux"],
        )
        common.update(kw)
        if cfg["geometry"] == "edge_array":
            return solve_edge_array_charging(ar, **common)
        if cfg["geometry"] == "trench":
            return solve_trench_charging(
                ar,
                poly_mode=cfg["poly_mode"],
                poly_bias_V=cfg["poly_bias_V"],
                **common,
            )
        raise ValueError(f"unknown PETCH_CHARGING_GEOMETRY={cfg['geometry']}")

    floor_flux = []
    vcenter = []
    foot_e = []
    foot_flux = []
    vpoly = []
    isurv = []
    esurv = []
    residual = []

    for i, ar in enumerate(HG_AR):
        t0 = time.time()
        r = run_solver(float(ar), seed=100 + i)
        if cfg["geometry"] == "edge_array" or cfg["poly_mode"] == "edge_open":
            e = r["foot_ion_Emean_edge"]
            f = r["foot_ion_flux_edge"]
            vp = r["V_poly_neighbor"]
        else:
            e = r["foot_ion_Emean"]
            f = r["foot_ion_flux"]
            vp = r["V_poly"]
        ti = r["diag"]["trace"]["last_ion"]
        te = r["diag"]["trace"]["last_electron"]
        res = r["diag"].get("residual", {})
        res_max = max(abs(float(v)) for v in res.values()) if res else np.nan
        floor_flux.append(r["floor_flux"])
        vcenter.append(r["V_floor_center"])
        foot_e.append(e)
        foot_flux.append(f)
        vpoly.append(vp)
        isurv.append(ti["survivor_frac"] if ti else np.nan)
        esurv.append(te["survivor_frac"] if te else np.nan)
        residual.append(res_max)
        print(
            f"AR {ar:3.1f}: flux={r['floor_flux']:.3f} HG={HG_FLUX[i]:.2f} "
            f"Vc={r['V_floor_center']:5.1f} footE={e:5.1f} HG_E={HG_FOOT_E[i]:.1f} "
            f"footFlux={f:.3f} Vpoly={vp:5.1f} HG_Vp={HG_VPOLY[i]:.0f} "
            f"surv_i/e={isurv[-1]:.4f}/{esurv[-1]:.4f} res={res_max:.3f} "
            f"dt={time.time() - t0:.1f}s",
            flush=True,
        )

    floor_flux = np.array(floor_flux)
    foot_e = np.array(foot_e)
    foot_flux = np.array(foot_flux)
    vpoly = np.array(vpoly)
    isurv = np.array(isurv)
    esurv = np.array(esurv)
    residual = np.array(residual)

    rmse = float(np.sqrt(np.mean((floor_flux - HG_FLUX) ** 2)))
    rel_e = np.abs(foot_e - HG_FOOT_E) / HG_FOOT_E
    rel_v = np.abs(vpoly - HG_VPOLY) / HG_VPOLY
    deep = HG_AR >= 1.6
    flux_ratio = float(foot_flux[deep].max() / max(foot_flux[deep].min(), 1e-12))
    survivor_max = float(max(np.nanmax(isurv), np.nanmax(esurv)))
    residual_max = float(np.nanmax(residual))

    print("=== SUMMARY ===", flush=True)
    print(f"floor RMSE={rmse:.3f} gate<=0.050 {'PASS' if rmse <= 0.05 else 'fail'}", flush=True)
    print(
        f"survivor max={survivor_max:.4f} gate<0.001 "
        f"{'PASS' if survivor_max < 0.001 else 'fail'}",
        flush=True,
    )
    print(
        f"current residual max={residual_max:.3f} gate<0.080 "
        f"{'PASS' if residual_max < 0.08 else 'fail'}",
        flush=True,
    )
    print(
        f"foot energy max rel err={rel_e.max():.3f}, rising={bool(foot_e[-1] > foot_e[0])} "
        f"gate<=0.30+rising {'PASS' if (rel_e <= 0.30).all() and foot_e[-1] > foot_e[0] else 'fail'}",
        flush=True,
    )
    print(
        f"foot flux ratio AR>=1.6={flux_ratio:.2f} gate<=2.0 "
        f"{'PASS' if flux_ratio <= 2.0 else 'fail'}",
        flush=True,
    )
    print(
        f"poly potential max rel err={rel_v.max():.3f}, rising={bool(vpoly[-1] > vpoly[0])} "
        f"gate<=0.30+rising {'PASS' if (rel_v <= 0.30).all() and vpoly[-1] > vpoly[0] else 'fail'}",
        flush=True,
    )

    r300 = run_solver(4.0, seed=199, V_dc=300.0, V_rf=30.0)
    print(
        f"Matsui AR4 @300eV floor flux={r300['floor_flux']:.3f} gate>0.1 "
        f"{'PASS' if r300['floor_flux'] > 0.1 else 'fail'}",
        flush=True,
    )
    f0 = [floor_balance(a, n=120000)[1] for a in (1.0, 4.0)]
    m0 = floor_balance(4.0, V_dc=300.0, n=120000)[1]
    print(
        f"0D closure AR1={f0[0]:.3f} AR4={f0[1]:.3f} Matsui={m0:.3f} "
        f"{'PASS' if f0[0] > f0[1] and m0 > 0.05 else 'fail'}",
        flush=True,
    )


if __name__ == "__main__":
    main()
