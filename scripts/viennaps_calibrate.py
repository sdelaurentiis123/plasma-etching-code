#!/usr/bin/env python3
"""Run on a Linux+NVIDIA box (ViennaPS installs there; no arm64 wheel on the Mac).

Generates the EXACT ViennaPS SF6O2 reference at the benchmark conditions so we can:
  (1) replace the cached scalar reference with real depth + ARDE + full surface meshes
      (unlocks Chamfer / sidewall), and
  (2) calibrate our mc_flux normalization (the Phase-0 root cause) against ground truth.

Also dumps the installed SF6O2 default parameters (verify the Belen constants) and, if the
API exposes it, probes the open-field flux normalization.

Outputs: viennaps_reference.json (scalars + params) and viennaps_surfaces.npz (x,y nodes/width).
"""
import json
import numpy as np

DX, W, H, MASK, DUR = 0.25, 20.0, 24.0, 2.0, 3.0
WIDTHS = [4.0, 6.0, 8.0, 12.0]

import viennaps as ps
try:
    import viennaps.d2 as vps
except Exception:
    vps = ps  # older/newer layout fallback

print("ViennaPS version:", getattr(ps, "__version__", "?"))
print("gpuAvailable:", ps.gpuAvailable() if hasattr(ps, "gpuAvailable") else "n/a")

ps.Logger.setLogLevel(ps.LogLevel.ERROR)
ps.Length.setUnit("micrometer")
ps.Time.setUnit("min")
Mat = ps.Material


def dump_params():
    try:
        p = vps.SF6O2Etching.defaultParameters()
        out = {}
        for attr in dir(p):
            if attr.startswith("_"):
                continue
            try:
                v = getattr(p, attr)
                if isinstance(v, (int, float, bool, str)):
                    out[attr] = v
            except Exception:
                pass
        # nested structs (Ions / Substrate / Passivation / Mask) if present
        for grp in ("Ions", "Substrate", "Passivation", "Mask"):
            g = getattr(p, grp, None)
            if g is not None:
                for attr in dir(g):
                    if attr.startswith("_"):
                        continue
                    try:
                        v = getattr(g, attr)
                        if isinstance(v, (int, float, bool, str)):
                            out[f"{grp}.{attr}"] = v
                    except Exception:
                        pass
        return out
    except Exception as e:
        return {"error": str(e)}


def vps_etch(width):
    d = vps.Domain()
    vps.MakeTrench(domain=d, gridDelta=DX, xExtent=W, yExtent=H, trenchWidth=width,
                   trenchDepth=MASK, taperingAngle=0.0, baseHeight=0.0,
                   periodicBoundary=False, makeMask=True, material=Mat.Si).apply()
    model = vps.SF6O2Etching(vps.SF6O2Etching.defaultParameters())
    # ViennaPS 4.x auto-selects the GPU/OptiX engine when a GPU is present; force the CPU
    # (Embree) disk engine — identical physics, and avoids the OptiX driver-version mismatch.
    p = vps.Process()
    p.setDomain(d)
    p.setProcessModel(model)
    p.setProcessDuration(DUR)
    try:
        p.setFluxEngineType(ps.FluxEngineType.CPU_DISK)
    except Exception as e:
        print("  (setFluxEngineType failed, using default:", e, ")")
    p.apply()
    n = np.array(d.getSurfaceMesh().getNodes())
    return n[:, 0], n[:, 1]


def depth_centre(x, y, half=1.5):
    sel = y < -0.02                       # substrate top at y=0; etched region is y<0
    xc, dep = x[sel], -y[sel]
    c = np.abs(xc) < half
    return float(dep[c].max()) if c.any() else 0.0


def main():
    params = dump_params()
    print("\nSF6O2 default parameters (installed version):")
    for k in sorted(params):
        print(f"  {k} = {params[k]}")

    depths = {}
    surfaces = {}
    for w in WIDTHS:
        x, y = vps_etch(w)
        depths[str(w)] = depth_centre(x, y)
        surfaces[f"x_{w}"] = x
        surfaces[f"y_{w}"] = y
        print(f"  ViennaPS width {w}: depth = {depths[str(w)]:.3f} um  (nodes {len(x)})")

    ref = dict(widths=WIDTHS, depth=depths,
               depth8=depths["8.0"],
               arde=[depths[str(w)] for w in WIDTHS],
               params=params, DX=DX, W=W, H=H, MASK=MASK, DUR=DUR,
               version=str(getattr(ps, "__version__", "?")))
    with open("viennaps_reference.json", "w") as f:
        json.dump(ref, f, indent=2)
    np.savez("viennaps_surfaces.npz", **surfaces)
    print("\nwrote viennaps_reference.json + viennaps_surfaces.npz")


if __name__ == "__main__":
    main()
