"""petch high-level API — a clean, ViennaPS-shaped public interface.

    import petch
    dom    = petch.Domain.hole(extent=14, dx=0.25, diameter=6, mask=2, depth=18)
    model  = petch.SF6O2()
    result = petch.Process(dom, model, duration=3.0).run(steps=40)
    result.depth          # center etch depth (um)
    result.save("etch.vtk")

Mirrors ViennaPS's Domain / Process so its users migrate in minutes; underneath it calls the modular
`run_etch_3d` with the auto-fast (GPU) / portable (CPU) config. The low-level `run_etch_3d`, `Flags`,
and `PAR` stay public for full control.
"""
import time
import numpy as np

from .params import PAR, Flags
from . import threed as _t3

# friendly API name -> internal PAR key
_ALIASES = {"meanE": "Emean", "sigmaE": "Esig", "iedf": "ied_mode"}


class Domain:
    """Etch geometry: a circular hole or a line trench in a masked substrate. Units: micrometres.

    `extent` = lateral domain size (square), `dx` = grid spacing, `mask` = mask thickness,
    `depth` = etchable substrate depth (how far down the feature can go)."""

    def __init__(self, geo):
        self._geo = geo                                  # kwargs for run_etch_3d geometry

    @classmethod
    def hole(cls, extent, dx, diameter, mask=2.0, depth=16.0, headroom=4.0):
        sub_top = float(depth)
        return cls(dict(Lx=extent, Ly=extent, Lz=sub_top + mask + headroom, dx=dx,
                        trench_width=diameter, mask_th=mask, sub_top=sub_top, hole=True))

    @classmethod
    def trench(cls, extent, dx, width, mask=2.0, depth=16.0, headroom=4.0):
        sub_top = float(depth)
        return cls(dict(Lx=extent, Ly=extent, Lz=sub_top + mask + headroom, dx=dx,
                        trench_width=width, mask_th=mask, sub_top=sub_top, hole=False))

    @property
    def feature_width(self):
        return self._geo['trench_width']


class SF6O2:
    """SF6/O2 plasma-etch process model — the same physics as ViennaPS `SF6O2Etching` (Belen coupled
    F/O coverages + chemical/sputter/ion-enhanced etch, Steinbruchel sqrt(E) yields, ViennaPS angular
    yields, MC neutral transport, flux smoothing). Override any plasma parameter by keyword:

        SF6O2(ionFlux=12, Fflux=1800, Oflux=100, meanE=100, sigmaE=10, iedf="gauss")

    `iedf`: "mean" | "gauss" (=ViennaPS Gaussian IEDF) | "bimodal" (real RF-sheath, beyond ViennaPS).
    `rate_scale` calibrates absolute etch rate to a specific tool (ViennaPS `unitConversion` analog)."""

    def __init__(self, rate_scale=0.025, iedf="gauss", **overrides):
        self.par = dict(PAR)
        self.par['rate_scale'] = rate_scale
        self.par['ied_mode'] = iedf
        for k, v in overrides.items():
            self.par[_ALIASES.get(k, k)] = v
        # physics flags = the FULL faithful-ViennaPS config: belen coupled coverages + ViennaPS angular
        # yields + coverage-dependent sticking + faithful ion reflection (sticking=0, coned-cosine,
        # energy loss -> funnels ions to the deep floor; the deep-AR ARDE term, gated on ion_reflection).
        # warm_start_coverage + sobol are accuracy-neutral speedups. This is the config all the accuracy
        # validation was run with; do NOT drop ion_reflection (it ~3x's the deep-feature rate -> rate_scale
        # is calibrated for the faithful path).
        self.flags = dict(chemistry="belen", yield_angular="viennaps", coverage_sticking=True,
                          ion_reflection=True, warm_start_coverage=True, sampling="sobol")


class Result:
    """Outcome of a Process.run() — etch metrics + the surface/level-set, with VTK/OBJ export."""

    def __init__(self, geo, wall_time):
        self._geo = geo
        self.wall_time = wall_time

    @property
    def depth(self):
        """Center etch depth (um)."""
        return _t3.center_depth_3d(self._geo)

    @property
    def max_depth(self):
        return _t3.max_depth_3d(self._geo)

    @property
    def aspect_ratio(self):
        return self.max_depth / self._geo['trench_width']

    @property
    def phi(self):
        """The signed-distance level-set array (nx, ny, nz)."""
        return self._geo['phi']

    @property
    def mesh(self):
        """Surface triangle mesh as (verts, faces)."""
        v, f, _, _ = _t3.extract_mesh_3d(self._geo['phi'], self._geo['dx'])
        return v, f

    def plot(self, path=None, show=False):
        """Quick matplotlib cross-section (the etch profile through the feature centre). Returns the
        Axes; pass `path` to save a PNG. The fast way to SEE the result."""
        import matplotlib
        if path is not None and not show:
            matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        g = self._geo
        jc = g['phi'].shape[1] // 2                       # mid-y slice (x-z plane through the feature)
        phi = g['phi'][:, jc, :].T                        # (nz, nx)
        fig, ax = plt.subplots(figsize=(5, 4))
        ax.contourf(g['xs'], g['zs'], (phi < 0).astype(float), levels=[0.5, 1.5], colors=["#cfe8ff"])
        ax.contour(g['xs'], g['zs'], phi, levels=[0.0], colors="k", linewidths=1.2)
        ax.axhline(g['sub_top'], color="0.6", lw=0.6, ls="--")
        ax.set_xlabel("x (µm)"); ax.set_ylabel("z (µm)")
        ax.set_title(f"depth {self.depth:.1f} µm  AR {self.aspect_ratio:.1f}  ({self.wall_time:.1f}s)")
        ax.set_aspect("equal")
        if path is not None:
            fig.tight_layout(); fig.savefig(path, dpi=130)
        if show:
            plt.show()
        return ax

    def save(self, path):
        """Write the surface mesh to `path` (.vtk for ParaView/ViennaPS, or .obj)."""
        v, f = self.mesh
        if str(path).lower().endswith(".obj"):
            with open(path, "w") as fh:
                for p in v:
                    fh.write(f"v {p[0]} {p[1]} {p[2]}\n")
                for t in f:
                    fh.write(f"f {t[0]+1} {t[1]+1} {t[2]+1}\n")
        else:                                            # legacy VTK PolyData (ParaView/ViennaPS read it)
            with open(path, "w") as fh:
                fh.write("# vtk DataFile Version 3.0\npetch surface\nASCII\nDATASET POLYDATA\n")
                fh.write(f"POINTS {len(v)} float\n")
                for p in v:
                    fh.write(f"{p[0]} {p[1]} {p[2]}\n")
                fh.write(f"POLYGONS {len(f)} {4*len(f)}\n")
                for t in f:
                    fh.write(f"3 {int(t[0])} {int(t[1])} {int(t[2])}\n")
        return path


class Process:
    """Run a process `model` on a `domain` for `duration` (minutes). On a GPU the full fast pipeline
    auto-enables (~14x faster than ViennaPS); on CPU it falls back to the portable path. `steps` sets
    the time discretization (dt = duration/steps; more steps = finer)."""

    def __init__(self, domain, model, duration, n_ion=30000, n_neu=30000, transport="mc"):
        self.domain = domain
        self.model = model
        self.duration = float(duration)
        self.n_ion = n_ion
        self.n_neu = n_neu
        self.transport = transport

    def run(self, steps=40, verbose=False):
        flags = Flags(neutral_transport=self.transport, **self.model.flags)
        t0 = time.time()
        geo = _t3.run_etch_3d(t_end=self.duration, n_steps=steps,
                              par=dict(self.model.par), flags=flags,
                              n_ion=self.n_ion, n_neu=self.n_neu, reinit_method="fsm",
                              verbose=verbose, **self.domain._geo)
        return Result(geo, time.time() - t0)
