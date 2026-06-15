# petch public API — scope / proposal

Goal: a clean, documented public API for the open-source release. Design rules:
1. **Familiar** — mirror ViennaPS's `Domain` / `Process` shape so ViennaPS users migrate in minutes.
2. **One obvious path** — fast+accurate by default (auto-GPU), no flag incantations to get 14×.
3. **Modular escape hatches** — every internal knob still reachable for research/debugging.
4. **Portable** — same code runs GPU (fast) or CPU (slow); never crashes off-CUDA.

---

## Layer 1 — high-level (what 90% of users write)

```python
import petch

# geometry (ViennaPS MakeHole/MakeTrench analog)
dom = petch.Domain.hole(extent=14.0, dx=0.25, diameter=6.0, mask=2.0, depth=18.0)
#     petch.Domain.trench(extent=14.0, dx=0.25, width=6.0,  mask=2.0, depth=18.0)

# process model (ViennaPS SF6O2Etching analog) -- our belen model = ViennaPS physics
model = petch.SF6O2()                      # defaults match ViennaPS; override any: SF6O2(ionFlux=12, ...)

# run (auto-fast on GPU, portable on CPU)
result = petch.Process(dom, model, duration=3.0).run(steps=40)

result.depth          # center etch depth (µm)
result.max_depth      # deepest point
result.aspect_ratio   # depth / feature width
result.mesh           # (verts, faces) surface triangle mesh
result.phi            # the level-set array (nx,ny,nz)
result.wall_time      # seconds
result.save("etch.vtk")        # surface mesh to VTK (ViennaPS/ParaView compatible)
```

Compare to ViennaPS (near 1:1, so the migration story writes itself):
```python
# ViennaPS                                  # petch
d = v3.Domain()                             dom = petch.Domain.hole(...)
v3.MakeHole(domain=d, ...).apply()
m = v3.SF6O2Etching(params)                 model = petch.SF6O2(...)
p = v3.Process(); p.setDomain(d)            proc = petch.Process(dom, model, duration=3.0)
p.setProcessModel(m); p.setProcessDuration(3.0)
p.apply()                                   result = proc.run(steps=40)
```

## Layer 2 — config objects (override physics / numerics, still clean)

```python
model = petch.SF6O2(ionFlux=12.0, Fflux=1800.0, Oflux=100.0, meanE=100.0, sigmaE=10.0,
                    iedf="gauss")           # "mean" | "gauss" (=ViennaPS) | "bimodal" (beyond)
proc  = petch.Process(dom, model, duration=3.0,
                      n_ion=30000, n_neu=30000,        # ray budget
                      transport="mc",                  # "mc" | "radiosity" (Knudsen conductance)
                      device="auto")                   # "auto" | "cuda" | "cpu"
```

## Layer 3 — full modular control (research / current low-level)

The existing primitives stay public and unchanged (nothing breaks):
```python
from petch import threed as t3
geo = t3.run_etch_3d(trench_width=6.0, dx=0.25, n_steps=40, ..., par=dict(petch.PAR), flags=petch.Flags(...))
```
Layer 1/2 are thin wrappers over this. The `Flags` (each GPU speedup, chemistry model, transport) and the
`PAR` dict remain the source of truth.

## Differentiable (the unique capability)

```python
# gradient of an outcome (e.g. depth) w.r.t. a process parameter, through the simulator
grad = petch.gradient(dom, model, duration=3.0, wrt="meanE", of="depth")
# inverse design: find the recipe that hits a target profile
recipe = petch.inverse_design(dom, target_depth=10.0, wrt=["meanE", "ionFlux"])
```

---

## What to build (implementation checklist)

- [ ] `petch.Domain` with `.hole()` / `.trench()` classmethods → wraps `make_trench_3d`, holds geo dict.
- [ ] `petch.SF6O2(**overrides)` → returns a params object (wraps `PAR` + sets belen/viennaps/iedf flags).
- [ ] `petch.Process(dom, model, duration, **opts).run(steps)` → calls `run_etch_3d`, returns `Result`.
- [ ] `petch.Result` → `.depth/.max_depth/.aspect_ratio/.mesh/.phi/.wall_time/.save()`.
- [ ] `device="auto"` resolves to cuda if available (already the internal default).
- [ ] `petch.gradient()` / `petch.inverse_design()` → wrap the existing `wp.Tape` autodiff demo.
- [ ] Keep Layer 3 (`run_etch_3d`, `Flags`, `PAR`) exported and documented as the power path.
- [ ] `__init__.py` exports the Layer-1/2 names; docstrings on each.

## Open questions for review
1. `Domain.hole(depth=...)` — is `depth` the substrate depth (Lz - sub_top region) or the initial mask
   hole depth? Propose: `depth` = how deep the substrate is (etchable region); `mask` = mask thickness.
2. Units: keep µm + min (ViennaPS convention) at the API surface? Propose: yes.
3. `duration` vs `steps`: ViennaPS uses physical duration; we use steps × dt. Propose: `duration` is
   physical, `steps` is the discretization (more steps = finer/accurate), `dt = duration/steps`.
4. Rate calibration: ViennaPS's absolute rate vs our `rate_scale`. Propose: `SF6O2` ships a sane default
   `rate_scale`; expose `model.rate_scale` for matching a specific tool.
