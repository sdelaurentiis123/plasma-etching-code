#!/usr/bin/env python3
"""Minimal petch example: etch a 3D hole and print the profile.

On an NVIDIA GPU this runs the full GPU pipeline (~14x faster than ViennaPS) automatically.
On CPU / Apple Silicon the SAME code runs the portable numpy/skimage path (slower).

    PETCH_DEVICE=cuda python examples/etch_hole.py     # GPU
    PETCH_DEVICE=cpu  python examples/etch_hole.py     # CPU (or just omit; defaults to cpu)
"""
import petch

# High-level, ViennaPS-shaped API. SF6O2() carries the full faithful config (Belen coupled coverages,
# ViennaPS angular yields, coverage-dependent sticking, faithful ion reflection). rate_scale sets the
# absolute etch rate (the per-tool calibration knob); it does not change the ARDE shape.
dom    = petch.Domain.hole(extent=14, dx=0.25, diameter=6, mask=2, depth=18)
model  = petch.SF6O2()
result = petch.Process(dom, model, duration=3.0).run(steps=40)

print(f"depth        : {result.depth:.2f} µm")
print(f"max depth    : {result.max_depth:.2f} µm")
print(f"aspect ratio : {result.aspect_ratio:.2f}")
print(f"wall clock   : {result.wall_time:.2f} s   (ViennaPS-GPU on the same hole: ~22 s)")

result.save("etch_hole.vtk")          # ParaView / ViennaPS-readable surface mesh
# result.plot("etch_hole.png")        # quick x-z cross-section through the feature centre
print("wrote etch_hole.vtk")
