#!/usr/bin/env python3
"""Etch-evolution movies: watch a trench and a hole deepen over time (centre x-z cross-section).
Runs petch locally (CPU ok), captures frames, writes GIF + MP4 to viz/. Pure human-facing viz."""
import os, sys
os.environ.setdefault("PETCH_DEVICE", "cpu")
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import animation
from scipy.ndimage import gaussian_filter, label
import petch
from petch import threed as t3


def clean_gas(phi, sigma=0.8):
    """Display-only cleanup: smooth the level set, keep only the etched (gas) region connected to the
    open top -- removes MC-noise specks and detached filaments so the profile reads cleanly."""
    sm = gaussian_filter(phi, sigma)
    gas = sm < 0
    lbl, n = label(gas)
    if n == 0:
        return gas.astype(float), sm
    top = lbl[:, -1]                                    # labels touching the open top (last z column)
    keep = set(int(v) for v in np.unique(top) if v > 0)
    mask = np.isin(lbl, list(keep)) if keep else (lbl == np.bincount(lbl.flat)[1:].argmax() + 1)
    return mask.astype(float), sm

OUT = os.path.join(os.path.dirname(__file__), "..", "viz")
os.makedirs(OUT, exist_ok=True)


def run(hole, dx=0.10, W=1.0, sub=10.0, t_end=2.6, steps=52):
    GEO = dict(Lx=3.6, Ly=(3.6 if hole else 0.6), Lz=2*dx+sub+0.4, dx=dx, trench_width=W,
               mask_th=2*dx, sub_top=sub, hole=hole)
    p = dict(petch.PAR); p['rate_scale'] = 0.12
    if not hole:
        p['periodic_y'] = 1
    fl = petch.Flags(chemistry="belen", yield_angular="viennaps", coverage_sticking=True,
                     warm_start_coverage=True, sampling="sobol", ion_reflection=True)
    g = t3.run_etch_3d(t_end=t_end, n_steps=steps, par=p, flags=fl, n_ion=14000, n_neu=14000,
                       reinit_method="fsm", verbose=True, record_depth_every=1, record_frames=True, **GEO)
    return g


def animate(g, title, stub):
    frames = g['frames']; xs, zs = g['xs'], g['zs']; W = g['trench_width']; sub = g['sub_top']
    X = np.meshgrid(xs - g['Lx']/2, zs, indexing='ij')[0]       # centre x on the feature
    D = np.meshgrid(xs - g['Lx']/2, sub - zs, indexing='ij')[1]  # depth below the wafer surface
    dmax = sub                                                 # show full etchable depth
    fig, ax = plt.subplots(figsize=(5.0, 5.8))
    def draw(i):
        ax.clear()
        f = frames[i]; gas, sm = clean_gas(f['phi_xz'])        # cleaned gas mask
        ax.set_facecolor("#9a8c7a")                            # solid silicon background
        ax.contourf(X, D, gas, levels=[0.5, 1.5], colors=["#bfe0ff"])                        # etched gas
        ax.contour(X, D, gas, levels=[0.5], colors="k", linewidths=1.4)                      # surface
        ax.axhline(0.0, color="0.3", lw=0.8, ls="--")          # wafer surface (depth 0)
        ax.set_xlim(-g['Lx']/2, g['Lx']/2); ax.set_ylim(0, dmax); ax.invert_yaxis()         # surface at top
        ax.set_aspect("equal"); ax.set_xlabel("x (µm)"); ax.set_ylabel("depth below surface (µm)")
        ar = f['depth'] / W
        ax.set_title(f"{title}\nt = {f['t']:.2f} min    depth {f['depth']:.2f} µm    AR {ar:.1f}", fontsize=11)
    anim = animation.FuncAnimation(fig, draw, frames=len(frames), interval=120)
    gif = os.path.join(OUT, f"{stub}.gif"); mp4 = os.path.join(OUT, f"{stub}.mp4")
    anim.save(gif, writer=animation.PillowWriter(fps=8)); print("wrote", gif, flush=True)
    try:
        anim.save(mp4, writer=animation.FFMpegWriter(fps=12, bitrate=1800)); print("wrote", mp4, flush=True)
    except Exception as e:
        print("mp4 skipped:", e, flush=True)
    plt.close(fig)
    filmstrip(g, title, stub)


def filmstrip(g, title, stub):
    """Static gallery: the etched shape at 5 evenly-spaced times (for docs / print)."""
    frames = g['frames']; xs, zs = g['xs'], g['zs']; W = g['trench_width']; sub = g['sub_top']
    X = np.meshgrid(xs - g['Lx']/2, zs, indexing='ij')[0]
    D = np.meshgrid(xs - g['Lx']/2, sub - zs, indexing='ij')[1]
    idx = np.linspace(0, len(frames)-1, 5).astype(int)
    fig, axes = plt.subplots(1, 5, figsize=(13, 4.4), sharey=True)
    for ax, i in zip(axes, idx):
        f = frames[i]; gas, _ = clean_gas(f['phi_xz'])
        ax.set_facecolor("#9a8c7a")
        ax.contourf(X, D, gas, levels=[0.5, 1.5], colors=["#bfe0ff"])
        ax.contour(X, D, gas, levels=[0.5], colors="k", linewidths=1.2)
        ax.axhline(0.0, color="0.3", lw=0.7, ls="--")
        ax.set_xlim(-g['Lx']/2, g['Lx']/2); ax.set_ylim(0, sub); ax.invert_yaxis()
        ax.set_aspect("equal"); ax.set_xlabel("x (µm)")
        ax.set_title(f"t={f['t']:.1f}m  AR {f['depth']/W:.1f}", fontsize=10)
    axes[0].set_ylabel("depth below surface (µm)")
    fig.suptitle(title, fontsize=13, fontweight="bold")
    plt.tight_layout()
    p = os.path.join(OUT, f"{stub}_filmstrip.png"); plt.savefig(p, dpi=140); print("wrote", p, flush=True)
    plt.close(fig)


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "both"
    if which in ("trench", "both"):
        print("=== trench etch ===", flush=True)
        animate(run(False), "SF6/O2 trench etch (W=1.0 µm)", "etch_trench")
    if which in ("hole", "both"):
        print("=== hole etch ===", flush=True)
        animate(run(True), "SF6/O2 hole etch (Ø=1.0 µm)", "etch_hole")
    print("done", flush=True)
