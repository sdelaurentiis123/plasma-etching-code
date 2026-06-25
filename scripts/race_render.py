#!/usr/bin/env python3
"""Render the petch-vs-ViennaPS real-time etch race from race.pkl: two panels etching the same hole,
played in WALL-CLOCK time (scaled to a watchable length but preserving the true speed ratio), so petch
finishes while ViennaPS is still grinding. Writes viz/race.gif + viz/race.mp4."""
import os, sys, pickle, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import animation
from scipy.ndimage import gaussian_filter, label

D = pickle.load(open(sys.argv[1] if len(sys.argv) > 1 else "race.pkl", "rb"))
W, sub, EXT = D['W'], D['sub_top'], D['EXT']
xs, zs = D['xs'], D['zs']
pw, vw = D['petch_wall'], D['vps_wall']
MOVIE_LEN, FPS = 9.0, 18                                  # scaled movie seconds; ratio preserved
nfr = int(MOVIE_LEN * FPS)


def petch_gas(phi):
    sm = gaussian_filter(phi, 0.8); gas = sm < 0
    lbl, n = label(gas)
    if n == 0:
        return gas.astype(float)
    keep = set(int(v) for v in np.unique(lbl[:, -1]) if v > 0)
    return np.isin(lbl, list(keep)).astype(float)


def vps_depthfill(x, z, xg):
    """ViennaPS centre-slice surface nodes -> filled cavity: floor depth per x-bin."""
    depth = np.zeros_like(xg)                              # depth below surface (>=0)
    for i, xc in enumerate(xg):
        m = np.abs(x - xc) < D['DX']
        if m.any():
            depth[i] = max(0.0, -z[m].min())               # deepest node near this x (z<0 = etched)
    return depth


Xp = np.meshgrid(xs - EXT/2, sub - zs, indexing='ij')      # petch (x, depth)
Xpx, Dp = Xp[0], Xp[1]
xg = np.linspace(-EXT/2, EXT/2, 160)
fig, (axp, axv) = plt.subplots(1, 2, figsize=(9.5, 5.6), sharey=True)
dmax = sub


def frame(k):
    t_movie = k / FPS
    real_t = t_movie / MOVIE_LEN * vw                      # wall-clock seconds this frame represents
    for ax, wall, label_ in [(axp, pw, "petch (GPU)"), (axv, vw, "ViennaPS (GPU)")]:
        ax.clear(); ax.set_facecolor("#9a8c7a")
        prog = min(real_t / wall, 1.0)
        if ax is axp:
            i = min(int(prog * (len(D['petch_frames']) - 1)), len(D['petch_frames']) - 1)
            gas = petch_gas(D['petch_frames'][i]['phi'])
            ax.contourf(Xpx, Dp, gas, levels=[0.5, 1.5], colors=["#bfe0ff"])
            ax.contour(Xpx, Dp, gas, levels=[0.5], colors="k", linewidths=1.3)
            dep = D['petch_frames'][i]['depth']
        else:
            i = min(int(prog * (len(D['vps_frames']) - 1)), len(D['vps_frames']) - 1)
            fr, x, z = D['vps_frames'][i]
            dep_prof = vps_depthfill(x, z, xg)
            ax.fill_between(xg, 0, dep_prof, color="#bfe0ff", zorder=2)
            ax.plot(xg, dep_prof, "k", lw=1.3, zorder=3)
            dep = float(dep_prof.max())
        ax.axhline(0, color="0.3", lw=0.8, ls="--")
        ax.set_xlim(-EXT/2, EXT/2); ax.set_ylim(0, dmax); ax.invert_yaxis(); ax.set_aspect("equal")
        ax.set_xlabel("x (µm)")
        done = " ✓ DONE" if prog >= 1.0 else ""
        ax.set_title(f"{label_}{done}\ndepth {dep:.2f} µm   AR {dep/W:.1f}", fontsize=11,
                     color=("#1a7a3a" if prog >= 1.0 else "black"))
    axp.set_ylabel("depth below surface (µm)")
    fig.suptitle(f"Same hole, both on RTX 3090 — wall clock {real_t:5.1f} s     "
                 f"(petch {pw:.1f}s  vs  ViennaPS {vw:.1f}s  =  {vw/pw:.0f}× faster)",
                 fontsize=12, fontweight="bold")


anim = animation.FuncAnimation(fig, frame, frames=nfr, interval=1000/FPS)
out = os.path.join(os.path.dirname(__file__), "..", "viz"); os.makedirs(out, exist_ok=True)
gif = os.path.join(out, "race.gif"); mp4 = os.path.join(out, "race.mp4")
anim.save(gif, writer=animation.PillowWriter(fps=FPS)); print("wrote", gif)
try:
    anim.save(mp4, writer=animation.FFMpegWriter(fps=FPS, bitrate=2200)); print("wrote", mp4)
except Exception as e:
    print("mp4 skipped:", e)
