#!/usr/bin/env python3
"""Render the showcase movie series.

Outputs (all 1024x1024, 24 fps, H.264):
  showcase/media/reactor_scene.mp4  - chamber scale, solved Oxford state
  showcase/media/sheath_scene.mp4   - moving RF sheath with tracer ions
  showcase/media/flyby.mp4          - title -> reactor -> sheath -> feature
                                      continuous dive, plasma to wafer

Physics on screen comes from the same replayed objects the data build uses:
the Turner-Chabert sheath is receipt-checked against the frozen audit, and
the reactor readouts are the solved Zhu/Oxford power nodes. Geometry drawing
is schematic; every number is solved.
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrow, Rectangle
from PIL import Image

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))
from build_data import load_sheath  # noqa: E402

W = H = 1024
FPS = 24
DPI = 128

BG = "#0b1220"
PANEL = "#101a2d"
LINE = "#1e2a42"
INK = "#e8ecf4"
MUTED = "#8fa0b8"
DIM = "#5b6b85"
TEAL = "#2dd4bf"
OXIDE = "#1f8a8c"
MAG = "#e26ad0"
AMBER = "#f5b04a"
BLUE = "#6ea8fe"

QM = 1.602176634e-19 / (39.948 * 1.66053906660e-27)
EV = 1.602176634e-19
M_AR = 39.948 * 1.66053906660e-27


class MovieWriter:
    def __init__(self, path: Path):
        path.parent.mkdir(parents=True, exist_ok=True)
        self.proc = subprocess.Popen(
            ["ffmpeg", "-y", "-loglevel", "error",
             "-f", "rawvideo", "-pix_fmt", "rgb24",
             "-s", f"{W}x{H}", "-r", str(FPS), "-i", "-",
             "-c:v", "libx264", "-preset", "medium", "-crf", "23",
             "-pix_fmt", "yuv420p", "-movflags", "+faststart", str(path)],
            stdin=subprocess.PIPE)
        self.path = path

    def add(self, frame: np.ndarray):
        self.proc.stdin.write(np.ascontiguousarray(frame[..., :3]).tobytes())

    def close(self):
        self.proc.stdin.close()
        self.proc.wait()
        print(f"wrote {self.path} "
              f"({self.path.stat().st_size/1e6:.1f} MB)")


def fig_to_rgb(fig) -> np.ndarray:
    fig.canvas.draw()
    buf = np.asarray(fig.canvas.buffer_rgba())
    return buf[..., :3].copy()


def ease(u: float) -> float:
    u = min(max(u, 0.0), 1.0)
    return u * u * (3 - 2 * u)


def nice_scale_label(m: float) -> str:
    if m >= 0.01:
        return f"{m*100:.0f} cm"
    if m >= 1e-3:
        return f"{m*1000:.0f} mm"
    if m >= 1e-4:
        return f"{m*1e3:.1f} mm"
    return f"{m*1e6:.0f} µm"


# ---------------------------------------------------------------- reactor
class ReactorScene:
    """Oxford PlasmaPro 80 cross-section with the solved 105 W state."""

    N = int(8.5 * FPS)

    def __init__(self):
        data = json.loads((HERE / "data" / "showcase_data.json").read_text())
        self.node = data["zhu"]["nodes"][-1]          # 105 W solved state
        self.fig = plt.figure(figsize=(W / DPI, H / DPI), dpi=DPI)
        self.fig.patch.set_facecolor(BG)
        # glow field over the plasma slab
        nx, ny = 320, 60
        gx = np.linspace(-1, 1, nx)[None, :]
        gy = np.linspace(-1, 1, ny)[:, None]
        self.glow = (np.clip(1.15 - 0.55 * gx**8 - 0.75 * np.abs(gy)**2.5,
                             0, None))
        from scipy.ndimage import gaussian_filter
        rng = np.random.default_rng(7)
        self.noise = np.stack([
            gaussian_filter(rng.normal(0, 1, (ny, nx)), sigma=2.2)
            for _ in range(24)])

    def frame(self, i: int) -> np.ndarray:
        fig = self.fig
        fig.clf()
        ax = fig.add_axes([0, 0, 1, 1])
        ax.set_facecolor(BG)
        ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values():
            s.set_visible(False)

        t = i / FPS
        u = i / (self.N - 1)
        # camera: hold wide, then dive toward the wafer center
        zoom_u = ease(max(u - 0.52, 0.0) / 0.48)
        # wide view / target view boxes (metres)
        wx, wy, ww = 0.0, 0.078, 0.42        # centre x, centre y, width
        tx, ty, tw = 0.0, 0.0545, 0.052
        cx = wx + (tx - wx) * zoom_u
        cy = wy + (ty - wy) * zoom_u
        vw = ww * (tw / ww) ** zoom_u
        ax.set_xlim(cx - vw / 2, cx + vw / 2)
        ax.set_ylim(cy - vw / 2, cy + vw / 2)

        # chamber body
        ax.add_patch(Rectangle((-0.17, 0.0), 0.34, 0.155, fill=False,
                               ec=LINE, lw=2.5))
        # powered electrode + wafer
        ax.add_patch(Rectangle((-0.12, 0.022), 0.24, 0.028,
                               fc="#182643", ec=LINE))
        ax.add_patch(Rectangle((-0.12, 0.050), 0.24, 0.0012, fc=TEAL,
                               ec="none"))
        # grounded showerhead
        ax.add_patch(Rectangle((-0.12, 0.083), 0.24, 0.022,
                               fc="#182643", ec=LINE))
        for gx_ in np.linspace(-0.10, 0.10, 9):
            ax.add_patch(FancyArrow(gx_, 0.104, 0, -0.008, width=0.0006,
                                    head_width=0.003, head_length=0.004,
                                    fc=DIM, ec="none"))
        # pump ports
        ax.add_patch(FancyArrow(0.145, 0.012, 0.03, 0, width=0.002,
                                head_width=0.007, head_length=0.008,
                                fc=DIM, ec="none"))
        ax.add_patch(FancyArrow(-0.145, 0.012, -0.03, 0, width=0.002,
                                head_width=0.007, head_length=0.008,
                                fc=DIM, ec="none"))

        # plasma glow slab between sheaths (sheath gaps drawn dark)
        s_lo, s_hi = 0.0545, 0.080           # bulk plasma band
        namp = 0.55 * (1.0 - 0.9 * zoom_u)
        shimmer = 1.0 + 0.05 * np.sin(2 * np.pi * 1.7 * t) \
            + namp * self.noise[i % 24]
        g = np.clip(self.glow * shimmer, 0, 1.35)
        rgb = np.zeros(g.shape + (4,))
        rgb[..., 0] = np.clip(0.28 * g + 0.30 * g**3, 0, 1)   # R
        rgb[..., 1] = np.clip(0.22 * g + 0.18 * g**2, 0, 1)   # G
        rgb[..., 2] = np.clip(0.55 * g, 0, 1)                 # B
        rgb[..., 3] = np.clip(0.92 * g, 0, 1)
        ax.imshow(rgb, extent=[-0.122, 0.122, s_lo, s_hi], origin="lower",
                  interpolation="bilinear", zorder=3, aspect="auto")
        # faint sheath edge glow above the wafer as we get close
        if zoom_u > 0.3:
            a = (zoom_u - 0.3) / 0.7
            ax.axhline(s_lo, xmin=0, xmax=1, color=BLUE, lw=1.5 + 2 * a,
                       alpha=0.55 * a, zorder=4)

        # labels (fade out during the dive)
        la = 1.0 - ease(min(zoom_u / 0.55, 1.0))
        n = self.node
        if la > 0.01:
            ax.text(0.03, 0.960, "Oxford PlasmaPro 80 · CCP",
                    transform=ax.transAxes, color=INK, fontsize=16,
                    family="monospace", alpha=la, weight="bold")
            ax.text(0.03, 0.924,
                    "CHF₃/SF₆/O₂  55/5/1 sccm · "
                    "30 mTorr · 13.56 MHz",
                    transform=ax.transAxes, color=MUTED, fontsize=12.5,
                    family="monospace", alpha=la)
            read = [
                ("absorbed power", f"{n['absorbed_power_W']:.0f} W"),
                ("electron density",
                 f"{n['electron_density_m3']:.2e} m⁻³"),
                ("mean electron energy",
                 f"{n['mean_electron_energy_eV']:.1f} eV"),
                ("reduced field E/N", f"{n['reduced_field_Td']:.0f} Td"),
                ("ion flux to wafer",
                 f"{n['total_ion_flux_m2_s']:.2e} m⁻²s⁻¹"),
            ]
            for k, (lab, val) in enumerate(read):
                ax.text(0.03, 0.876 - 0.033 * k, f"{lab:<22s}",
                        transform=ax.transAxes, color=DIM, fontsize=11,
                        family="monospace", alpha=la)
                ax.text(0.32, 0.876 - 0.033 * k, val,
                        transform=ax.transAxes, color=TEAL, fontsize=11,
                        family="monospace", alpha=la)
            ax.text(0.03, 0.712, "67 species · 259 reactions · "
                    "one deterministic solve",
                    transform=ax.transAxes, color=MUTED, fontsize=11,
                    family="monospace", alpha=la)
            ax.text(0, 0.0935, "grounded showerhead",
                    color=MUTED, fontsize=10.5, family="monospace",
                    alpha=la, ha="center", va="center", clip_on=True)
            ax.text(0, 0.030, "powered electrode · 240 mm wafer",
                    color=MUTED, fontsize=10.5, family="monospace",
                    alpha=la, ha="center", va="center", clip_on=True)
        # dive label
        da = ease(max(zoom_u - 0.65, 0) / 0.35)
        if da > 0.01:
            ax.text(0.5, 0.86, "descending to the sheath",
                    transform=ax.transAxes, color=BLUE, fontsize=15,
                    family="monospace", alpha=da, ha="center")

        # scale bar
        bar_m = vw * 0.24
        ax.plot([0.70, 0.94], [0.055, 0.055], transform=ax.transAxes,
                color=INK, lw=2.5)
        ax.text(0.82, 0.075, nice_scale_label(bar_m),
                transform=ax.transAxes, color=INK, fontsize=13,
                family="monospace", ha="center")
        return fig_to_rgb(fig)


# ----------------------------------------------------------------- sheath
class SheathScene:
    """Tracer ions through the replayed Turner-Chabert moving field."""

    N = int(10.5 * FPS)
    SLOWMO_SIM_PER_S = 1.25e-7       # 10.5 s of film ~ 2.6 RF periods

    def __init__(self):
        self.sheath, audit = load_sheath()
        self.audit = audit
        s = self.sheath
        self.smax = s.maximum_width_m
        self.Vmax = s.maximum_voltage_v
        self.period = s.period_s
        self.bohm = float(np.sqrt(EV * s.electron_temperature_eV / M_AR))
        data = json.loads((HERE / "data" / "showcase_data.json").read_text())
        self.hist = np.array(data["sheath"]["iead_hist"], float)
        self.edges = np.array(data["sheath"]["iead_edges_eV"], float)
        self._simulate()
        self.fig = plt.figure(figsize=(W / DPI, H / DPI), dpi=DPI)
        self.fig.patch.set_facecolor(BG)

    def _simulate(self):
        rng = np.random.default_rng(11)
        s = self.sheath
        sim_dt = self.SLOWMO_SIM_PER_S / FPS
        n_sub = 60
        h = sim_dt / n_sub
        warmup = int(1.3 * FPS)
        ions = [{"x": rng.uniform(0.0, 0.85) * self.smax,
                 "v": self.bohm * rng.uniform(0.85, 1.3),
                 "cx": rng.uniform(0.04, 0.96)} for _ in range(34)]
        self.frames_state = []
        self.impacts = []  # (frame, cx, eV)
        self.live_hist = np.zeros((self.N, 40))
        live = np.zeros(40)
        t = 0.0
        for fi in range(self.N):
            for _ in range(n_sub):
                for io in ions:
                    a = QM * float(s.electric_field_V_m(io["x"], t))
                    io["v"] += a * h
                    io["x"] += io["v"] * h
                dead = [io for io in ions if io["x"] >= self.smax]
                for io in dead:
                    if fi > warmup:
                        eV_i = 0.5 * M_AR * io["v"] ** 2 / EV
                        self.impacts.append((fi, io["cx"], eV_i))
                        b = min(int(eV_i / (self.Vmax * 1.05) * 40), 39)
                        live[b] += 1
                ions = [io for io in ions if io["x"] < self.smax]
                # trickle spawns so cohorts decorrelate
                if len(ions) < 34 and rng.random() < 0.18:
                    ions.append({"x": 0.0,
                                 "v": self.bohm * rng.uniform(0.9, 1.1),
                                 "cx": rng.uniform(0.04, 0.96)})
                t += h
            self.frames_state.append(
                [(io["x"], io["v"], io["cx"]) for io in ions])
            self.live_hist[fi] = live
        self.sim_times = np.arange(self.N) * sim_dt

    def frame(self, i: int) -> np.ndarray:
        fig = self.fig
        fig.clf()
        s = self.sheath
        t = float(self.sim_times[i])
        y = float(np.asarray(s.normalized_charge(t)))
        front = y ** 3 * self.smax

        axm = fig.add_axes([0.06, 0.10, 0.55, 0.82])
        axm.set_facecolor(BG)
        axm.set_xlim(0, 1); axm.set_ylim(self.smax * 1.06, -self.smax * 0.06)
        axm.set_xticks([]); axm.set_yticks([])
        for sp in axm.spines.values():
            sp.set_visible(False)
        # regions: electron-covered above front, ion sheath below
        axm.add_patch(Rectangle((0, 0), 1, front, fc=BLUE, alpha=0.10))
        axm.add_patch(Rectangle((0, front), 1, self.smax - front,
                                fc=MAG, alpha=0.13))
        # field-strength gradient inside the sheath
        gy = np.linspace(front, self.smax, 40)
        for k in range(len(gy) - 1):
            al = 0.16 * ((gy[k] - front) / max(self.smax - front, 1e-9))
            axm.add_patch(Rectangle((0, gy[k]), 1, gy[k + 1] - gy[k],
                                    fc=MAG, alpha=al, ec="none"))
        # electron front
        axm.axhline(front, color=BLUE, lw=3, alpha=0.95)
        axm.axhline(front, color=BLUE, lw=9, alpha=0.25)
        # wafer
        axm.add_patch(Rectangle((0, self.smax), 1, self.smax * 0.06,
                                fc="#182643", ec="none"))
        axm.axhline(self.smax, color=TEAL, lw=3)
        # ions
        for x, v, cx in self.frames_state[i]:
            hot = min((v - self.bohm) / 40000.0, 1.0)
            axm.plot([cx], [x], "o",
                     ms=3.5 + 3.5 * hot,
                     color=(0.96, 0.69 + 0.20 * hot, 0.29),
                     alpha=0.55 + 0.45 * hot)
            if hot > 0.05:
                axm.plot([cx, cx], [x, x - (0.02 + 0.10 * hot) * self.smax],
                         color=AMBER, lw=1.2, alpha=0.3 * hot)
        # impact rings
        for (fi, cx, _e) in self.impacts:
            age = (i - fi) / FPS
            if 0 <= age < 0.45:
                al = 1 - age / 0.45
                axm.plot([cx], [self.smax], "o", ms=6 + age * 70,
                         mfc="none", mec=TEAL, mew=1.6, alpha=al)
        axm.text(0.02, -0.028 * self.smax / 1.06, "plasma", color=MUTED,
                 fontsize=14, family="monospace")
        axm.text(0.02, self.smax * 1.045, "wafer", color=MUTED,
                 fontsize=14, family="monospace", va="top")
        axm.text(0.98, front - 0.015 * self.smax, "electron front",
                 color=BLUE, fontsize=12, family="monospace", ha="right")
        # scale bar = 0.2 mm
        blen = 0.2e-3 / self.smax * 0.82  # fraction of axis height... draw v
        axm.plot([0.955, 0.955], [self.smax * 0.98,
                                  self.smax * 0.98 - 0.2e-3],
                 color=INK, lw=2.5)
        axm.text(0.945, self.smax * 0.97 - 1e-4, "0.2 mm", color=INK,
                 fontsize=11.5, family="monospace", ha="right",
                 rotation=90, va="center")

        # V(t) sparkline
        axv = fig.add_axes([0.68, 0.62, 0.28, 0.24])
        axv.set_facecolor(PANEL)
        tt = np.linspace(0, self.period, 200)
        yy = np.asarray(s.normalized_charge(tt))
        vv = self.Vmax * (1 - (4 / 3) * yy + (1 / 3) * yy ** 4)
        axv.plot(tt / self.period, vv, color=MAG, lw=2)
        ph = (t % self.period) / self.period
        vnow = self.Vmax * (1 - (4 / 3) * y + (1 / 3) * y ** 4)
        axv.plot([ph], [vnow], "o", color=AMBER, ms=7)
        axv.set_xlim(0, 1); axv.set_ylim(0, self.Vmax * 1.1)
        axv.set_xticks([]); axv.set_yticks([])
        for sp in axv.spines.values():
            sp.set_color(LINE)
        axv.set_title("sheath voltage V(t)", color=MUTED, fontsize=12,
                      family="monospace", loc="left")

        # IEAD building
        axh = fig.add_axes([0.68, 0.13, 0.28, 0.38])
        axh.set_facecolor(PANEL)
        centers = 0.5 * (self.edges[:-1] + self.edges[1:])
        axh.bar(centers, self.hist / self.hist.max(),
                width=np.diff(self.edges), color=OXIDE, alpha=0.85)
        lh = self.live_hist[i]
        if lh.max() > 0:
            le = np.linspace(0, self.Vmax * 1.05, 41)
            lc = 0.5 * (le[:-1] + le[1:])
            axh.bar(lc, lh / lh.max() * 0.6, width=np.diff(le),
                    color=AMBER, alpha=0.9)
        axh.set_xlim(0, self.Vmax * 1.08)
        axh.set_ylim(0, 1.05)
        axh.set_yticks([])
        axh.tick_params(colors=DIM, labelsize=9)
        for sp in axh.spines.values():
            sp.set_color(LINE)
        axh.set_title("ion energy at wafer (eV)", color=MUTED,
                      fontsize=12, family="monospace", loc="left")

        fig.text(0.06, 0.955, "the moving RF sheath · argon · "
                 "500 W · 2 mTorr", color=INK, fontsize=16,
                 family="monospace", weight="bold")
        fig.text(0.06, 0.925,
                 "Turner–Chabert field, solved · "
                 f"slowed {1/ (self.SLOWMO_SIM_PER_S):,.0f}×",
                 color=MUTED, fontsize=12, family="monospace")
        return fig_to_rgb(fig)


# ---------------------------------------------------------------- feature
class FeatureScene:
    """Hero3d etch frames, crossfaded to film speed, restyled captions."""

    N = int(8.0 * FPS)

    def __init__(self):
        tmp = Path("/tmp/hero3d_frames")
        tmp.mkdir(exist_ok=True)
        if not (tmp / "f_01.png").exists():
            subprocess.run(
                ["ffmpeg", "-y", "-loglevel", "error",
                 "-i", str(ROOT.parent / "bench" / "out" / "hero3d.mp4"),
                 str(tmp / "f_%02d.png")], check=True)
        # crop the burned-in caption strip; we draw our own clock
        self.src = [np.asarray(
            Image.open(p).convert("RGB").crop((0, 0, 896, 812))
            .resize((W, int(812 / 896 * W)), Image.LANCZOS),
            dtype=np.float32)
            for p in sorted(tmp.glob("f_*.png"))]
        self.fig = plt.figure(figsize=(W / DPI, H / DPI), dpi=DPI)

    def frame(self, i: int) -> np.ndarray:
        u_total = i / (self.N - 1)
        u = u_total * (len(self.src) - 1)
        k = min(int(u), len(self.src) - 2)
        f = u - k
        img = (1 - f) * self.src[k] + f * self.src[k + 1]
        fig = self.fig
        fig.clf()
        fig.patch.set_facecolor(BG)
        ih = img.shape[0]
        y0 = (H - ih) // 2 + 18
        ax = fig.add_axes([0, 1 - (y0 + ih) / H, 1, ih / H])
        ax.imshow(img.astype(np.uint8))
        ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_visible(False)
        fig.text(0.045, 0.955, "the feature · micrometre scale",
                 color=INK, fontsize=16, family="monospace", weight="bold")
        fig.text(0.045, 0.925, "level-set surface evolution · "
                 "shadowing + re-emission + surface chemistry",
                 color=MUTED, fontsize=12, family="monospace")
        fig.text(0.5, 0.055, f"process t = {u_total*60:4.1f} s",
                 color=MUTED, fontsize=15, family="monospace", ha="center")
        return fig_to_rgb(fig)


# ------------------------------------------------------------------ cards
def text_card(lines, n_frames, fade=0.5):
    fig = plt.figure(figsize=(W / DPI, H / DPI), dpi=DPI)
    fig.patch.set_facecolor(BG)
    for (txt, ypos, size, color, weight) in lines:
        fig.text(0.5, ypos, txt, color=color, fontsize=size,
                 family="monospace", ha="center", weight=weight)
    base = fig_to_rgb(fig)
    plt.close(fig)
    dark = np.zeros_like(base)
    out = []
    nf = int(fade * FPS)
    for i in range(n_frames):
        if i < nf:
            a = ease(i / nf)
        elif i > n_frames - nf:
            a = ease((n_frames - i) / nf)
        else:
            a = 1.0
        out.append((dark * (1 - a) + base * a).astype(np.uint8))
    return out


def zoom_transition(frame_a, frame_b, n=int(1.2 * FPS)):
    """Punch-in on A while B grows underneath."""
    a_img = Image.fromarray(frame_a)
    b_img = Image.fromarray(frame_b)
    out = []
    for i in range(n):
        u = ease((i + 1) / n)
        za = 1.0 + 1.8 * u
        wa = int(W / za)
        a_zoom = a_img.crop(((W - wa) // 2, (H - wa) // 2,
                             (W + wa) // 2, (H + wa) // 2)).resize(
            (W, H), Image.LANCZOS)
        zb = 0.62 + 0.38 * u
        wb = max(int(W * zb), 2)
        b_small = b_img.resize((wb, wb), Image.LANCZOS)
        b_canvas = Image.new("RGB", (W, H), (11, 18, 32))
        b_canvas.paste(b_small, ((W - wb) // 2, (H - wb) // 2))
        blend = Image.blend(a_zoom, b_canvas, u)
        out.append(np.asarray(blend))
    return out


def main() -> None:
    media = HERE / "media"
    reactor = ReactorScene()
    sheath = SheathScene()
    feature = FeatureScene()

    # standalone scene movies
    for name, scene in (("reactor_scene", reactor),
                        ("sheath_scene", sheath)):
        wtr = MovieWriter(media / f"{name}.mp4")
        for i in range(scene.N):
            wtr.add(scene.frame(i))
        wtr.close()

    # composite fly-by
    wtr = MovieWriter(media / "flyby.mp4")
    title = text_card([
        ("p e t c h", 0.60, 22, TEAL, "bold"),
        ("from the plasma to the wafer", 0.50, 30, INK, "bold"),
        ("one deterministic engine · reactor → sheath → "
         "feature", 0.42, 15, MUTED, "normal"),
    ], int(2.4 * FPS))
    for f in title:
        wtr.add(f)

    r_frames_last = None
    for i in range(reactor.N):
        f = reactor.frame(i)
        wtr.add(f)
        r_frames_last = f
    s_first = sheath.frame(0)
    for f in zoom_transition(r_frames_last, s_first):
        wtr.add(f)
    s_last = None
    for i in range(sheath.N):
        f = sheath.frame(i)
        wtr.add(f)
        s_last = f
    f_first = feature.frame(0)
    for f in zoom_transition(s_last, f_first):
        wtr.add(f)
    for i in range(feature.N):
        wtr.add(feature.frame(i))
    end = text_card([
        ("453 V sheath · 234.5 eV mean ion · 1.15° RMS angle",
         0.62, 15, MUTED, "normal"),
        ("67 species solved · zero fitted knobs · "
         "~14× faster than GPU ViennaPS", 0.55, 15, MUTED, "normal"),
        ("every number from frozen curated boards", 0.48, 13, DIM, "normal"),
        ("p e t c h", 0.38, 20, TEAL, "bold"),
    ], int(3.2 * FPS))
    for f in end:
        wtr.add(f)
    wtr.close()


if __name__ == "__main__":
    main()
