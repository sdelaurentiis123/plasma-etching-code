#!/usr/bin/env python3
"""Assemble the de Boer AR>20 divergence curves (2026-07-07 campaign) from the arm npz files
(scripts/deboer_divergence_arm.py) and render viz/deboer_divergence.png.

Method notes (the two artifacts this fixes):
1. r0 NORMALIZATION BIAS: long runs record depth coarsely; the first record already sits at AR~3,
   so max-rate-at-AR<2 is understated -> the whole normalized curve inflates. Fix: r0 = median
   smoothed rate over 0.5<AR<2 from a FINE-cadence short run of the same config+throttle.
2. GRID QUANTIZATION: depth moves in dx=0.25 quanta -> raw np.gradient rates are a 6 um/t comb.
   Fix: boxcar the DEPTH history (n=9 records) before differentiating.

The long-run curve is only used ABOVE the fine run's reach; the fine curve (weighted like 2 seeds)
covers the low-AR knee. de Boer experiment: nr 1.0/0.43/0.29/0.20 @ AR 0/10/20/40.

Usage: deboer_divergence_plot.py [dir-with-arm_*.npz]   (expects arm_{k0,kt,d0,dt}_{fine,s0,s1}.npz)
"""
import os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SC = sys.argv[1] if len(sys.argv) > 1 else "."
EXP_AR = np.array([0., 10., 20., 40.]); EXP_R = np.array([1., 0.43, 0.29, 0.20])


def rate_curve(path, smooth_n=9):
    d = np.load(path); t, dd, W = d['t'], d['depth'], float(d['W'])
    if len(dd) > smooth_n:                       # boxcar the DEPTH history, then differentiate
        pad = smooth_n // 2
        dp = np.concatenate([dd[:1].repeat(pad), dd, dd[-1:].repeat(pad)])
        dd = np.convolve(dp, np.ones(smooth_n) / smooth_n, mode='valid')
    return dd / W, np.gradient(dd, t)


def curve(fine, longs, arN):
    arf, rf = rate_curve(os.path.join(SC, f"arm_{fine}.npz"))
    m = (arf > 0.4) & (arf < 2.0)
    r0 = float(np.median(rf[m])) if m.any() else float(rf.max())
    acc = np.zeros_like(arN); cnt = np.zeros_like(arN); amax = arf.max()
    s = np.interp(arN, arf, np.clip(rf / r0, 0, 1.4), left=1.0, right=np.nan)
    ok = ~np.isnan(s); acc[ok] += 2 * s[ok]; cnt[ok] += 2
    for tg in longs:
        arl, rl = rate_curve(os.path.join(SC, f"arm_{tg}.npz"))
        s = np.interp(arN, arl[2:], np.clip(rl[2:] / r0, 0, 1.4), left=np.nan, right=np.nan)
        s[arN <= arf.max() - 1.0] = np.nan       # long runs only contribute beyond the fine reach
        ok = ~np.isnan(s); acc[ok] += s[ok]; cnt[ok] += 1; amax = max(amax, arl.max())
    return np.where(cnt > 0, acc / np.maximum(cnt, 1), np.nan), float(amax), r0


arN = np.linspace(0.5, 46, 183)
out = {}
for name, fine, longs in [("KNEE", "k0_fine", ["k0_s0", "k0_s1"]),
                          ("KNEE+THR", "kt_fine", ["kt_s0", "kt_s1"]),
                          ("DEFAULT", "d0_fine", ["d0_s0", "d0_s1"]),
                          ("DEFAULT+THR", "dt_fine", ["dt_s0", "dt_s1"])]:
    try:
        out[name] = curve(fine, longs, arN)
    except FileNotFoundError:
        print(f"  ({name}: arms missing, skipped)")

probe = np.array([2, 5, 8, 10, 15, 20, 25, 30, 35, 40, 44])
exp = np.interp(probe, EXP_AR, EXP_R)
print(f"{'AR':>5} {'exp':>6}" + "".join(f"{k:>13}" for k in out))
for i, p in enumerate(probe):
    row = f"{p:5.0f} {exp[i]:6.3f}"
    for k, (nr, _, _) in out.items():
        ok = ~np.isnan(nr); row += f"{np.interp(p, arN[ok], nr[ok], right=np.nan):13.3f}"
    print(row)
for k, (nr, amax, r0) in out.items():
    ok = ~np.isnan(nr); a = arN[ok]; m = (a >= 5) & (a <= 44)
    e = np.interp(a[m], EXP_AR, EXP_R)
    print(f"  {k:12s} maxAR {amax:5.1f}  r0 {r0:6.2f} um/t  RMSE(AR5-44) {np.sqrt(np.mean((nr[ok][m]-e)**2)):.3f}")
np.savez(os.path.join(SC, "deboer_divergence.npz"), arN=arN, exp_ar=EXP_AR, exp_r=EXP_R,
         **{f"nr_{k.replace('+','_')}": v[0] for k, v in out.items()},
         **{f"amax_{k.replace('+','_')}": np.float64(v[1]) for k, v in out.items()})

# ---- figure ----
def cl(nr, n=7):
    ok = ~np.isnan(nr); ar, v = arN[ok], nr[ok]
    if len(v) >= n: v = np.convolve(np.clip(v, 0, 1.1), np.ones(n) / n, mode='same')
    return ar, np.clip(v, 0, 1.1)

fig, ax = plt.subplots(figsize=(9.6, 5.9))
xr = np.linspace(0, 46, 300)
ax.plot(xr, np.interp(xr, EXP_AR, EXP_R), 'k-', lw=2.8, label='de Boer / Blauw EXPERIMENT (sustains ~0.20 floor)')
ax.plot(EXP_AR, EXP_R, 'ko', ms=8, zorder=6)
STYLE = {"KNEE": ('-', '#9534c0', 2.3, 1.0), "KNEE+THR": ('--', '#0a8f3c', 2.3, 1.0),
         "DEFAULT": ('-', '#1f77b4', 1.6, 0.8), "DEFAULT+THR": ('--', '#e07b00', 1.6, 0.8)}
LBL = {"KNEE": lambda a: f"petch KNEE (matched process params)  reaches AR{a:.0f}",
       "KNEE+THR": lambda a: 'KNEE + charging floor throttle Q(AR)  [opt-in]',
       "DEFAULT": lambda a: 'petch DEFAULT (ViennaPS regime)',
       "DEFAULT+THR": lambda a: 'DEFAULT + throttle'}
for k, (nr, amax, _) in out.items():
    ls, c, lw, al = STYLE[k]
    ax.plot(*cl(nr), ls, color=c, lw=lw, alpha=al, label=LBL[k](amax))
ax.axvspan(10, 46, color='#c0392b', alpha=0.05)
ax.axvline(10, color='#c0392b', ls=':', lw=1.3)
ax.text(10.4, 1.00, 'divergence onset AR~10', color='#c0392b', fontsize=9)
ax.set_xlim(0, 46); ax.set_ylim(0, 1.09)
ax.set_xlabel('aspect ratio  (depth / width)', fontsize=12)
ax.set_ylabel('normalized bottom etch rate', fontsize=12)
ax.set_title('de Boer deep-trench ARDE to AR46: the SIMULATION floor collapses (all configs ~0.05-0.09 by AR30)\n'
             'while the EXPERIMENT sustains ~0.20; the charging throttle is inert at deep AR (neutral-limited floor)',
             fontsize=10.2)
ax.legend(loc='upper right', fontsize=9); ax.grid(alpha=0.3)
ax.annotate('experiment: factor 2.7 decay AR8->40', xy=(40, 0.21), xytext=(24, 0.44), fontsize=9,
            color='#333', arrowprops=dict(arrowstyle='->', color='#888'))
ax.annotate('simulation: factor ~11 decay AR8->40\n(floor collapse; throttle does not change it)',
            xy=(35, 0.06), xytext=(18, 0.62), fontsize=9, color='#6a1b9a',
            arrowprops=dict(arrowstyle='->', color='#9534c0'))
fig.tight_layout()
root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
outp = os.path.join(root, 'viz', 'deboer_divergence.png')
fig.savefig(outp, dpi=150)
print('wrote', outp)
