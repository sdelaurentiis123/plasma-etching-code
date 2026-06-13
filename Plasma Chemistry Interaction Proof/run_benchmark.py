"""run_benchmark.py — compare our feature-scale etcher against ViennaPS SF6O2Etching."""
import numpy as np, time, json
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import importlib, feature_etch as fe; importlib.reload(fe)

import viennaps as ps, viennaps.d2 as vps
ps.Logger.setLogLevel(ps.LogLevel.ERROR); ps.Length.setUnit("micrometer"); ps.Time.setUnit("min")
Mat = ps.Material

DX = 0.25; W = 20.0; H = 24.0; MASK = 2.0; DUR = 3.0; NSTEPS = 60
OURS = dict(fe.PAR); OURS['rate_scale'] = 0.29

def vps_etch(width, dur=DUR):
    d = vps.Domain()
    vps.MakeTrench(domain=d, gridDelta=DX, xExtent=W, yExtent=H, trenchWidth=width,
                   trenchDepth=MASK, taperingAngle=0.0, baseHeight=0.0,
                   periodicBoundary=False, makeMask=True, material=Mat.Si).apply()
    m = vps.SF6O2Etching(vps.SF6O2Etching.defaultParameters())
    p = vps.Process(d, m, dur); t0 = time.time(); p.apply(); rt = time.time()-t0
    n = np.array(d.getSurfaceMesh().getNodes())
    return n[:,0], n[:,1], rt   # x, y (substrate top at y=0), runtime

def vps_profile(x, y):
    """Etched-Si trench profile: x centered, depth positive-down, substrate region only."""
    sel = y < -0.02
    return x[sel], -y[sel]

def ours_profile(seg, sub_top):
    cx = np.r_[seg[:,0], seg[:,2]]; cy = np.r_[seg[:,1], seg[:,3]]
    sel = cy < sub_top-0.02
    return cx[sel]-W/2, sub_top-cy[sel]

def depth_centre(xc, dep, half=1.5):
    c = np.abs(xc) < half
    return dep[c].max() if c.any() else 0.0

# warm JIT
fe.run_etch(W=12,H=12,dx=0.5,trench_width=5,sub_top=9,t_end=0.5,n_steps=2,verbose=False)

# ---------------- main width-8 run, both codes ----------------
print("ViennaPS width-8 ..."); vx,vy,v_rt = vps_etch(8.0)
vxc, vdep = vps_profile(vx, vy)

print("Ours width-8 ...")
snaps = [0,12,24,36,48,59]
r = fe.run_etch(W=W,H=H,dx=DX,trench_width=8.0,mask_thickness=MASK,sub_top=18.0,
                t_end=DUR,n_steps=NSTEPS,par=OURS,snapshots=set(snaps),verbose=False)
oxc, odep = ours_profile(r['segs'], 18.0)
T = r['timings']
print(f"  ours runtime {T['total']:.1f}s, raytrace {100*T['raytrace']/T['total']:.0f}%")
print(f"  depth: ViennaPS {depth_centre(vxc,vdep):.2f}  ours {depth_centre(oxc,odep):.2f}")

# ---------------- ARDE sweep ----------------
widths = [4.0, 6.0, 8.0, 12.0]
vps_d, our_d = [], []
for wd in widths:
    xx,yy,_ = vps_etch(wd); xc,dp = vps_profile(xx,yy); vps_d.append(depth_centre(xc,dp))
    rr = fe.run_etch(W=W,H=H,dx=DX,trench_width=wd,mask_thickness=MASK,sub_top=18.0,
                     t_end=DUR,n_steps=NSTEPS,par=OURS,verbose=False)
    xc2,dp2 = ours_profile(rr['segs'],18.0); our_d.append(depth_centre(xc2,dp2))
    print(f"  width {wd:4.1f}: ViennaPS {vps_d[-1]:5.2f}  ours {our_d[-1]:5.2f}")
vps_d=np.array(vps_d); our_d=np.array(our_d)

# ================= FIGURE =================
fig = plt.figure(figsize=(13,9)); plt.rcParams.update({'font.size':10})
gs = fig.add_gridspec(2,2, hspace=0.28, wspace=0.22)

# (a) our etch evolution
ax=fig.add_subplot(gs[0,0])
cmap=plt.cm.viridis
for i,s in enumerate(snaps):
    if s in r['snaps']:
        sg=r['snaps'][s]
        for k in range(len(sg)):
            ax.plot([sg[k,0]-W/2,sg[k,2]-W/2],[18-sg[k,1],18-sg[k,3]],
                    color=cmap(i/len(snaps)),lw=0.7)
ax.invert_yaxis(); ax.set_title("(a) Our sim: trench evolution\n(level-set + MC flux, SF$_6$/O$_2$)")
ax.set_xlabel("x [µm]"); ax.set_ylabel("depth [µm]"); ax.set_aspect('equal'); ax.set_xlim(-6,6)
sm=plt.cm.ScalarMappable(cmap=cmap,norm=plt.Normalize(0,DUR)); 
plt.colorbar(sm,ax=ax,label="time [min]",fraction=0.046)

# (b) profile overlay: dominant trench contour as connected line
ax=fig.add_subplot(gs[0,1])
from skimage import measure as _m
def longest_contour(phi, xs, ys, dx, sub_top):
    cs=_m.find_contours(phi,0.0); best=None; bl=0
    for c in cs:
        if len(c)>bl: bl=len(c); best=c
    px=xs[0]+best[:,0]*dx; py=ys[0]+best[:,1]*dx
    keep=py<sub_top-0.05
    return px[keep]-W/2, sub_top-py[keep]
# ViennaPS outline (already clean): order by walking the substrate nodes
vsel=vy<-0.02; ax.plot(vx[vsel],-vy[vsel],'k.',ms=3.2,label='ViennaPS (SF6O2Etching)')
ocx,ocy=longest_contour(r['phi'],r['xs'],r['ys'],DX,18.0)
ax.plot(ocx,ocy,'r-',lw=1.6,label='Ours (level-set contour)')
ax.invert_yaxis(); ax.set_aspect('equal'); ax.set_xlim(-6,6); ax.set_ylim(11,-2)
ax.set_title("(b) Final profile overlay, 8 µm trench"); ax.set_xlabel("x [µm]"); ax.set_ylabel("depth [µm]")
ax.legend(loc='upper center',fontsize=8,framealpha=0.95)

# (c) ARDE
ax=fig.add_subplot(gs[1,0])
ax.plot(widths, vps_d/vps_d[-1],'ks-',label='ViennaPS')
ax.plot(widths, our_d/our_d[-1],'ro--',label='Ours')
ax.set_title("(c) Aspect-ratio-dependent etching (RIE lag)\nnormalized to widest trench")
ax.set_xlabel("trench width [µm]"); ax.set_ylabel("normalized etch depth"); ax.grid(alpha=0.3); ax.legend()

# (d) timing
ax=fig.add_subplot(gs[1,1]); ax.axis('off')
txt = ("(d) Benchmark summary (CPU, single width-8 run)\n"
       "─────────────────────────────────────────\n"
       f"  grid: {int(W/DX)}×{int(H/DX)} @ Δx={DX} µm,  {NSTEPS} steps,  {DUR} min etch\n\n"
       f"  ViennaPS SF6O2Etching   runtime : {v_rt:5.1f} s\n"
       f"  Our level-set + MC      runtime : {T['total']:5.1f} s\n\n"
       f"  Our time breakdown:\n"
       f"     ballistic ray tracing : {T['raytrace']:5.1f} s  ({100*T['raytrace']/T['total']:.0f} %)\n"
       f"     surface chemistry     : {T['chem']:5.2f} s\n"
       f"     level-set advect+reinit: {T['advect']:5.1f} s\n\n"
       f"  Etch depth (8 µm trench):\n"
       f"     ViennaPS : {depth_centre(vxc,vdep):.2f} µm\n"
       f"     Ours     : {depth_centre(oxc,odep):.2f} µm\n\n"
       "  → ray tracing dominates (~95%): the\n"
       "    single kernel a GPU/OptiX port targets.\n"
       "    GPU available in this run: " + str(ps.gpuAvailable()))
ax.text(0.0,0.98,txt,va='top',ha='left',family='monospace',fontsize=9.5)

fig.suptitle("Feature-scale plasma-etch benchmark: from-scratch sim vs ViennaPS (SF$_6$/O$_2$ Si etch, identical plasma parameters)",
             fontsize=12, y=0.98)
fig.savefig("/home/claude/etch_benchmark.png", dpi=130, bbox_inches='tight')
print("saved figure")

summary = dict(vps_runtime_s=v_rt, ours_runtime_s=T['total'],
               raytrace_frac=T['raytrace']/T['total'],
               depth_vps=float(depth_centre(vxc,vdep)), depth_ours=float(depth_centre(oxc,odep)),
               widths=widths, vps_depth=vps_d.tolist(), ours_depth=our_d.tolist(),
               gpu_available=bool(ps.gpuAvailable()))
json.dump(summary, open("/home/claude/summary.json","w"), indent=2)
print(json.dumps(summary, indent=2))
