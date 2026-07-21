# Research: AMR + deep-feature (200:1) strategy for petch — literature and recommendation (2026-07-21)

Scope: how to make the uniform-Eulerian-grid petch engine (level-set phi, dx=10nm authority,
conservative surface-state remap on marching-cubes triangles, deterministic exact 2-D diffuse
exchange, GPU ballistic transport) tractable and accurate for 200:1 holes/trenches
(40–90 nm opening x 8–18 um depth). Sources cited inline. NOT committed; working research doc.

---

## 0. Executive answer — what to build first

The literature is unambiguous on one point: **nobody solves 200:1 by brute-forcing a fine
uniform 3-D grid + all-pairs exchange; every credible approach splits the feature into a
resolved region (mouth + etch front) and a reduced "tube" region (1-D radiosity / Clausing
conductance), and/or grades the mesh along the tube axis.** Full octree AMR is the *last*
thing to build, not the first — for a straight tube it buys less than an axis-graded grid
plus a 1-D transport reduction, and it costs the most in reinit/remap/GPU complexity.

Build order (details and costs in §5):

1. **Stage 0 (days): axisymmetric hole mode + Clausing gates.** Ring-radiosity exact
   exchange (coaxial-disk view factors = the axisymmetric twin of today's crossed-strings
   2-D operator). Free validation: Santeler/Clausing transmission values (1.30% at 100:1,
   0.66% at 200:1). Aligns with the pinned partner deliverable ("axisymmetric hole mode =
   biggest lever for 200:1 dynamics", NEXT_STEPS.md §3).
2. **Stage 1 (1–2 wk): axis-graded (tall-cell) grid.** Fine lateral dx (2.5–5 nm — the
   10 nm authority *cannot* resolve a 40 nm opening at the >=12 cells/opening gate), coarse
   axial dz (10–20 nm) in the straight tube. Precedent: graphics "tall cell" free-surface
   grids (Irving et al. 2006; Chentanez & Müller 2011).
3. **Stage 2 (2–4 wk): domain decomposition mouth/tube/front with matched BCs.** Full exact
   exchange in the mouth region and within ~5–10 CD of the front; 1-D radiosity
   (Kokkoris-style, view factors from coaxial disks / crossed strings) in the straight
   middle; fail-closed **convexity monitor** flips a tube slab back to the full solver when
   bowing/necking appears. This is exactly the published Kokkoris/Manstetten architecture.
4. **Stage 3 (as needed for 3-D/LER): hierarchical radiosity** (Hanrahan-style O(N) block
   matrix) over the triangle exchange — mandatory before any full-3-D 200:1 diffuse solve
   (all-pairs at 2.5 nm is ~10^10 pairs; see §1).
5. **Stage 4 (only if bowed/non-straight deep geometry at fine resolution forces it):
   two-level nested/block-structured narrow-band refinement** (Lenz/ViennaLS-style, or
   AMReX-style blocks), *not* a general octree. Keep the interface band single-resolution;
   coarsen only the far field — this sidesteps the known cross-level reinit/mass-loss
   failure modes (§2).

---

## 1. Cost anatomy at 200:1 (what actually binds)

Geometry: opening W = 40–90 nm, depth D = 8–18 um, AR 200. Take W=40 nm, D=8 um.

- **The 10 nm authority under-resolves the cross-section.** 40 nm opening = 4 cells. The
  repo's own ARDE gate (tests/test_arde_transport.py lineage) required >=~12 cells across
  the opening; the R4 quadrature lesson says resolution/sampling is physics-grade. So
  **lateral refinement to 2.5–5 nm is mandatory at 200:1 regardless of any AMR decision.**
  This is the single most load-bearing fact: the problem is not "coarsen the tube" first,
  it is "refine the cross-section without paying for it along 8 um of axis".
- **phi storage is NOT the bottleneck.** Dense 2-D at 2.5 nm, one 120 nm-pitch cell:
  ~48 x 3200 = 1.5e5 cells. Dense 3-D: 48 x 48 x 3200 ~ 7.4e6 cells (~30 MB/field) —
  fits on any GPU; a narrow band (surface cells x band width ~ 1e6 nodes) is MBs. Sparse
  level-set structures (HRLE/VDB, §2–3) make even multi-hole domains cheap.
- **Diffuse-exchange pair count IS the bottleneck.**
  - 2-D trench, wall segments at 2.5 nm: N ~ 2 x 3200 = 6400 segments -> ~4e7 ordered pairs
    per rebuild. Feasible but wasteful; with axial grading (segments ~ dz = 15 nm on straight
    walls) N ~ 1100 -> ~1e6 pairs — trivial.
  - 3-D hole, triangles at 2.5 nm: N ~ (pi*40/2.5) x (8000/2.5) ~ 1.6e5 triangles ->
    ~2.6e10 pairs. **All-pairs is dead in 3-D at 200:1.** Hierarchical radiosity (O(N) blocks,
    Hanrahan et al. 1991) or the 1-D tube reduction is not an optimization, it is an
    enabling requirement.
- **Ballistic (GPU ray) transport** scales with rays x bounce depth; a 200:1 tube multiplies
  wall-bounce path lengths ~AR. The tube reduction (§4) removes both the diffuse N^2 and the
  long ballistic corridors from the resolved region.

---

## 2. AMR for level-set interface evolution — what is proven, what breaks

**Octree narrow-band level sets (graphics lineage).**
- Losasso, Gibou & Fedkiw, "Simulating water and smoke with an octree data structure",
  SIGGRAPH 2004 ([ACM](https://dl.acm.org/doi/10.1145/1015706.1015745)): unrestricted octree,
  semi-Lagrangian advection of velocity and level set, symmetric-positive-definite octree
  Poisson discretization. Known weakness: semi-Lagrangian advection is **non-conservative**
  (volume loss); graphics papers of that era paired it with particle level sets to patch mass
  loss — acceptable for film, not for a metrology-grade etch ledger.
- Min & Gibou, "A second order accurate level set method on non-graded adaptive Cartesian
  grids", JCP 225 (2007) 300–321
  ([PDF](https://www.ljll.fr/~frey/papers/levelsets/Min%20C.,%20Gibou%20F.,%20A%20second%20order%20accurate%20level%20set%20method%20on%20non%20graded%20adaptive%20cartesian%20grids.pdf)):
  the reference for *sharp-interface* quality on quad/octrees — second-order accurate
  advection on **non-graded** trees (no 2:1 balance requirement), locally third-order
  reinitialization, "negligible mass loss" demonstrated. This is the octree scheme to copy
  if Stage 4 ever becomes a full tree.
- Mirzadeh, Guittet, Burstedde & Gibou, "Parallel level-set methods on adaptive tree-based
  grids", JCP 322 (2016) 345–364
  ([ScienceDirect](https://www.sciencedirect.com/science/article/abs/pii/S002199911630242X)):
  the same machinery parallelized on **p4est** forests-of-octrees (MPI, CPU); demonstrates
  scalable tree regrid + semi-Lagrangian LS + reinit. p4est itself: Burstedde, Wilcox &
  Ghattas, SIAM J. Sci. Comput. 33 (2011)
  ([PDF](https://ins.uni-bonn.de/media/public/publication-media/BursteddeWilcoxGhattas11.pdf?pk=583)).
  Note: this ecosystem is MPI/CPU-proven, **not** GPU-proven for the level-set part.
- Reinitialization across levels: the eikonal solvers that are GPU-proven are
  uniform-grid — Detrixhe, Gibou & Min, "A parallel fast sweeping method for the Eikonal
  equation", JCP 237 (2013) 46–55
  ([ScienceDirect](https://www.sciencedirect.com/science/article/abs/pii/S002199911200722X)).
  On trees, reinit stencils that straddle coarse-fine boundaries perturb the interface at
  O(dx_coarse) unless the whole interface band sits at the finest level — which is exactly
  the standard mitigation: **keep a uniform finest-band around the zero level set; coarsen
  only the far field.** For petch this is cheap because the far field of a 200:1 tube is
  precisely the featureless region.

**Sparse voxel structures (VDB lineage).**
- Museth, "VDB: High-Resolution Sparse Volumes with Dynamic Topology", ACM TOG 32(3) 2013
  ([PDF](https://www.museth.org/Ken/Publications_files/Museth_TOG13.pdf)): the production
  standard for narrow-band level sets with dynamic topology — O(1) random access, CPU.
  Crucially VDB is **sparse but single-resolution at the leaf level** — it solves *memory*,
  not *resolution grading*.
- NanoVDB (Museth, SIGGRAPH 2021 talk;
  [overview](https://history.siggraph.org/learning/nanovdb-a-gpu-friendly-and-portable-vdb-data-structure-for-real-time-rendering-and-simulation-by-museth/),
  [ASWF note](https://www.digitalmediaworld.tv/vfx/3110-academy-software-foundation-openvdb-adds-gpu-support-with-nanovdb)):
  GPU-portable **linearized, topology-static snapshot** of VDB. Good for GPU transport
  queries against a frozen surface; **not** a GPU dynamic-topology level-set solver. What
  breaks on GPU: narrow-band rebuild (topology mutation) still round-trips to CPU/OpenVDB.
- Verdict: VDB/NanoVDB-style sparsity maps well onto petch's existing narrow-band +
  GPU-mesh pipeline; a full GPU octree LS solver has no proven off-the-shelf exemplar.

**Block-structured AMR frameworks.**
- AMReX: Zhang et al., IJHPCA 2021
  ([SAGE](https://journals.sagepub.com/doi/10.1177/10943420211022811)); GPU-portable
  (CUDA/HIP/SYCL) block-structured AMR with embedded-boundary (cut-cell) geometry, and —
  the key discipline petch would have to import — **flux registers to restore conservation
  at coarse-fine interfaces** ("the area- and time-weighted fluxes from level l and l+1 faces
  do not agree at the interface, resulting in a loss of conservation… remedied by the
  EBFluxRegister refluxing" — [AMReX docs](https://amrex-codes.github.io/amrex/docs_html/AmrCore.html),
  [EB docs](https://amrex-codes.github.io/amrex/docs_html/EB.html)). AMReX is the proof that
  **2-level block refinement on GPU is mature engineering**, while octree-node schemes on GPU
  are research.

**What breaks — summary for our contract.**
1. Semi-Lagrangian LS advection on trees: non-conservative; volume drift unless
   Min-Gibou-quality interpolation + reinit (they report negligible loss; verify per-step
   with our ledger, fail-closed).
2. Reinit across refinement levels: interface-adjacent stencils must never cross a level
   jump — enforce "interface band = finest level" as an invariant.
3. Conservation at coarse-fine boundaries for any flux-form quantity (our surface-state
   ledgers): needs explicit flux-register-style synchronization (AMReX pattern).
4. GPU: block-structured 2-level = proven (AMReX); dynamic octree LS = unproven; VDB
   topology updates = CPU.

---

## 3. What etch/deposition simulators actually do

- **ViennaLS / ViennaPS (TU Wien, open source).** Sparse-field level set on a
  **hierarchical run-length-encoded (HRLE) data structure** — sparse in memory, but a
  *single* grid resolution (Ertl's parallel HRLE lineage; ViennaLS
  [GitHub](https://github.com/ViennaTools/ViennaLS); ViennaPS SoftwareX 2025
  ([ScienceDirect](https://www.sciencedirect.com/science/article/pii/S2352711025004194)).
  Flux by Monte Carlo ray tracing (Embree CPU; OptiX GPU since v3.4.0, "experimental",
  disk/line/level-set-voxel primitives, **9–34x GPU speedup** — Riedel diploma thesis 2026,
  [reposiTUm](https://repositum.tuwien.at/handle/20.500.12708/226294)). So the closest
  open-source competitor is: sparse-but-uniform LS + GPU MC rays — i.e., petch's current
  architecture class. They do NOT have deep-tube reductions in the main line.
- **TU Wien hierarchical-grid research line (the published "feature-scale AMR etch" work).**
  Lenz, "Hierarchical Grid Algorithms for Topography Simulation", PhD TU Wien 2023
  ([reposiTUm](https://repositum.tuwien.at/handle/20.500.12708/188148)): curvature-based
  feature detection on the zero level set drives **local sub-grid refinement** ("detected
  features … used to … locally adapt the resolution of the hierarchical grid"); reported
  gains are modest — up to **58% runtime improvement** on SiGe selective-epi, 15% on MC flux
  via simplified meshes. Companion paper: "Curvature Based Feature Detection for Hierarchical
  Grid Refinement in TCAD Topography Simulations", Solid-State Electronics 2022
  ([ScienceDirect](https://www.sciencedirect.com/science/article/pii/S0038110122000302),
  [IEEE SISPAD version](https://ieeexplore.ieee.org/document/9560690)). Takeaway: published
  feature-scale LS-AMR exists, is two-level/nested (not octree), and yields <2x — supporting
  the judgment that grid AMR alone is not where the 200:1 win is.
- **Synopsys Sentaurus Topography (2-D/3-D).** Commercial; level-set surface evolution on a
  uniform rectilinear grid; documentation explicitly frames resolution as a global
  accuracy/memory tradeoff (no public AMR)
  ([Synopsys page](https://www.synopsys.com/manufacturing/tcad/process-simulation/sentaurus-topography.html),
  [tutorial](https://ghzphy.github.io/Sentaurus_Training/stopo/stopo_1.html)).
- **Coventor SEMulator3D (now Lam Research).** **Uniform voxel** "virtual fabrication"
  engine (voxel + surface-evolution hybrid); scales by throwing memory/compute at uniform
  voxels, is process-emulation (calibrated behavioral) rather than transport-physics
  ([Lam page](https://www.lamresearch.com/product/semulator3d/),
  [SemiWiki overview](https://semiwiki.com/x-subscriber/coventor/2430-semulator3d-a-virtual-fab-platform/)).
- **Kushner's MCFPM (Michigan).** **Uniform rectilinear cell-based mesh** (2-D and 3-D),
  material-identity cells, MC pseudoparticles; resolution pushed very fine on small domains
  (their HAR twisting studies quote sub-nm/atomistic-scale cells)
  ([MCFPM page](https://cpseg.eecs.umich.edu/Projects/MCFPM/MCFPM.htm),
  [HAR twisting project](https://cpseg.eecs.umich.edu/Projects/TWISTING/Twisting_Ebm_v02.html)).
  No AMR; they pay uniform-grid cost and accept stochastic noise (the exact axis petch
  differentiates on).
- **HAR-specific simulators.** ViPER (Kokkoris et al.): HAR silicon etch profile simulator
  built around the 1-D flux framework
  ([ResearchGate](https://www.researchgate.net/publication/258159988_ViPER_Simulation_software_for_high_aspect_ratio_plasma_etching_of_silicon)).
  Lam-affiliated analyses of HAR transport are quasi-1-D Clausing-type (§4).

**Conclusion:** no production etch simulator runs interface-adapted octrees. The field's
answer to deep features is (a) uniform + sparse + GPU rays (Vienna), (b) uniform voxel +
money (Coventor/MCFPM), or (c) quasi-1-D transport reduction (Kokkoris/ViPER, Lam). A
correct-by-construction deterministic exchange + tube reduction + graded grid would be
differentiated against all three.

---

## 4. The 1-D-tube exploitation: conductance, slice reduction, validity

**Exact/reference transmission probabilities (Clausing factors), circular tube.**
- Clausing 1932 (Ann. Physik 12, 961; reprinted JVST 8 (1971) 636) posed the integral
  equation; modern reference values: Cole's upper/lower bounding sequences for l/r from 0.1
  to 1000 (agreeing to many digits; "published Clausing factors… differ from one another in
  the 5th or even 3rd significant figure for long tubes; Cole appears most precise")
  ([Vacuum 2002 comparison](https://www.sciencedirect.com/science/article/abs/pii/S0042207X02002257);
  Cole 1977, Prog. Astronaut. Aeronaut. 51, 261
  ([ADS](https://ui.adsabs.harvard.edu/abs/1977PrAA...51..261C/abstract))).
- Practical closed forms (CERN "Vacuum Technology for Ion Sources", arXiv:1404.0960
  ([PDF](https://arxiv.org/pdf/1404.0960)):
  - **Santeler** (JVST A 4 (1986) 348), error <0.7% at all L/R:
    tau = 1 / (1 + (3L/8R) * (1 + 1/(3(1 + L/(7R))))).
  - Long-tube limit: tau -> 8R/(3L) = **4/(3*AR)** with AR = depth/diameter.
  - Dushman interpolation: tau ~ 1/(1 + 3L/8R).
  - Berman series: 0.1% accuracy regimes for short (l/r<0.5) and long (l/r>20) tubes
    ([Vacuum 1978 approx. calc](https://www.sciencedirect.com/science/article/abs/pii/S0042207X78802437)).
- Anchor numbers for gates (Santeler): AR=10 -> 10.9%; AR=50 -> 2.54%; **AR=100 -> 1.30%**
  (matches the industry figure cited in NEXT_STEPS.md and by Panagopoulos & Lill); 
  **AR=200 -> 0.66%**. Slots/trenches: transmission falls only ~ (W/L)ln(L/W) (2-D slit,
  DeMarcus/Berman lineage; rectangular/slit tables in
  [Vacuum 2013 polygon-tube study](https://www.sciencedirect.com/science/article/abs/pii/S0042207X12004964)) —
  trenches starve much more slowly than holes; the hole is the hard case.
- Direct etch-context authority: Panagopoulos & Lill, "Neutral transport during etching of
  high aspect ratio features", JVST A 41, 033006 (2023)
  ([AIP](https://pubs.aip.org/avs/jva/article/41/3/033006/2877892)): neutral flux at the
  bottom of a 100:1 hole = **1.3%** of mouth flux; validates MC transport against Clausing
  analytics as a function of AR, profile shape, and sticking. Also: Coburn & Winters'
  conductance ARDE model and Gottscho's neutral-ion synergy extension (summarized in the
  ARDE literature, e.g. [JVST A 35, 05C301 (2017)](https://pubs.aip.org/avs/jva/article/35/5/05C301/244889);
  multiscale confirmation [PSST 2023](https://iopscience.iop.org/article/10.1088/1361-6595/acdc4f)).

**Slice-marching / quasi-1-D coupled models (the published version of "1b").**
- Kokkoris, Boudouvis & Gogolides, "Integrated framework for the flux calculation of
  neutral species inside trenches and holes during plasma etching" (JVST A 2006) and the
  TU Wien continuation: Manstetten et al., "Framework to model neutral particle flux in
  convex high aspect ratio structures using one-dimensional radiosity", Solid-State
  Electronics 2016 ([ScienceDirect](https://www.sciencedirect.com/science/article/abs/pii/S0038110116301964));
  Manstetten PhD, "Efficient Flux Calculations for Topography Simulation", TU Wien 2018
  ([thesis](https://www.iue.tuwien.ac.at/phd/manstetten/)). Mechanics: discretize the
  feature into axial rings/slices; **view factors between coaxial disks** (hole) or via
  **crossed strings** (trench — literally today's petch 2-D operator); constant flux +
  sticking per element; reformulated "receiving" radiosity handles fully absorbing elements;
  **agrees with 3-D MC ray tracing across ARs for convex geometries** at negligible cost;
  Manstetten reports overall topography-simulation speedups >14x from the combined
  acceleration set. This is precisely the proposed petch tube segment, already validated in
  the literature.
- ALD conformality lineage (independent confirmation of quasi-1-D validity in tubes):
  Ylilammi et al. diffusion-reaction model for lateral HAR channels; Thiele-modulus
  analysis; extended across Kn regimes
  ([PCCP 2024](https://pubs.rsc.org/en/Content/ArticleLanding/2024/CP/D4CP00131A);
  [ballistic-vs-diffusion comparison, Aalto](https://aaltodoc.aalto.fi/items/c9037e60-d5fe-4020-aa1a-d70594729f9e)).
  Useful as the *cheapest* tube closure (1-D reaction-diffusion with Knudsen diffusivity)
  below the fidelity of 1-D radiosity.

**Hierarchical radiosity (for when the resolved region is 3-D and large).**
- Hanrahan, Salzman & Aupperle, "A rapid hierarchical radiosity algorithm", SIGGRAPH 1991
  ([ACM](https://dl.acm.org/doi/abs/10.1145/122718.122740)): adaptive block decomposition of
  the form-factor matrix to a user-supplied error bound; **O(n) interaction blocks** after
  an O(k^2) top level; uniform-precision form factors. The natural upgrade path for petch's
  exact exchange: keep exact near-field pairs, cluster far-field tube-wall panels — on a
  near-cylinder the far field is smooth, so clustering error is small and certifiable
  (error-bounded refinement literature:
  [error-bounding algorithms for HR](https://inria.hal.science/inria-00098630/file/HS98.pdf)).

**Validity bounds of the tube reduction (when it breaks).**
1. **Convexity is the hard requirement** (stated by both Kokkoris and Manstetten): 1-D
   radiosity assumes every element sees every other unoccluded. **Bowing** (barrel), 
   **necking/overhang**, and twist create self-shadowing -> the slab must fall back to the
   full 2-D/3-D solver. Practical monitor: per-slab convexity of the cross-section radius
   profile r(z) (sign of d2r/dz2 crossings + overhang detection on the MC triangles);
   fail-closed to the resolved solver — same design language as the extrusion certification.
2. **Tapering is fine** — conical-tube transmission probabilities are classical (Cole 1977;
   [MC vs analytic for conical tubes](https://www.researchgate.net/publication/242364926));
   slow taper keeps view-factor positivity; the 1-D radiosity handles r(z) variation
   natively via coaxial-disk factors.
3. **Angular resolution at interfaces.** A conductance scalar (Clausing factor) assumes a
   diffuse (cosine) re-emitted distribution; deep-tube *ballistic* species develop strongly
   forward-peaked angular distributions. The matched BC between resolved region and tube
   segment must pass an angular decomposition (at minimum: direct/collimated + diffuse
   partition; petch already separates ballistic and diffuse operators, which is exactly the
   right split — Clausing/1-D radiosity applies to the diffuse channel only).
4. **Re-emission physics richer than constant sticking** (charging-dependent sticking,
   polymer ledgers) works unchanged: 1-D radiosity supports per-ring sticking/sources; only
   *non-local* view-factor changes (shape) break it, not surface-state variation.

---

## 5. Anisotropic/graded axial meshes + the conservative remap contract

- **Graded/anisotropic grids for level sets are routine when grading is smooth and the
  fine direction resolves the interface curvature** (anisotropic-mesh level-set literature:
  [topology-optimization on anisotropic graded meshes](https://www.sciencedirect.com/science/article/pii/S0096300323000723);
  mesh-adapted LS reviews, e.g. [Dapogny/Allaire shape-opt LS remeshing](https://dapogny.org/publis/shapeopfinal.pdf)).
  Caveats: (a) reinit on stretched cells must use per-axis metric (signed distance in
  *physical* space, Godunov Hamiltonian with dx != dz — mechanical change to the existing
  uniform reinit); (b) upwind stencils see anisotropic truncation error — acceptable when
  the surface is near-parallel to the coarse axis (straight walls: normal is lateral, so the
  *fine* axis resolves the normal; exactly our case); (c) never let the etch FRONT live in
  tall cells — cap dz near the front (the front region stays isotropic-fine).
- **Direct precedent for "long featureless column = tall cells":** Irving, Guendelman,
  Losasso & Fedkiw, "Efficient simulation of large bodies of water by coupling two and three
  dimensional techniques", SIGGRAPH 2006
  ([ACM](https://dl.acm.org/doi/10.1145/1141911.1141959)) — 3-D fine cells near the free
  surface, tall thin cells below; and Chentanez & Müller, "Real-time Eulerian water
  simulation using a restricted tall cell grid", SIGGRAPH 2011
  ([ACM](https://dl.acm.org/doi/10.1145/1964921.1964977)) — the "restricted" two-layer
  version (uniform fine band + tall cells) that is the cleanest template: petch's tube walls
  are the water column, the etch front is the free surface.
- **Conservative surface-state remap across resolution changes.** The proven practice is
  intersection-based ("supermesh"/common-refinement) transfer: build the intersection of old
  and new discretizations and transfer by exact area-weighted integrals (Galerkin
  projection) — Farrell & Maddison, "Conservative interpolation between volume meshes by
  local Galerkin projection", CMAME 2011
  ([PDF](https://www.ljll.fr/~frey/papers/divers/Farrell%20P.E.,%20Conservative%20interpolation%20between%20volume%20meshes%20by%20local%20Galerkin%20projection.pdf);
  supermesh construction [paper](https://www.researchgate.net/publication/223420260_Conservative_interpolation_between_unstructured_meshes_via_supermesh_construction);
  FV extension [CMAME 2011](https://www.sciencedirect.com/science/article/abs/pii/S0045782511001666)).
  For petch's marching-cubes-triangle-resident state: project old triangles onto new ones
  along the surface normal within the band, clip (2-D polygon intersection in the local
  tangent frame), transfer state by intersection area; ledger total adsorbate/polymer before
  and after; fail-closed on ledger drift. This is the same contract as today's step-to-step
  remap — a resolution change is just a bigger remap, and supermesh clipping is the
  conservative way to do it. AMReX's flux registers (§2) are the volume-field analog if
  cell-resident fields ever cross levels.

---

## 6. Recommended architecture for petch, staged, with honest costs

**Stage 0 — Axisymmetric hole mode + Clausing gates (days; do first).**
Ring-discretized exact exchange with analytic coaxial-disk view factors (the axisymmetric
twin of crossed-strings; formulas in the view-factor catalogs used by Manstetten). Cost:
N_rings ~ depth/dz ~ 500–1500; N^2 ~ 1e6 — trivial, deterministic, zero sampling noise.
Unlocks: 200:1 static starvation ladder (partner deliverable) with *analytic* gates
(Santeler curve; 1.30% @100:1, 0.66% @200:1); today's extruded-2-D operator stays the
authority for trenches, ring mode becomes the authority for holes. NOTE: today's
extruded-2-D exact exchange *cannot* represent a hole — extrusion gives a trench; at 200:1
the hole/trench difference is the whole game (4/(3AR) vs ~(W/L)lnL/W starvation).

**Stage 1 — Axis-graded (tall-cell) grid (1–2 weeks).**
Lateral dx 2.5–5 nm (>=12 cells/opening gate), axial dz graded 5 nm (front + mouth bands)
-> 15–20 nm (straight tube). Implementation: per-axis metric in advection/reinit + graded-dz
marching-cubes; state remap across regrade via supermesh clipping (§5). Cost at 200:1 2-D:
~1e5–2e5 cells, wall segments ~1–2e3, exchange pairs ~1e6 — the *entire 2-D 200:1 problem
becomes interactive*. Risk: low (no level jumps, single structured grid, smooth grading).

**Stage 2 — Domain decomposition mouth/tube/front with matched BCs (2–4 weeks).**
Three segments: (i) mouth (isotropic-fine, full exchange + ballistic rays — reflects the
above-wafer view factor and mask effects), (ii) tube (1-D radiosity rings/slices, Kokkoris/
Manstetten; diffuse channel only; ballistic species propagated analytically through the
straight segment with angular binning), (iii) front (isotropic-fine full solver, extends
5–10 CD above the instantaneous bottom and *marches down* with it — "slice marching").
Matched BCs: angular-partitioned flux (collimated + cosine-diffuse at minimum) passed at
both interfaces; reciprocity ledger across the seams. Fail-closed convexity monitor per tube
slab (§4) reverts slabs to resolved mode. Cost: the resolved regions are AR<~10 problems
regardless of total AR — **total cost becomes ~independent of depth** (the 1-D segment is
O(N_rings) per solve). This stage is what makes 200:1 *dynamics* (not just statics) cheap.
Risk: medium — the seams are new certified surfaces; validation ladder §7 covers them.

**Stage 3 — Hierarchical radiosity for resolved-3-D (when 3-D/LER at depth is needed).**
Hanrahan-style block clustering over the triangle exchange with a declared error bound;
exact pairs near-field, clustered far-field. Required before any full-3-D 200:1 diffuse
solve (2.6e10 pairs otherwise, §1). Fits the "adaptive rays #1" NEXT_STEPS item as its
deterministic-exchange counterpart.

**Stage 4 — Two-level nested narrow-band refinement (only on demand).**
If deep *non-straight* geometry (bowed 3-D profiles, twist) must be resolved at 2.5 nm over
microns, adopt the restricted two-level pattern: uniform-finest band at the interface,
coarse far field, block-structured (AMReX-style flux-register discipline; Lenz-style
curvature-triggered placement; Min-Gibou numerics if it ever becomes a true tree). Do NOT
build a general octree first — no etch-industry precedent (§3), GPU path unproven (§2), and
Stages 0–3 remove the tube cost that motivates it.

---

## 7. Validation ladder (each stage gates on the previous)

1. **Clausing-limit gates (static).** Straight cylinder, unit sticking at bottom, diffuse
   walls: bottom flux vs AR must track Santeler/Cole (<=1% dev; anchors 10.9%@10, 2.54%@50,
   1.30%@100, 0.66%@200). Trench twin vs 2-D slit values. Run vs AR sweep, both ring mode
   and (for trench) today's extruded operator: this is also the ring-mode acceptance test.
2. **Manufactured tube solutions.** Uniform-sticking cylinder radiosity has semi-analytic
   axial flux profiles (Kokkoris/Manstetten reproduce 3-D MC); manufacture r(z) tapers and
   verify 1-D radiosity vs the full exact-exchange solver on the same geometry. Reciprocity
   + row-sum (escape) ledgers on every operator, per step (existing petch discipline).
3. **Seam consistency (Stage 2).** Shallow-vs-deep: for AR 5–20 run (a) full resolved solver,
   (b) decomposed mouth/tube/front; profiles and fluxes must agree within declared tol
   everywhere, *including transient* (short dynamic etches). Then hold (b) fixed and extend
   to AR 50/100/200 where (a) is unaffordable — the physics extrapolation is now carried by
   gates 1–2.
4. **Grading consistency (Stage 1).** Same physics on uniform-fine vs graded grid at AR<=20:
   profile Hausdorff distance within tol; state-ledger drift zero to roundoff across
   regrades; grid ladder in the *lateral* axis (>=12 cells/opening; 10/5/2.5 nm).
5. **Convexity-monitor fault injection.** Manufactured bowed/necked profiles must trip the
   monitor and revert slabs; verify the reverted solve matches full-resolved reference.
6. **End-to-end anchor.** de Boer deep-trench series (existing) + the 1.3%@100:1 industry
   figure (Panagopoulos-Lill) as the hole-mode anchor; later Nozawa/charging protocols per
   NEXT_STEPS.

---

## 8. Key sources

AMR/level set: [Losasso-Gibou-Fedkiw 2004](https://dl.acm.org/doi/10.1145/1015706.1015745) ·
[Min-Gibou 2007](https://www.ljll.fr/~frey/papers/levelsets/Min%20C.,%20Gibou%20F.,%20A%20second%20order%20accurate%20level%20set%20method%20on%20non%20graded%20adaptive%20cartesian%20grids.pdf) ·
[Mirzadeh et al. 2016 (p4est)](https://www.sciencedirect.com/science/article/abs/pii/S002199911630242X) ·
[Burstedde-Wilcox-Ghattas 2011](https://ins.uni-bonn.de/media/public/publication-media/BursteddeWilcoxGhattas11.pdf?pk=583) ·
[Detrixhe-Gibou-Min 2013](https://www.sciencedirect.com/science/article/abs/pii/S002199911200722X) ·
[Museth VDB 2013](https://www.museth.org/Ken/Publications_files/Museth_TOG13.pdf) ·
[NanoVDB](https://history.siggraph.org/learning/nanovdb-a-gpu-friendly-and-portable-vdb-data-structure-for-real-time-rendering-and-simulation-by-museth/) ·
[AMReX 2021](https://journals.sagepub.com/doi/10.1177/10943420211022811) ·
[AMReX flux registers](https://amrex-codes.github.io/amrex/docs_html/AmrCore.html)

Simulators: [ViennaPS SoftwareX 2025](https://www.sciencedirect.com/science/article/pii/S2352711025004194) ·
[ViennaLS](https://github.com/ViennaTools/ViennaLS) ·
[ViennaPS GPU ray tracing thesis 2026](https://repositum.tuwien.at/handle/20.500.12708/226294) ·
[Lenz hierarchical-grid PhD 2023](https://repositum.tuwien.at/handle/20.500.12708/188148) ·
[Curvature-refinement SSE 2022](https://www.sciencedirect.com/science/article/pii/S0038110122000302) ·
[Sentaurus Topography](https://www.synopsys.com/manufacturing/tcad/process-simulation/sentaurus-topography.html) ·
[SEMulator3D](https://www.lamresearch.com/product/semulator3d/) ·
[MCFPM](https://cpseg.eecs.umich.edu/Projects/MCFPM/MCFPM.htm) ·
[MCFPM HAR twisting](https://cpseg.eecs.umich.edu/Projects/TWISTING/Twisting_Ebm_v02.html)

Tube transport: [CERN vacuum notes (Santeler formula)](https://arxiv.org/pdf/1404.0960) ·
[Clausing-factor comparison, Vacuum 2002](https://www.sciencedirect.com/science/article/abs/pii/S0042207X02002257) ·
[Cole 1977](https://ui.adsabs.harvard.edu/abs/1977PrAA...51..261C/abstract) ·
[Panagopoulos-Lill JVSTA 2023](https://pubs.aip.org/avs/jva/article/41/3/033006/2877892) ·
[ARDE 3-D JVSTA 2017](https://pubs.aip.org/avs/jva/article/35/5/05C301/244889) ·
[Manstetten 1-D radiosity SSE 2016](https://www.sciencedirect.com/science/article/abs/pii/S0038110116301964) ·
[Manstetten PhD 2018](https://www.iue.tuwien.ac.at/phd/manstetten/) ·
[ViPER](https://www.researchgate.net/publication/258159988_ViPER_Simulation_software_for_high_aspect_ratio_plasma_etching_of_silicon) ·
[Ylilammi-lineage PCCP 2024](https://pubs.rsc.org/en/Content/ArticleLanding/2024/CP/D4CP00131A) ·
[Hanrahan et al. 1991](https://dl.acm.org/doi/abs/10.1145/122718.122740)

Meshes/remap: [Irving et al. 2006 tall cells](https://dl.acm.org/doi/10.1145/1141911.1141959) ·
[Chentanez-Müller 2011 restricted tall cells](https://dl.acm.org/doi/10.1145/1964921.1964977) ·
[Farrell-Maddison 2011 supermesh Galerkin](https://www.ljll.fr/~frey/papers/divers/Farrell%20P.E.,%20Conservative%20interpolation%20between%20volume%20meshes%20by%20local%20Galerkin%20projection.pdf) ·
[Supermesh construction](https://www.researchgate.net/publication/223420260_Conservative_interpolation_between_unstructured_meshes_via_supermesh_construction) ·
[Anisotropic graded LS meshes](https://www.sciencedirect.com/science/article/pii/S0096300323000723)
