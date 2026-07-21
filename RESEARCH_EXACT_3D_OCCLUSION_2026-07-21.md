# Exact-occlusion 3-D form factors: research report (Phase-2+ reference)

(Deep-research agent output, 2026-07-21. Companion to
FORM_FACTOR_IMPLEMENTATION_GUIDE_2026-07-21.md and to
`src/petch/deterministic_exchange_2d.py` — the 2-D `analytic_occlusion` operator whose 3-D
analogue this report designs. Sources verified by web search today; NISTIR 6925 read in full.)

## Verdict (architectural)

The 2-D operator generalizes along exactly one axis and it is the right one:
**exact per-source-point visible-region computation (blocker projection + polygon boolean +
closed-form Baum point-to-polygon per visible piece) as the inner integral, and a certified
2-D quadrature over the source triangle whose panels are split at analytically enumerated
visual events (EV lines exactly; EEE conics optionally) as the outer integral.**

This is *not* exotic: it is what NIST View3D already does for the inner part (shadow
projection + convex clipping + the exact Hottel angle formula), minus the event-split
certification — and Walton's own "Problem Case 2" (NISTIR 6925, Fig. 15) documents the exact
failure mode we found today in 2-D: a grazing sliver of visibility that **no integration
point of the outer quadrature ever sees, so adaptive subdivision never triggers and F is
silently computed as 0**. Sampled visibility (rays or Gauss-classified cells) has the same
disease with the same silence. The load-bearing asymmetry:

> With an EXACT inner occlusion operator, an incomplete event list can only make the outer
> quadrature converge slower (and the receipt says so). With SAMPLED visibility, a missed
> sliver biases the answer with no receipt. Exactness must live in the inner integral;
> certification must live in the outer one. Never put sampling in the visibility classifier.

Full global visibility structures (complete discontinuity mesh, visibility skeleton,
3-D visibility complex) are the theory source for the event enumeration but are **not worth
building** for meshes rebuilt every etch step (worst-case O(n^4) size; see §1, §3). Use
their event catalogue *locally, per pair*.

For 200:1 axisymmetric holes there is a reduction (§5) that collapses the whole 3-D exact
operator to the 2-D operator's cost structure, with occlusion events given in closed form
(arccos of a rational function of the profile-kink radii). That is the 80/20 and should be
built first.

---

## 1. Exact occlusion in 3-D radiosity: what the graphics literature actually provides

- **Nishita–Nakamae (SIGGRAPH '85; and TOG "Shading models for point and linear sources")**:
  first exact penumbra/umbra decomposition for polygonal area sources — shadow volumes per
  (source, blocker) pair, penumbra = convex-hull volumes, umbra = intersection volumes;
  illumination at a point computed by clipping the source polygon against blocker
  projections and applying the exact contour-integral (Lambert) formula. This is the
  per-point EXACT visible-polygon primitive we need, published 1985.
  [dl.acm.org/doi/10.1145/282918.282938](https://dl.acm.org/doi/10.1145/282918.282938),
  [history.siggraph.org](https://history.siggraph.org/learning/the-visibility-skeleton-a-powerful-and-efficient-multi-purpose-global-visibility-tool-by-durand-drettakis-and-puech/)
- **Heckbert 1992 (3rd Eurographics Workshop on Rendering)**: *incomplete* discontinuity
  meshing — builds only the EV wedges (source edge x blocker vertex, source vertex x blocker
  edge) and inserts them as mesh boundaries. Cheap, catches most derivative discontinuities.
  [cs.cmu.edu/~ph/discon.ps.gz](https://www.cs.cmu.edu/~ph/discon.ps.gz)
- **Lischinski–Tampieri–Greenberg 1992 (IEEE CG&A 12(6):25–39)**: discontinuity meshing for
  radiosity; classifies D0/D1/D2 discontinuities of the radiance function along event
  curves; mesh elements bounded by event curves + piecewise-quadratic interpolation.
  Key fact for us: across EV/EEE events the point-to-source factor is C0-continuous
  (derivative discontinuities only), except D0 jumps only in contact/coplanar
  configurations. So an exact inner integrand is continuous in the source point, and
  event-splitting the outer quadrature restores high-order convergence — precisely the 2-D
  lesson. [dl.acm.org/doi/10.1109/38.163622](https://dl.acm.org/doi/10.1109/38.163622)
- **Drettakis–Fiume (SIGGRAPH '94)** and **Stewart–Ghali (SIGGRAPH '94)**: *complete*
  discontinuity meshing via **backprojection** — for each face of the discontinuity mesh on
  the receiver, the topological view of the source (the "backprojection") is constant; the
  visible polygon varies continuously and analytically within a face. This is the theorem
  that makes the outer quadrature certifiable: inside an event-free panel the integrand is
  smooth (piecewise-analytic). [dl.acm.org/doi/pdf/10.1145/192161.192207](https://dl.acm.org/doi/pdf/10.1145/192161.192207),
  [Drettakis–Sillion follow-up](https://inria.hal.science/inria-00510111/file/DrettakisSillion96.pdf)
- **Durand–Drettakis–Puech, Visibility Skeleton (SIGGRAPH '97)** and **3D Visibility
  Complex (ACM TOG 21(2):176–206, 2002)**: the full catalogue of visual events for polygon
  scenes — extremal stabbing lines (VV, VEE, EEEE) as skeleton nodes, EV and EEE line
  swaths as arcs. Worst-case size O(n^4); later experimental work (Zhang thesis; Demouth
  et al., "On the Size of the 3D Visibility Skeleton: Experimental Results") measures
  worst-case Θ(n²k²) for k polytopes with n total edges, observed size ~C·k·sqrt(nk) and
  build time O(n^{3/2} k log k) in practice.
  [Sig97 PDF](https://inria.hal.science/inria-00510107/file/DurandDrettakisPuechSig97.pdf),
  [TOG 2002](https://dl.acm.org/doi/10.1145/508357.508362),
  [Springer experimental](https://link.springer.com/chapter/10.1007/978-3-540-87744-8_67),
  [thesis](https://theses.hal.science/tel-00431464)
- **Stark–Riesenfeld, "Exact Illumination in Polygonal Environments using Vertex Tracing"
  (Rendering Techniques 2000)**: irradiance from partially occluded uniform polygonal
  emitters computed EXACTLY from the set of *apparent vertices* only — no explicit polygon
  clipping, low overhead. A serious implementation alternative to boolean clipping for the
  inner integral. [Springer](https://link.springer.com/chapter/10.1007/978-3-7091-6303-0_14),
  [follow-up PDF](https://www-old.cs.utah.edu/gdc/publications/papers/mstark01a.pdf)
- **Overbeck et al., "A Real-time Beam Tracer with Application to Exact Soft Shadows"
  (EGSR 2007)**: modern evidence that exact per-point visible-polygon computation
  (beam clipping) runs at interactive rates on ~1e4–1e5-triangle scenes.
  [dl.acm.org](https://dl.acm.org/doi/10.5555/2383847.2383861)

**Which give exact per-source-point visible polygons with polygonal blockers?**
Nishita–Nakamae (projection + clipping), backprojection methods (maintain it across the
mesh), vertex tracing (implicitly, via apparent vertices), beam tracing (explicitly).
Cost per point: O(K log K + I) for K candidate blocker triangles and I shadow-edge
intersections; the global structures cost O(n²)–O(n^4) and are amortized only for static
scenes — wrong trade for a mesh rebuilt every timestep.

## 2. The inner integral: exact visible polygons from a point + Baum closed form

Pipeline per source point p and target triangle T (3-D analogue of
`_point_segment_exchange`):

1. **Cull**: back-face both ways (kernel sign), horizon-clip T against p's tangent plane
   (Phase-2 rule 2), BVH/hull prefilter of candidate blockers (analogue of
   `_candidate_blockers`; View3D's cone/cylinder + minimal-box tests, NISTIR 6925 pp. 9–10).
2. **Project** each candidate blocker triangle from p onto plane(T) (projective map;
   handle blockers crossing the p-side plane by pre-clipping the blocker against the plane
   through p parallel to plane(T) — the projective image of a clipped blocker is a convex
   polygon or an unbounded wedge; keep it as a half-plane set in the unbounded case).
3. **Boolean**: visible set = T \ union(shadow polygons). Two robust routes:
   - *Arrangement + point classification* (RECOMMENDED — it is the 2-D code's
     "cut, then classify each elementary cell once at an interior point by the exact
     predicate" trick lifted one dimension): insert all shadow-polygon edges into a 2-D
     arrangement restricted to T, then classify each face by ONE exact point-in-shadow
     test (segment p→interior-point vs blocker triangles, exact predicates). Misplaced
     events cost only extra faces, never wrong classification.
   - *Direct clipping*: Weiler–Atherton (handles concave-with-holes;
     [Weiler & Atherton, SIGGRAPH '77](https://dl.acm.org/doi/10.1145/965141.563896)) or
     Greiner–Hormann/Martinez family; View3D's appendix does successive convex splitting
     in homogeneous coordinates (NISTIR 6925 Appendix, eqs. 12–17) — simple and fast but
     convex-only and tolerance-based (`~1e-5 * area` zero tests).
4. **Evaluate**: for each visible piece, Baum/Hottel point-to-polygon closed form
   F = (1/2pi) sum_i atan2(|R_i x R_{i+1}|, R_i . R_{i+1}) (n . unit(R_i x R_{i+1}))
   (Phase-2 rule 1; [Baum–Rushmeier–Winget SIGGRAPH '89, doi
   10.1145/74333.74367](https://doi.org/10.1145/74333.74367); van Oosterom–Strackee form).
   Sum of pieces = EXACT point-to-target factor under occlusion. Zero rays fired.

**Robustness practice**: Shewchuk adaptive exact predicates (orient2d/orient3d,
incircle) for every orientation/side test — the classifier and the arrangement must agree,
which is exactly what tolerance-cascade code fails at in grazing configurations.
[Shewchuk, Discrete & Comput. Geom. 18:305–363, 1997](https://link.springer.com/article/10.1007/PL00009321),
[code](https://people.eecs.berkeley.edu/~jrs/papers/robustr.pdf),
maintained ports: [libigl-predicates](https://github.com/libigl/libigl-predicates),
[georust/robust](https://github.com/georust/robust).
Degeneracies (blocker edge collinear with p and a target vertex; coplanar blocker/target;
shared mesh edges) — handle as in 2-D: shared-feature pairs get V=1-or-culled analytically
(Phase-2 rule 5), everything else measure-zero-classified by the interior-point rule.
Complexity per point: O(K) projections + O((K+I) log K) arrangement + O(K) Baum edges;
in a trench/hole mesh K is 10–100 after culling, so ~1–10 µs/point in compiled code
(Overbeck's beam tracer supports this scale of optimism; pure Python will not).

## 3. The outer integral: certified 2-D quadrature with analytic event curves

The visibility-interval structure on T changes, as p moves in the source triangle S, only
across **visual event surfaces** (aspect-graph theory):

- **EV events**: p collinear-planar with (blocker vertex, target edge) or (blocker edge,
  target vertex) — each is a PLANE in p-space; its trace on plane(S) is a LINE. Count:
  O(K·(v_B·e_T + e_B·v_T)) ≈ O(9K + 9K) lines per pair. These are the exact analogue of
  `_panel_events`' linear roots and are enumerable with rational predicates.
- **EEE events**: p on the ruled QUADRIC spanned by three edges (two blocker edges + one
  target edge, or three blocker edges); trace on plane(S) is a CONIC arc. Count O(K³)
  worst-case triples, but only triples whose swath actually crosses S matter.
  [Gigus–Malik, IEEE TPAMI 1990](https://dl.acm.org/doi/abs/10.1109/34.44399);
  [Plantinga–Dyer, IJCV 5:137–160, 1990](https://link.springer.com/article/10.1007/BF00054919)
  (aspect-graph sizes O(n^6) orthographic / O(n^9) perspective — a warning against global
  enumeration, not against per-pair local use).
- Facing/horizon clips: linear in p — same as the 2-D facing-clip roots.

**Recommended scheme (mirrors the 2-D architecture + its fallback):**
split S by the EV lines only (cheap, exact, catches the overwhelming majority of kinks —
Heckbert's "incomplete discontinuity meshing" finding); run certified adaptive triangle
quadrature (Simpson/cubature pair with embedded error, per-panel budgets summed into the
pair receipt exactly as `_analytic_pair_exchange` does); do NOT enumerate EEE conics in
v1 — inside an EV-panel the only possible surprise is an EEE kink, which the adaptive
pair-difference detects and refines toward, *without bias*, because the inner value is
exact everywhere. Escalation ladder if receipts refuse to close: (a) split along the
suspected EEE conic (root-find the tangency along the refinement direction — 1-D,
certified), (b) fall back to the conservative 2-D-cell shadow refinement
(`adaptive_refinement` analogue) for that pair. Keep the RuntimeError-on-budget-exhaustion
contract.

Store the symmetric quantity S_ij = A_i F_ij computed once per pair (double integral is
symmetric) — reciprocity by construction, as with `H_ij`; escape = open-boundary row
deficit; same closure gates.

## 4. Prior art in heat transfer / transport codes

- **View3D (Walton, NIST)** — read in full today
  ([NISTIR 6925, 2002](https://nvlpubs.nist.gov/nistpubs/Legacy/IR/nistir6925.pdf);
  [Walton 1986, NBSIR 86-3463]; code: [view3d.sourceforge.net](https://view3d.sourceforge.net/),
  [github.com/jasondegraw/View3D](https://github.com/jasondegraw/View3D)):
  - Obstructed pairs: **1AI with shadow projection** — project obstruction shadows from
    integration points on the *nearer* surface onto plane of the other, clip the target to
    its unshaded portion (convex-polygon processing in homogeneous coordinates, Appendix),
    then the EXACT Hottel angle formula (eq. 6–7, atan2-stabilized) around the unshaded
    contour. So View3D's inner is exact-per-point, like ours.
  - Outer: adaptive Gaussian integration — 9- vs 16-pt (parallelogram) or 7- vs 13-pt
    (triangle) forms compared; if |AF^[hi] − AF^[lo]| ≥ eps·A_min, split into 4 congruent
    subsurfaces, recurse (eq. 9, Fig. 12). NOT event-aware.
  - Culling gauntlet: exclude never-obstructing surfaces; per-row reduced lists; cylinder/
    cone test on enclosing spheres; minimal-box test; behind-plane classifications;
    projection-direction choice; obstruction grouping ("View3D-O") — 2/3 of runtime was
    obstruction processing in their BB2510 benchmark (4350 surfaces, 762 s → 264 s grouped,
    866 MHz Pentium).
  - Validation style to copy: **Shapiro's analytic obstructed case** (two opposed unit
    squares, back-to-back half-squares between; F_12 = 0.11562061 exactly; View3D hits
    −6e-8 error at eps=1e-4 with 125 points — Table 1), rowsum(=closure) error as the
    figure of merit across 174–4350-surface enclosures (rowsum < 1e-3 at eps=1e-4),
    2AI-with-blockage MC comparison converging like sqrt (Table 2).
  - Their documented failure = our finding: Fig. 15 grazing sliver where all Gauss points
    are blocked → F set to 0, adaptive refinement never triggers, rowsums stall (Table 3,
    shifts 0.01–0.04). Event-split panels are precisely the cure.
- **Mitalas–Stephenson 1966** (1LI: one contour integral analytic, one numeric) and
  **Schröder–Hanrahan 1993** closed form (0LI; dilogarithms) — unobstructed oracles only;
  0LI cost blows up on skewed edges (28–218 s / 100k factors vs 2.5 s adaptive, NISTIR
  Fig. 5) — keep as TEST ORACLE (Phase-2 rule 4).
  [multires.caltech.edu/pubs/ffpaper.pdf](https://www.multires.caltech.edu/pubs/ffpaper.pdf)
- **Modest, Radiative Heat Transfer (3rd ed., ch. 4 view factors)**: for obstructed cases
  recommends view-factor algebra / unit-sphere (Nusselt) projection / contour integration
  where possible and otherwise statistical (MC) sampling — i.e., the textbook default IS
  sampled visibility; no certified deterministic occlusion operator is standard in the heat
  transfer literature. [3rd ed. PDF](https://ostad.nit.ac.ir/payaidea/ospic/file8971.pdf);
  comparative study: Emery et al., ASME J. Heat Transfer 113:413–422 (1991).
- **Modern codes**: Chaparral (Sandia; adaptive 0LI/2AI-with-blockage rays, hemicube) —
  View3D beat it on accuracy and often speed (NISTIR Figs. 16–18);
  [pyviewfactor](https://pypi.org/project/pyviewfactor/) (double contour integral +
  ray-based obstruction tests — sampled visibility again); OpenGL hemicube-style obstructed
  view factors [(Building Simulation 2014)](https://www.tandfonline.com/doi/full/10.1080/19401493.2014.917700)
  (raster visibility = biased at slivers by construction);
  [MOOSE heat transfer module](https://mooseframework.inl.gov/modules/heat_transfer/)
  (angular-quadrature ray tracing). **No production code found that does event-certified
  exact occlusion; the pieces all exist separately (View3D inner, graphics event theory).**
- **Clausing-factor literature** (free-molecular tube transmission = our hole problem with
  diffuse walls): Clausing's Fredholm integral equation of the 2nd kind on the tube wall;
  modern high-accuracy solutions via singularity-subtracted Gauss–Legendre
  ([JVST 4:360 (1967)](https://pubs.aip.org/avs/jvst/article/4/6/360/248712/),
  [JVST A 25:758 (2007), "Efficient numerical solution of the Clausing problem"](https://pubs.aip.org/avs/jva/article-abstract/25/4/758/101816/));
  analytic ring–ring kernels: Leuenberger–Person concentric-cylinder view factors
  ([AIAA J., doi 10.2514/3.6964](https://arc.aiaa.org/doi/10.2514/3.6964)); finite-length
  ring/annulus enclosure factors
  ([JQSRT 2012](https://www.sciencedirect.com/science/article/abs/pii/S0022407312003196)).
  These supply closed-form kernels AND reference transmission probabilities for the
  axisymmetric validation gates.

## 5. The azimuthally-symmetric reduction (the 80/20 for 200:1 holes)

Let the hole be a body of revolution with piecewise-linear profile r(z) (stack of conical
frusta; kink circles at radii r_k, heights z_k; bottom disk; annular rings at ledges).
Mesh = axial bands; assume azimuthally uniform flux (valid for symmetric drift/charging).
Then the EXACT 3-D occlusion operator collapses to closed form:

- Take source point p = (r_s, 0, z_s), target point q = (r_t cos DPhi, r_t sin DPhi, z_t).
  The connecting chord has squared cylinder-radius
  rho(u)^2 = (1−u)^2 r_s^2 + u^2 r_t^2 + 2u(1−u) r_s r_t cos DPhi,  z(u) linear in u.
  rho(u) is convex; r(z(u)) is piecewise linear. So g(u) = rho(u) − r(z(u)) is convex on
  each frustum sub-interval and attains its maximum at sub-interval endpoints, i.e. AT THE
  KINK PARAMETERS u_k = (z_k − z_s)/(z_t − z_s), which do not depend on DPhi.
- **Exact blocking test**: chord blocked  iff  exists kink k with z_k strictly between
  z_s and z_t and rho(u_k) > r_k. Since rho(u_k)^2 is increasing in cos DPhi, each kink
  contributes a closed-form threshold
  C_k = (r_k^2 − (1−u_k)^2 r_s^2 − u_k^2 r_t^2) / (2 u_k (1−u_k) r_s r_t),
  and the visible azimuth set is the single symmetric interval
  cos DPhi <= min_k C_k, i.e. DPhi in [DPhi*, 2pi − DPhi*], DPhi* = arccos(clamp(min_k C_k)).
  **Occlusion events are arccos of rational functions of the kink radii — the exact
  analogue of the 2-D "blocker endpoint projection" cuts.** A straight or monotonically
  tapered hole has no occlusion (all C_k >= 1); events appear exactly with necking/bowing.
- Facing clips (cos theta > 0 on conical normals) add further cos DPhi thresholds — also
  closed form.
- Operator assembly: inner-inner = kernel integral over the VISIBLE DPhi interval
  (algebraic function of cos DPhi; closed form in the cylindrical cases per
  Leuenberger–Person, else certified 1-D quadrature); outer = certified 2-D quadrature over
  generator coordinates (s, t) with panels split where any C_k(s,t) crosses ±1 or the
  min-attaining k changes — smooth algebraic event curves in the (s,t) square. Cost is
  within a small constant of the existing 2-D operator; exactness is inherited, receipts
  identical in kind. This is the first thing to build (it also directly serves the R4/
  charging deep-hole frontier and the Clausing validation ladder).

## 6. Complexity and scaling for petch (N ~ 1e3–1e4 triangles, rebuilt every step)

| approach | per-pair cost | full build (N=2e3, ~10% facing pairs) | bias receipt |
|---|---|---|---|
| current 2-D analytic (trench mean-profile) | ~events + O(10–100) exact inner evals | seconds (CPU, measured; 54x over refinement) | certified, no visibility bias |
| axisymmetric exact 3-D (§5) | ~2-D cost x small const | seconds–tens of s | certified, no visibility bias |
| general 3-D exact (this report) | 30–300 outer pts x 1–10 µs exact inner | ~2e5 pairs x 0.1–1 ms ≈ 20–200 s compiled, parallel | certified, no visibility bias |
| Gauss/ray-sampled visibility (Phase-2 rule 6) | 5–64 rays/pair | ~3–10x cheaper than exact | UNRECEIPTED sliver bias (today's finding; Walton Fig. 15) |
| visibility skeleton / complete disc. mesh | global O(n²k²) worst, ~k·sqrt(nk) obs. | minutes+, not incremental | exact but wrong amortization for moving meshes |

Levers that keep general-3-D honest at 1e4 faces: View3D's culling gauntlet + obstruction
grouping (their measured 2/3 obstruction-processing share and 762→264 s win say this is
where the time goes); Phase-2 rule 9 caching (re-clip only pairs whose blockers moved);
EV-line splitting only on pairs whose first-pass receipt fails; hierarchy
(Hanrahan-style) only past ~1e4. Periodic-image handling carries over unchanged
(Phase-2 rule 7) — wrapped shadow projections instead of wrapped shadow rays.

## 7. Validation plan (mirrors today's 2-D campaign)

1. **Analytic gates**: Shapiro obstructed case (F = 0.11562061; NISTIR 6925 Table 1);
   parallel/perpendicular square catalogues; Schröder–Hanrahan 0LI oracle on unoccluded
   pairs; cylinder/ring closed forms (Leuenberger–Person; JQSRT 2012 annulus set).
2. **Grazing-sliver adversarial suite** (the point of the whole exercise): Walton Problem
   Case 2 geometry (Fig. 15) + randomized near-tangent blocker families; gate = exact
   operator returns F > 0 with closing receipt where Gauss-visibility classifiers return
   0/overcount; document the sampled-vs-exact gap exactly as the 2-D module tests do.
3. **Dense-reference convergence**: tolerance ladder (1e-3 → 1e-7) must show monotone
   receipt closure; refusal must raise, never silently degrade (contract carried over from
   `_analytic_pair_exchange`).
4. **MC cross-check**: stratified cos-weighted ray pairs with CI receipts; agreement within
   3 sigma on every pair class (open, partial, sliver).
5. **Closure/reciprocity gates**: S_ij symmetric by construction; rowsum + escape = 1
   within tolerance (Walton's rowsum figure of merit, threshold 1e-3 at eps=1e-4 as the
   floor, our 5e-13 style for the final operator); bounded renormalization concession as in
   the 2-D builder.
6. **Dimensional consistency**: long-extruded trench mesh → recover
   `deterministic_exchange_2d` per-unit-depth exchange in the L→inf limit;
   axisymmetric operator → Clausing transmission probabilities for straight tubes
   (JVST 1967/2007 reference solutions) and → the 2-D operator on the (r,z) section
   in the slit limit.

## What to build first (ordered)

1. **Axisymmetric exact-occlusion operator** (§5) in a new module beside
   `deterministic_exchange_2d.py` — same dataclass contract (symmetric exchange, escape
   row, receipts, fingerprint, RuntimeError on refusal). Unlocks 200:1 holes now.
2. **Exact inner kernel for triangles**: Baum point-to-polygon with horizon clip +
   arrangement-and-classify occlusion with Shewchuk predicates; gate vs MC + Shapiro.
3. **EV-split certified outer quadrature** per pair, with conservative-refinement fallback;
   port the 2-D error-budget bookkeeping verbatim.
4. **Culling/caching layer** (View3D gauntlet + move-tolerance link cache) before any
   hierarchy work; hierarchy only if profiles show >1e4-face meshes dominating.
5. Only then consider EEE conic enumeration and skeleton-style acceleration — driven by
   receipts, not speculation.
