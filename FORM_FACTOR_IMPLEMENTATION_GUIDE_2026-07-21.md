# Form factors for 3-D triangle radiosity: implementation guide (Phase 2 reference)

(Deep-research agent output, 2026-07-21. Companion to ENGINEERING_PLAN_NOISE_FREE_TRANSPORT.md.)

Converged architecture: **F_ij = F_ij^analytic-unoccluded x V_ij(ratio-estimated)** — outer
Gauss quadrature over the receiver triangle x exact point-to-polygon (Lambert/Baum) inner
factor; stratified ray-pair visibility; first-shell periodic images with wrapped shadow rays;
back-face cull + disc-formula thresholds; Hanrahan '91 hierarchy only past ~1e4 faces. This is
what NIST View3D and hierarchical-radiosity practice converged on.

## Key rules (load-bearing)

1. Point-to-polygon (Baum/Hottel, exact): F = (1/2pi) sum_i atan2(|R_i x R_{i+1}|,
   R_i . R_{i+1}) * n . unit(R_i x R_{i+1}). NEVER acos — atan2 of cross/dot (Walton NISTIR
   6925 stability note; van Oosterom-Strackee pattern). Guard collinear edge pairs (skip,
   limit 0). Clamp F to [0,1].
2. ALWAYS horizon-clip the source polygon against the receiver's tangent plane (lift by
   eps) — signed-cosine contributions otherwise go negative on meshes. Trivial culls first.
3. Point in polygon plane -> F=0; interpenetrating in-plane -> treat as touching, F=0.
4. Polygon-to-polygon: Schroder-Hanrahan closed form (dilogarithms, complex branch indices)
   is a TEST ORACLE, not production (View3D benchmarks: 28-218 s vs 2.5 s per 100k factors;
   Gauss-outer x analytic-inner reaches <1e-6 at N=4 with Gaussian points). Also
   Narayanaswamy IJHMT 2015 as a modern reference form.
5. Shared-edge/adjacent triangles (the common mesh case): interior triangle-Gauss points +
   analytic inner factor converge fine (<1e-6 at N=4); never pure double-line/double-area
   numeric integration (r->0 on the shared edge). Skip ray-cast visibility for shared-feature
   pairs (V=1 at reflex dihedral; back-face-culled otherwise). Optional: one forced
   subdivision of the strip along the shared edge. 1LI-a (analytic shared-edge pair term) is
   the View3D alternative.
6. Visibility ratio estimator: probe 4-5 rays (corners + centroid); all-visible V=1,
   all-blocked drop link; partial -> 16 stratified jittered pairs (32-64 only for
   energy-heavy links). Weight samples by the point kernel cos cos / r^2 to kill
   correlation bias at grazing angles. Offset origins/endpoints by eps*edge-length along
   normals; exclude endpoint triangles by ID, never by t-epsilon.
7. Periodic lateral domain: sum source-triangle images over lattice shifts; kernel decay is
   r^-4 for in-layer elements (cos factors ~ h/r each), so nearest image + first shell
   (+/-1; +/-2 for tall open structures) suffices; corridor pairs are self-occluded by
   replicas. Wrapped (minimum-image) shadow rays are equivalent to infinite images and
   truncate via occlusion; attribute each ray to its image cell. Do the image sum in
   area/point-to-polygon form (closed contours), never per unclosed edge. Reciprocity holds
   per image PAIR (m with -m).
8. Culling gauntlet: back-face (all vertex-pair directions) -> disc estimate
   F = A cos cos/(pi r^2 + A) with threshold (accumulate dropped mass into the open/escape
   row term) -> quadrature order by relative separation s = dist/(r_enc_i + r_enc_j): s>3 ->
   1-pt, 1<s<=3 -> 3-pt, else 7-pt + adaptive doubling check (|A F^[k+1] - A F^[k]| <
   eps A_min). Hierarchy (Hanrahan-Salzman-Aupperle refine oracle, BF refinement,
   O(N) links) only when N > ~1e4 faces.
9. Etch-specific: cache links + quadrature order + visibility class across time steps;
   re-evaluate only links whose triangles moved beyond tolerance; re-run visibility only for
   partial links near the moving front.
10. Post-assembly hygiene: F in [0,1]; reciprocity symmetrize F~ = (A_iF_ij + A_jF_ji)/2A_i;
    row closure sum_j F_ij + F_open <= 1 with small-deficit renormalization.

## Sources

Schroder-Hanrahan SIGGRAPH '93 (multires.caltech.edu/pubs/ffpaper.pdf); Walton NISTIR 6925 /
View3D (nvlpubs.nist.gov/nistpubs/Legacy/IR/nistir6925.pdf) — the single most useful
document; Baum-Rushmeier-Winget SIGGRAPH '89 (doi 10.1145/74333.74367); Wallace-Elmquist-
Haines SIGGRAPH '89 (doi 10.1145/74334.74366); Hanrahan-Salzman-Aupperle SIGGRAPH '91
(graphics.stanford.edu/papers/rad/); van Oosterom-Strackee IEEE TBME 1983; Manstetten TU Wien
dissertation 2018 (radiosity-form etch flux, Jacobi = re-emission orders); Narayanaswamy
IJHMT 2015; Ertl et al. Microelectronic Eng. 2009 (periodic lateral BC practice).
