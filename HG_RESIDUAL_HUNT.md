# The HG residual hunt (C13/C14, 2026-07-07): one suspect left standing

After the stack-inversion fix (edge and flux SNAP onto HG), the remaining deltas were floorV
(49 vs 33) and neighbor (14 vs 39). Systematic elimination with frozen-field ion traces
(derived source, corrected stack, AR4):

| frozen configuration                     | floor flux | verdict |
|---|---|---|
| floor 33, walls 7/7 (symmetric)          | 0.388 | pure energy rejection = bathtub P(E>33) ✓ ion physics clean |
| floor 33, walls 7/39 (published split)   | 0.000 | 32V cross-field sweeps EVERYTHING — HG's numbers are mutually inconsistent in a Laplace field if they describe one trench |
| floor 33, walls 39/39 (interior trench)  | 0.474 | symmetric + walls FUNNEL ions (flux UP) — interior-trench hypothesis DEAD |
| floor 33 + 60V foot horns (HG Fig 5)     | 0.389 | horns only shave corner columns — DEAD |
| floor 49, walls 7/7                      | 0.232 | = HG's flux. Our self-consistent state (floorV 49.1, flux 0.208) reproduces HG's PHYSICAL observable at a different voltage LABEL |

ELIMINATED: electron injection convention (matters, shifts floorV 49<->29, but no convention gives
the joint state); stack geometry (fixed, closed edge+flux); interior-trench measurement target;
foot-horn aperture squeeze; ion energy distribution (verified analytic); tracer dynamics (verified
against analytics).

LAST SUSPECT STANDING: the charge -> surface-potential MAP CONVENTION. HG compute surface
potentials by global Coulomb superposition with per-material epsilon + mirror images (JVST B p.75);
a charge at a gas/dielectric interface reads differently through an eps-weighted map than through
our per-cell Dirichlet Vs (interface factor 2/(1+eps_r) = 0.41 for SiO2; our observed label ratio
33/49 = 0.67 sits in the convention-dependent range). The physical observables agree:
flux 0.21-0.23 (HG 0.22), edge 7.8 (7), foot E 28.0 (28), foot horns ~61 (~60). The VOLTS
disagree in exactly the way a map convention would produce.

DECISIVE EXPERIMENT (queued, C14): implement HG's exact sigma->V map (field_model="hg_coulomb":
global eps-weighted Coulomb superposition + mirror images, per the C8 audit scoping) and read our
converged charge state through THEIR map. If it reads ~33 at the floor, the entire HG comparison
closes: full observable agreement + the voltage-label difference attributed to a documented
convention. The neighbor (14 vs 39) then gets re-examined under the same map with a full multi-line
array (their pattern has many lines; ours has 2.5).
