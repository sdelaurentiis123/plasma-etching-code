# Charging floor over-charge: root cause + first-principles fix plan (C8, 2026-07-06)

The last charging-accuracy gap: the AR-4 trench floor sits at **37 V** (converged, correct integrator)
vs Hwang-Giapis **33 V**. This note records the fully-diagnosed root cause (2 primary-source research
passes + a line-referenced code audit) and the concrete build to close it. SEE is NOT the fix (ruled out).

## Root cause: LOCAL vs GLOBAL charge -> surface-potential map (NOT the interior PDE)

Confirmed against the actual papers:
- **HG solve Laplace in the gas — same equation we do.** So Laplace-in-gas is not the bug.
  (HG JVST B 15,70 p.75: "the Laplace equation ∇²V=0 is solved ... FDM"; sheath solved separately.)
- The focusing ("electrostatics decreases the geometric shadowing", HG JAP 82,567) is a **LOCAL,
  in-trench** effect — the positive floor is an attractive well that bends near-mouth electrons onto
  the floor. HG: the dipole field "decays very fast ... felt only very close to the microstructure."
  (My earlier "long-range fringing field above the mouth" hypothesis was WRONG.)
- The gap is entirely in the **charge -> surface-potential map** (the Dirichlet values fed to Laplace):
  - **We map charge to potential LOCALLY:** `Vs[cell] += net[cell]` (charging_general.py:445). Each
    insulator cell's potential comes from its OWN accumulated current only.
  - **HG/Kushner use a GLOBAL, ε-aware map:** each surface cell's potential depends on ALL deposited
    charges through the dielectrics (HG: Coulomb superposition with per-material ε + mirror images;
    Kushner: full variable-ε Poisson `∇·(εE)=ρ` over gas+solids, ρ in the cells). The ε-jump
    `ε1E1-ε2E2=σ` falls out of the finite-volume form automatically.
  A local map flattens the lateral surface-potential profile, so the harmonic gas field has no
  inward-bending component -> electrons reach the floor by pure geometry -> floor climbs to 37 V to
  balance the geometric-only electron deficit. The two symptoms (over-charge AND e_traced == geometric)
  are one defect. Verified numerically: at fine grid + converged integrator, e_traced -> 0.124 = geometric
  exactly (C7).

## What is built (this session, committed, correct infrastructure)
- `GROUND` material (Dirichlet V=0) + `add_grounded_substrate(mat, ox_cells, sub_cells)` — extends the
  domain with an oxide dielectric stack on a grounded Si substrate (Kushner's bottom-ground route).
- Field solvers pin GROUND to 0 (laplace `apply_bc`, poisson sweep); `poisson_inside` excludes it.
- Backward-compatible: inert when no GROUND cells exist; 4/4 charging tests pass.

## What is NOT working yet (the hard part)
Turning on `field_model="poisson"` + substrate is **UNSTABLE**: electron-collecting walls run to
-1000s of V (Vmin -547..-5505), no focusing emerges (floorV still ~40). Cause: the charge update
`rho[insul] += scale*net` (charging_general.py:~443) accumulates **unboundedly** with no physical
scaling, and the decaying anneal freezes early overshoots. A naive rho clip did not bind. This is the
core of what makes Kushner MCFPM a decades-long codebase — it is real numerics work, not a one-liner.

## Update: under-relaxation stabilizes but does NOT focus (confirms the full build is needed)
Added `poisson_step` (under-relaxes the lagged charge update). It tames the runaway (Vmin -1922 -> -110
at poisson_step=0.03) but focusing still does NOT emerge (e_traced stays ~0.13 = geometric; floor either
stays 37 or under-converges). So stability is necessary but not sufficient — the focusing needs the
CHARGE-DEPOSITION physics below (interface-sigma + physical units + capacitance match + conductor
corner charge), not just a stable solve. This is a multi-hour numerics build.

## The concrete first-principles build to land floor = 33 V (Kushner route, recommended)
1. **Physical-unit charge.** Replace the `rho_coupling` fudge with `ρ·h²/ε₀` scaling (h = cell size in
   meters, ρ = accumulated particle charge × statistical weight / cell volume). This makes the σ->V map
   first-principles and self-consistently bounded (a given σ maps to a definite, physical V).
2. **Charge as a σ-SHEET on the interface cell**, not smeared over the dielectric body (audit P2):
   deposit `net` only on the gas-facing surface layer of each insulator, not `rho[insul]` (whole body).
3. **Capacitance-matched substrate** (Kushner p.031304-8): tune the bottom 2-3 substrate rows' ε so the
   feature-to-ground capacitance matches the true oxide thickness. Sets the floor's absolute potential.
4. **Stable charge dynamics.** Under-relaxed, non-decaying-to-zero step; monitor per-surface net->0.
   Consider implicit/damped update (the gain scale × d(net)/d(rho) must be <1 on high-flux walls).
5. **Uneven conductor equipotential** (HG JVST B p.75, "critically important"): the poly-Si line must be
   equipotential with charge PILED at the poly/SiO2 corner (not equidistributed) — needed for the field
   that focuses electrons and drives the neighbour rise.
6. Then **remove the band-aids** (audit P4): `insulator_e_focus`, `vf_focus_pot`, `open_wall_boost`, the
   `np.maximum` geometric electron floor, `insul_vguard` clip — they fake the focusing the correct field
   will now produce; leaving them on double-counts.

**Alternative (HG's exact route):** keep Laplace-in-gas, replace the local `Vs+=net` BC with a GLOBAL
ε-weighted Coulomb superposition (Green's function from all deposited charges + mirror images about the
side centerlines). Same physics, O(N²), GPU/differentiable-friendly, and stable (no runaway) — but more
code than turning on the (fixed) Poisson.

## Gate
AR 1->4: floor ion flux vs HG (0.59/0.40/0.29/0.22), floorV -> 33 at AR4, edge -> 7, neigh -> 39, with
ALL focus knobs OFF. e_traced must exceed geometric (0.124) — that is the focusing emerging from the field.

## Sources
- HG JVST B 15,70 (1997): https://authors.library.caltech.edu/records/ac5xn-zqb88 (methods, charge update)
- HG JAP 82,566 (1997): https://authors.library.caltech.edu/records/je8bd-j6v68 (floor=33 V benchmark)
- Huang/Kushner MCFPM, JVST A 37,031304 (2019): https://cpseg.eecs.umich.edu/pub/articles/JVSTA_37_031304_2019.pdf
