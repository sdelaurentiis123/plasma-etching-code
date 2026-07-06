# Full-Poisson charging solver — implementation-ready handoff plan (v2, 2026-07-07)

**Goal:** replace the convention-dependent Dirichlet charge→potential map with a PHYSICAL-UNITS
self-consistent variable-ε Poisson solve (Kushner MCFPM route). This makes voltage labels ABSOLUTE
— the ground truth that arbitrates the HG label question (C14: the same converged state reads
11–460 V through defensible variants of HG's under-specified σ→V kernel, while all physics
observables are invariant; a units-honest Poisson solve ends that ambiguity). It also buys:
thick-oxide capacitance physics (substrate coupling), the cryo surface-conductivity module in real
units, and the 3D path.

**Everything below is informed by what this campaign already built, broke, and proved. Read
HG_RESIDUAL_HUNT.md + PHYSICS.md §5 + FRONTIER_LOOP.md (cycle log) first.**

## 0. Physical scale (HG reference conditions — precomputed, use these)
- n = 1e12 cm⁻³, Te = 4 eV, Ti = 0.5 eV, M = 35.45 amu (Cl), V_s(t) = 37 + 30·sin(2π·400 kHz·t)
- Bohm speed u_B = √(eTe/M) = 3.3e3 m/s → ion flux J_i = e·n·u_B ≈ 0.53 A/m² (53 µA/cm²)
- λ_D = 743·√(Te/n[cm⁻³]) cm ≈ 15 µm ≫ W = 0.5 µm → **the gas is charge-free: Laplace-in-gas is
  EXACT in the gas.** The Poisson content is the DIELECTRIC interior + substrate coupling, not the
  gas. (Do not "fix" the gas solve.)
- Grid: W = 0.5 µm / 16 cells → h = 31.25 nm. HG oxide under the floor: 1.8 µm (58 cells).
- Charging timescale: σ(40 V across ~0.5 µm) ~ ε₀·8e7 V/m ≈ 7e-4 C/m² → τ ≈ σ/J_i ≈ 1.3 ms
  ≈ 500 RF cycles. Quasi-static accumulation ✓ (the loop's per-iteration Δσ maps to real time).
- ε_r: SiO2 3.9, photoresist 1.6, poly-Si = floating conductor, Si substrate = grounded conductor.

## 1. Where the code stands (all in src/petch/charging_general.py)
- `poisson(V, sweeps, omega)` — flux-conservative variable-ε red-black SOR, face-averaged ε
  (arithmetic mean), `∇·(ε∇φ) = -k·ρ` with the `rho_coupling` FUDGE k (to be replaced by real units).
  The ε-jump condition ε₁E₁−ε₂E₂=σ is automatic in this discretization (audit-verified).
- `GROUND` material + `add_grounded_substrate(mat, ox_cells, sub_cells)` — dielectric stack on a
  grounded substrate; solvers pin GROUND to 0. Built, tested inert-by-default.
- `poisson_step` — under-relaxation of the charge update (tames but does not fix the instability).
- `charge_update="log"` — quasi-Newton potential-space update (Laplace mode only; see P4 for the
  Poisson-mode analog).
- KNOWN-BROKEN naive path: `rho[insul] += scale*net[insul]` accumulates unboundedly (walls ran to
  −1000s of V; C8). Root causes: no physical scale, charge smeared over the dielectric body, no
  self-limiting dynamics.
- Geometry: use the CORRECTED HG stack (C13): `_build_edge_array_geometry(AR, poly_um = AR*0.5-0.54)`
  — PR fixed 0.54 µm, poly grows. The default builder is INVERTED vs HG; do not gate against HG
  with the default.

## 2. Build steps (P1–P7, in order; each with its own check)

### P1 — Physical units
Work in volts and real charge. Per launched macro-particle, the statistical weight is
  w = J_i · A_column · Δt_iter / (e · n_per_iter_per_column)
but the cleaner normalization: define the PER-ITERATION real fluence F_iter [ions/m² per iteration]
as a solver input (default: F_iter chosen so ~50 iterations ≈ one charging time τ). Deposited charge
per cell: Δσ_j [C/m²] = e·F_iter·(counts_j − e_counts_j)/(launched per column). RHS of the discrete
Poisson (2D, cell size h): the flux-conservative stencil solves Σ_faces ε_face(V_nb−V_c) = −q_c/ε₀
with q_c = σ_c·h (2D line-charge density per unit y). CHECK: a parallel-plate slab (uniform σ on a
plane above the grounded substrate through ε_r oxide of thickness d) must read V = σd/(ε₀ε_r) to <2%.

### P2 — σ-sheet on the interface cell
Deposit charge ONLY on the gas-facing surface layer (the cell where the particle lands), never
smear into the dielectric body. Bookkeep σ as a per-surface-cell array (not the volumetric rho
grid); inject into the Poisson RHS at those cells. Interior dielectric cells carry ρ=0 (no bulk
conduction at room T; the cryo module later moves σ laterally).

### P3 — Capacitance-matched substrate
Domain economy: instead of gridding the full 1.8 µm oxide (58 rows), grid ~15 rows and scale the
bottom 2–3 rows' ε so the feature→ground capacitance equals the true-thickness value
(Kushner JVST A 37, 031304 p.8, verbatim trick). CHECK: the P1 slab test with the matched stack
reproduces the full-thickness V to <5%.

### P4 — Stable charge dynamics (the part that killed C8 — do it this way)
The instability is a lagged stiff loop: charge → field → (next iter) fluxes. Two required pieces:
1. **Diagonal-capacitance Newton step**: precompute C_jj ≈ ΔV_j/Δσ_j numerically ONCE (unit charge
   on cell j → solve → read V_j; or the cheap estimate C_jj ~ ε₀ε_eff/h per unit area). Then update
   Δσ_j = C_jj⁻¹ · anneal · clip(Te·ln(Γi_j/Γe_j), ±2Te)  — the charge-space version of the proven
   log current-balance update (drives each cell toward local balance in VOLTS, which is what the
   fluxes respond to).
2. **Physical bound**: V_j may not exceed (V_dc+V_rf) nor go below −10·Te (Maxwellian tail reach)
   — enforce by capping Δσ when the solved V would cross (NOT a hard V clip; cap the CHARGE step).
Robbins–Monro anneal + Polyak tail-averaging as in the Laplace mode. CHECK: AR4 corrected-stack run
converges with Vmin ≥ −45 V, no runaway, and the scheme-independent observables reproduce the
Laplace-mode values (floor flux 0.21±0.03, edge 7.8±1.5, foot E 28±3) — the physics must be
map-invariant (we proved the fixed point is; this is the regression gate).

### P5 — Conductors
Floating poly lines: keep the whole-component potential-space update (equipotential Dirichlet with
Vc from the component log-balance — already correct and convention-free). Grounded substrate:
Dirichlet 0 (already wired). The conductor's induced surface charge is an OUTPUT (Gauss law over
its faces), not a state.

### P6 — Absolute-label readout (the payoff)
With P1–P5 green, the floor potential IS the physical answer in volts for HG's stated conditions.
Report it. Whatever it reads (∼49 per the ion's-eye argument, or between), it arbitrates the C14
convention question with units-honest electrostatics. Update HG_RESIDUAL_HUNT.md + PHYSICS.md §5
with the number. Then re-run `petch_floor_profile` (the notch table) through the Poisson mode —
the C11 AR4 notch-trend question (gates B/C) gets its first-principles retest with absolute Vf.

### P7 — Follow-ons (do not block on these)
- Cryo module re-wire: surface conductivity moves REAL σ laterally (units now exist).
- Multigrid replacing SOR only if profiling demands; then the GPU port (device-resident).
- 3D: the same discretization extends; the σ-sheet bookkeeping is already surface-based.

## 3. Pitfalls already paid for (do not re-pay)
1. Naive `rho += net` runaway (C8) — the reason P4 exists.
2. Anneal-frozen transients pinning clips (C6/C9) — the reason the log/Newton step exists.
3. Periodic-wrap tracer BCs (C9a) — mirror BCs are in; don't revert.
4. Stack inversion (C13) — use poly_um = AR·0.5−0.54 for HG work.
5. The `insul_vguard` V-clip is a Laplace-mode legacy — in Poisson mode bound the CHARGE step (P4.2).
6. Keep `field_model="laplace"` the DEFAULT until all P-checks pass (26-test suite must stay green;
   the Poisson mode is opt-in until gated).
7. HG's published labels are NOT a gate for absolute volts (proven under-specified + internally
   inconsistent; HG_RESIDUAL_HUNT.md) — gate on the P1/P3 analytic checks and the observable
   invariance, not on "33".

## 4. Effort estimate
P1–P2: ~1–2 h. P3: ~1 h. P4: the real work, ~2–4 h with the checks. P5: mostly done. P6: runs.
Total: a focused session. Everything is CPU-fine at these grids (~2 min/run at nit800).
