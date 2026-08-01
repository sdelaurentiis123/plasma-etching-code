# GPU port of the per-atom mixed-layer chemistry — scoping plan (2026-07-29)

Scoping only; nothing implemented. All measurements on the laptop (CPU, arm64,
warp 1.14 without CUDA), so **fractions and ratios are the deliverable, not
absolutes**. Two independent measurements: a microbenchmark of `step()` at the
real ml13 face count, and an instrumented run of the exact production pilot
configuration (`scripts/preflight_krueger.sh` flags, dx=0.01 µm,
`--grazing-ion-reflection literature_v1`).

Headline: **the dominant cost is not arithmetic, it is `np.add.at`** — 79–86 % of
`step()` on the atom-resolved path. A one-line CPU change (`np.bincount`) is
**bitwise identical** and ~28× faster on that operation, which is a ~4× chemistry
speedup available *before* any GPU work. The GPU port should be staged behind it.

---

## 1. Measured hot-loop breakdown

### 1a. Microbenchmark — inside `step()` (F = 1396 faces, the ml13 endpoint count)

Scratch script (not in repo). Times are per `step()` call, mean of 3.

| atoms | `step()` atom path | `step()` scalar path | 11× segment `add.at` | per-atom elementwise | face algebra ("rest") |
|---|---|---|---|---|---|
| 10 000 | 5.9 ms | 0.56 ms | **4.5 ms (76 %)** | 0.6 ms (10 %) | 0.8 ms (14 %) |
| 100 000 | 52.6 ms | 0.58 ms | **45.0 ms (86 %)** | 7.3 ms (14 %) | 0.4 ms (1 %) |
| 1 000 000 | 560.1 ms | 0.63 ms | **443.2 ms (79 %)** | 95.7 ms (17 %) | 21.3 ms (4 %) |

Cost is ~0.56 µs per atom and **linear in atom count**; the scalar (non-atom)
path is flat at ~0.6 ms, i.e. the face-resolved reservoir algebra is already
negligible. The atom path costs ~900× the scalar path at 10⁶ atoms.

Single-operation comparison, same inputs:

| operation | 10 k | 100 k | 1 M |
|---|---|---|---|
| `np.add.at(out, face, v)` | 0.40 ms | 4.09 ms | 40.28 ms |
| `np.bincount(face, weights=v, minlength=F)` | 0.02 ms | 0.14 ms | 1.43 ms |
| speedup | 26.1× | 28.8× | 28.2× |

`np.add.at` is numpy's unbuffered generic scatter and is known-slow; `bincount`
is a specialized C loop. **Both accumulate in input-array order**, so they agree
bitwise — verified below.

### 1b. Instrumented production run — real workload sizes

`mixed_layer.step` and `MixedLayerMechanism.advance` wrapped with counters, exact
production flags, first engine steps:

| advance # | material | faces | **atoms** | sub-steps | chemistry share of wall |
|---|---|---|---|---|---|
| 1 | mask (m1) | 32 | 3 719 104 | 1.0 | 6 % |
| 2 | oxide (m2) | 692 | 8 855 184 | 1.0 | 13 % |
| 3 | mask | 36 | 4 632 708 | 1.0 | 10 % |
| 4 | oxide | 700 | 5 029 782 | 1.0 | 12 % |

Two facts that dominate the design:

1. **Atom counts are 3.7–8.9 M per call**, an order of magnitude above the
   microbenchmark's top row. Atoms-per-face is enormous: **116 k on the 32-face
   mask, 12.8 k on the 692-face oxide.** Long segments are good news for a CSR
   design (plenty of work per segment) and bad news for naive atomics
   (contention on very few output slots).
2. **Sub-steps = 1.0 only because this is the dt-collapsed startup transient.**
   `advance()` uses `n_steps = ceil(duration_s / default_max_step_s)` with
   `default_max_step_s = 0.01` (mixed_layer_mechanism.py:106). At the production
   dt observed in ml13/ml14 (0.14–0.20 s) that is **17–20 sub-steps per engine
   step**, each of which currently redoes the entire per-atom block. The 6–13 %
   chemistry share measured here therefore **understates production by ~17×**;
   chemistry, not transport, is the production bottleneck.

Caveat to close on a box: the ml13/ml14 GPU-box logs show 13–22 s total wall per
engine step at production dt, which is not consistent with 17 sub-steps × the
per-atom cost measured here — either the box host is much faster on this loop or
production atom counts differ from the startup transient. **Stage 0's gate
includes re-running this instrument on the box** rather than extrapolating.

Also worth measuring: the choke-point y-strip symmetrization (64a71fa) replicates
each landed event across strip-equivalent faces, so it multiplies atom count by
the strip group size (≈3 for the quasi-2D slab). That is a ~3× input to every
number above and was introduced this week.

### 1c. Loop-invariant work inside the sub-step loop

`atom_face / atom_flux / atom_energy / atom_cos` are **constant across the
sub-step loop**; only the film thickness `d_fc` changes. Auditing the eleven
segment sums (mixed_layer.py:264–285) against that:

| kernel | depends on | invariant across sub-steps? |
|---|---|---|
| `kernel_sputter` | `atom_energy`, `atom_kress` | **yes** |
| `kernel_dexl` | `atom_energy` | **yes** |
| `kernel_act` | `atom_flux` | **yes** |
| `total_atom_flux` | `atom_flux` | **yes** |
| `kernel_mix`, `kernel_xl`, `kernel_complex`, `kernel_bare`, `kernel_ac`, `e_iface`, `eps_dep` | `atom_e_iface` ← `d_fc` | no |

**4 of 11 segment sums (36 %) are recomputed identically 17–20 times per engine
step.** Additionally, inside `interface_energy_eV`, the stopping-table lookup
`lam(atom_energy)` is invariant — only `exp(-d_fc[face]/lam)` varies — so most of
the per-atom `log`/interpolation work is hoistable too, leaving one `exp` per
atom per sub-step plus the `_deposited_energy` table lookups on `e_iface`.

---

## 2. Warp kernel design for the per-atom chemistry

### 2a. Precedent in this codebase

- `threed.py:210–224` `_edge_adjacency`: the existing sort-based pattern —
  int64 key encode → `cupy.argsort` on GPU (0.26 ms vs 66 ms host) → segment
  boundaries via `key_s[:-1] == key_s[1:]`. **The atom→face CSR build reuses this
  verbatim** with `key = atom_face`.
- `threed.py:258–263` `_smooth_scatter`: the atomics precedent (`wp.atomic_add`,
  float32, order-nondeterministic, accepted there because smoothing tolerates
  ~1e-5).
- `threed.py:1060` already uses `np.bincount` for the sky-view sum — the
  bincount idiom is established in this codebase.

### 2b. Proposed structure

Split the work by what changes:

**Prep — once per engine step (per material):**
1. `order = argsort(atom_face)` (stable) → CSR `offsets[F+1]`, permuted atom
   arrays. cupy on GPU, `np.argsort(kind='stable')` on host — same code shape as
   `_edge_adjacency`.
2. Per-atom invariants computed once: `kress`, `tpy(E; 0.9,20,500,0.5)`,
   `tpy(E; 0.3,8,500,0.5)`, `lam(E)`, `flux`.
3. The 4 invariant segment sums, once.

**Per sub-step (×17–20), device-resident, no host round-trip:**
- **Kernel A (map, one thread per atom):** `e_iface = E·exp(−d_fc[face]/lam)`;
  `eps = table(e_iface, cos)`; `tpy_bare`, `tpy_ac`. Writes the 7 state-dependent
  per-atom channels.
- **Kernel B (segmented reduction):** the 7 sums over CSR segments.
- **Kernel C (per-face):** existing reservoir algebra (4 % — port last or leave
  on host once A/B dominate).

**The load-bearing optimization is fusing A and B into a single pass**: one
thread walks its segment, computing all 7 channels into registers and
accumulating 7 running sums. This reads each atom's inputs **once** instead of
7×, converting a bandwidth-bound problem into a compute-bound one — and this
workload's segments (12.8 k–116 k atoms) are long enough to amortize everything.

### 2c. Reduction variants, ranked

| variant | parallelism | determinism | notes |
|---|---|---|---|
| **T1** CSR, one thread per face, sequential in original atom order | F = 32–700 threads | **bitwise identical to CPU** | ideal correctness reference; poor occupancy on a 3090 (82 SMs) |
| **T2** CSR, fixed-size chunks (e.g. 128 atoms), partials in chunk order, second fixed-order pass | ~A/128 ≈ 30–70 k threads | deterministic run-to-run; ≠ CPU bitwise (different association, ~1e-15 rel) | the performance path |
| **T3** raw `wp.atomic_add` | A threads | **non-deterministic** | rejected for the gated path; contention would be severe anyway (32–700 output slots for 4–9 M adds) |

Recommendation: **implement T1 first** (it preserves today's exact-equality gate
unchanged), then add T2 behind a flag once its own gate exists. T1's low
occupancy is partly offset by segments being huge — each of the 32 mask threads
does 116 k iterations of real work — but T2 is where the speed is.

---

## 3. Determinism contract

Today's contract is the strictest possible: `test_vectorized_step_matches_scalar_bitwise`
asserts **exact float equality** between N scalar calls and one vectorized call,
and the ledgers close to < 1e-9.

What each stage can honor:

- **Stage 0 (`bincount`)** — bitwise identical, *verified*: 5 trials × 300 000
  atoms into 1396 faces, `np.array_equal(add.at, bincount) == True`, max abs
  difference **0.0**. Both accumulate in input order. **No gate changes.**
- **Stage 1 (hoisting invariants)** — bitwise identical by construction
  (identical expressions, computed once instead of 17 times). **No gate changes.**
- **Stage 2 GPU/T1** — bitwise identical *if* FP64 is used and the segment walk
  follows original atom order; `exp`/`pow` must then also match libm bit-for-bit,
  which is **not guaranteed across CUDA vs host libm**. Realistic contract:
  **run-to-run bitwise reproducible on-device**, and vs CPU **≤ 1e-14 relative**
  on every state field, with ledger closure still < 1e-9.
- **Stage 2 GPU/T2** — same as T1 but ≤ 1e-13 relative vs CPU.

Replacement gate (proposed, to be signed off before Stage 2 lands, per the
receipts doctrine):
1. `test_gpu_chemistry_reproducible`: identical inputs → bitwise-identical
   outputs across 3 launches (catches atomics sneaking in).
2. `test_gpu_matches_cpu_within_tolerance`: all 9 state fields + etch/growth
   velocities within the stated relative tolerance over a mixed-condition
   ensemble.
3. `test_gpu_ledgers_close`: the existing < 1e-9 conservation gate, unchanged,
   run on device output.
4. The existing scalar==vector bitwise gate **stays** and continues to run on the
   CPU path, which remains the reference implementation.

---

## 4. Staged plan

| stage | change | gate | est. speedup (chemistry block) |
|---|---|---|---|
| **0** | `np.add.at` → `np.bincount` in `_segment` (CPU) | existing suite unchanged + explicit bitwise-equality test add.at vs bincount | **~4.2×** (560 → ~133 ms at 1 M atoms) |
| **0b** | re-run the instrument on the GPU box at production dt | records true sub-step count + atom counts + chemistry share; no code | — (measurement) |
| **1** | hoist loop-invariant work out of the sub-step loop (4 segment sums, `kress`, 2 `tpy`, `lam`); cache on a prep object, mirroring `threed.py`'s `prep` pattern | bitwise unchanged suite | ~1.3–1.5× on top |
| **2a** | box spike: warp FP64 arrays + `atomic_add` + `exp`/`pow` throughput on sm_86 | a throughput number, go/no-go on FP64 | — (risk retirement) |
| **2b** | CSR build (cupy argsort, `_edge_adjacency` pattern) + fused map/reduce kernel, variant T1, FP64, device-resident across sub-steps | new gates 1–3 above; CPU path untouched | 10–30× on the block (bandwidth-bound estimate: 7 channels × 8 B × 9 M atoms ≈ 0.5 GB/sub-step at ~900 GB/s) |
| **2c** | variant T2 (chunked) behind a flag if T1 occupancy proves inadequate | gates 1–3 at the looser tolerance | additional 2–5× |
| **3** | port Kernel C (per-face algebra) + keep state device-resident for the whole engine step, one readback | full-pilot A/B: identical profile to CPU within gate tolerance | removes the last per-sub-step transfer |

Sequencing note: **Stages 0 and 1 are pure CPU wins with zero gate risk and no
GPU box required** — they should land before any kernel is written, both because
they are free and because they shrink the problem the GPU port has to solve.

---

## 5. Risks

1. **FP64 on consumer GPUs (highest).** Every existing warp kernel in this repo
   uses `wp.array(dtype=float)` = **float32**, with host-side casts back to
   float64 (`threed.py:178, 380, 1444`). Transport tolerates that (flux smoothing
   matches to ~1e-5); **the chemistry cannot** — its gates are bitwise equality
   and 1e-9 ledger closure. On sm_86 (RTX 3090) FP64 arithmetic runs at 1/64 of
   FP32. Mitigating facts: the fused kernel is bandwidth-bound (FP64 costs ~2×
   bytes, not 64× time), and Stage 1 removes most `pow`/`log` from the loop,
   leaving one `exp` per atom per sub-step. Stage 2a exists specifically to
   measure this before committing.
2. **Warp FP64 API surface.** `wp.float64` arrays and `wp.atomic_add` on doubles
   need sm_60+ (satisfied) and warp support at the installed version (1.14 local
   / 1.15 on the box) — verify in 2a, not by assumption.
3. **Segment-length skew.** 116 k atoms/face on the 32-face mask vs 12.8 k on the
   oxide: T1 gives 32 threads for the mask material. If 2a shows this is
   unacceptable, T2 becomes mandatory rather than optional.
4. **Atom-count growth from symmetrization.** The choke-point y-strip
   redistribution multiplies events by strip-group size (~3×). Worth confirming
   whether the sums can instead be computed on one strip and broadcast — that
   would be a ~3× *algorithmic* win available on CPU, and it is arguably more
   faithful to the y-invariance contract than replicating events. **Flagged for a
   separate look; outside this scoping pass.**
5. **Determinism-contract change requires sign-off.** Moving off bitwise equality
   for the device path is a doctrine change; Stage 2 must not land until the
   replacement gate is agreed.
6. **Scope discipline.** Stage 0 alone is a ~4× chemistry win, bitwise safe, and
   one line. If GPU work is deprioritized, that stage still deserves to land.

---

## 6. Second-pass measurements (independent run, same day)

Added by a second scoping pass. Confirms §1–§5 and closes the open caveat in §1b.

### 6a. Stage 0 measured end-to-end, not extrapolated

§4 estimates Stage 0 at ~4.2× from the single-operation ratio. Measured directly
by patching `_segment` to `np.bincount` and timing the **whole `step()`** at the
real oxide workload (4 500 000 atoms, 704 faces, mean of 3):

```
step() as-is           : 2464.7 ms
step() with bincount    :  514.8 ms      ->  4.8x
bitwise identical state : True   (all 8 reservoir fields)
bitwise identical rates : True   (sif4_rate)
```

So the bitwise claim holds for the **full step result**, not only the isolated
scatter, and the speedup is 4.8× rather than 4.2×.

### 6b. The §1b caveat is resolved — the cost model reproduces the box

§1b flags that 17 sub-steps × the measured per-atom cost looked inconsistent with
box wall times and defers to a re-measurement. Projecting the number above:

```
2464.7 ms  x  17 sub-steps  x  2 materials  =  83.8 s per engine step
```

against the **~80 s/step** reported for ml13/ml14 with the cascade active. The
model is consistent; the chemistry block *is* essentially the whole production
step. Stage 0 moves that to **17.5 s**, and Stage 1 to ~14.5 s.

(The independent instrument in this pass measured 2.9–5.0 M atoms on 36/704 faces
versus §1b's 3.7–8.9 M on 32/692 — same order, same conclusion; atom count drifts
with surface area as the mesh evolves.)

### 6c. The profile inverts after Stage 0 — the GPU target is *not* the scatter

Re-profiling the components at 4.5 M atoms **after** the bincount swap (total
~515 ms):

| component | time | share of post-Stage-0 `step()` | invariance |
|---|---|---|---|
| `_deposited_energy` (2 table lookups + interp) | 175.5 ms | 34 % | state-dependent |
| `interface_energy_eV` (table lookup + `exp`) | 116.2 ms | 23 % | `lam` lookup invariant |
| 4 × `_threshold_power_yield` (`pow`) | 84.5 ms | 16 % | 2 of 4 invariant |
| **11 × `bincount` (the reduction)** | **74.9 ms** | **15 %** | 4 of 11 invariant |
| kress angular factor | 14.9 ms | 3 % | invariant |
| face algebra / clamps / ledgers | ~49 ms | 9 % | O(F) |

**Consequence for §2:** once Stage 0 lands, the segmented reduction is only 15 %
of the block and the port's real target is a **76 %-elementwise transcendental map
over ~4.5 M atoms — which needs no atomics at all.** That strengthens the fused
map/reduce recommendation in §2b (the map is the expensive half, so fusing it with
the reduction avoids a second pass over the largest arrays), and it *sharpens* risk
#1 rather than relaxing it: the work that lands on the GPU is dominated by `exp`,
`log` and `pow`, which on sm_86 have no FP64 SFU path and are software-emulated.
The bandwidth-bound consolation in risk #1 applies to the reduction, not to the
map. **Stage 2a's go/no-go measurement should therefore time the elementwise
kernel specifically, in FP64, before any CSR work is written.**

Hoistable-work total (§1c) measured at this scale: ~57 ms elementwise + ~27 ms
segment = **~84 ms of 515 ms (16 %) recomputed identically every sub-step**,
confirming Stage 1's ~1.2–1.5× estimate.
