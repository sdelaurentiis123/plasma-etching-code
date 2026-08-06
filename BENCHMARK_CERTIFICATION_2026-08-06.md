# Benchmark certification under the final mechanism (2026-08-06)

One table, every benchmark the engine claims, graded against the mechanism as
it stands at this commit. Honest MISSes are recorded as MISSes. Every row
names the mechanism version it was measured under, because the campaign's
central lesson is that a number measured under a superseded mechanism is not
evidence about the current one.

**Mechanism version at certification** — the constants and laws in force:

| element | value | provenance |
|---|---|---|
| SiO2 physical sputter | `Y = 0.0139(sqrt(E) - sqrt(18))` | Gray MIT thesis 1993, Table 5-1 |
| SiO2 ion-enhanced yield | `beta_e = 0.053(sqrt(E) - sqrt(4))` | Gray, Eq. (5-35) |
| F adsorption coefficient | `s0 = 0.02` | Gray, Table 5-10 + p.246 (co-regressed with `beta_e`) |
| polymer angular class | Kress `B = 9.3` | Krueger Appendix B, class 1 |
| oxide/mask angular class | `B = 1.7` (peak 1.31) | bounded by Cho 2000 / Schaepkens 1998 |
| chemical angular class | unity to 45 deg, roll-off | Chang & Sawin 1997 (class 2) |
| crosslink formation | deposition-driven, per-material bond count | Krueger JVST A 42, 043008 (2024) sec. III |
| crosslink breaking | ion/hot-neutral/photon | Krueger thesis sec. 2.2.3 |
| thermalized radical return (E8) | implemented, **default off** | Huang thesis L5714-5727 |
| beam | two-component core+tail, tail fraction swept | Kim JJAP 64, 05SP15 (2025) |

Suite at certification: **1198 passed, 1 skipped**.

---

## (c) Transport exactness — the engine's strongest claim

| benchmark | target | measured | verdict |
|---|---|---|---|
| Clausing transmission, AR 1 | Santeler closed form, 0.5% | within 0.5% | **PASS** |
| Clausing transmission, AR 10 | Santeler, 0.5% | within 0.5% | **PASS** |
| Clausing transmission, AR 50 | Santeler, 1% | within 1% | **PASS** |
| Clausing transmission, AR 100 | Santeler, 1%; industry ~1.3% bottom flux | `tau = 0.013 +/- 0.0008` | **PASS** |
| AR 50 bottom delivery vs industry | Lam: "only 2.5% of the incoming flux reaches the other end" | 0.025287 | **PASS (1.1%)** |
| band convergence | quadrature-independent | converged | **PASS** |

Gate file: `tests/test_axisymmetric_exchange_3d.py` (6 gates). This is the only
benchmark in the engine graded against an *exact analytic solution* rather than
another simulation, and it holds at every aspect ratio to 200:1.

## (d) Hole study — phases 1 and 2

| benchmark | target | measured | verdict |
|---|---|---|---|
| phase-1 regression gates | committed characterization | all pass | **PASS** |
| phase-2 evolution driver gates | 15 gates incl. Clausing limit | all pass | **PASS** |
| degenerate wide-hole limit | 0-D blanket rate, 1% | rel 5.6e-5 | **PASS** |
| enclosure vs gated benchmark | exact | rel 0 / 2.8e-14 / 8.9e-14 (AR 10/50/100) | **PASS** |
| cascade vs phase-1 delivery | bitwise | rel 0.0 | **PASS** |
| per-step neutral balance | < 1e-9 | <= 3.2e-15 | **PASS** |
| coupled rate-ARDE at AR 200 | *no published reference exists* | rate falls 69% (tail 0.65) | **REPORTED, not gradeable** |
| rate-ARDE at AR 50 | measured band 43-80% (Nguyen & Jansen 2020 / Huang) | 60% | **inside band** |

The AR-200 coupled curve (tail 0.65): rate 6.62 -> 2.04 nm/s over AR 1 -> 200,
against energetic delivery falling only 31% — chemistry more than doubles the
ARDE that transport alone produces. The core-only beam saturates past AR 16
(2.94 -> 2.88 nm/s, AR 32 -> 100) and is flagged not-usable for attributing
depth trends to transport.

## (e) Beam / IADF gates

| benchmark | target | measured | verdict |
|---|---|---|---|
| two-component IADF gates | 13 gates | all pass | **PASS** |
| Kim 2025 measured widths | core/tail ratio 3.6x | 3.599 | **PASS** |
| Krueger Fig-4 round trip | digitization band [0.822, 0.860] | 0.8412 | **PASS** |
| analytic vs numerical acceptance | 1e-10 | <= 1e-10 | **PASS** |
| flux conservation | 1e-12 | 1e-12 | **PASS** |
| Gray N1 floor | 0.25 +/- 0.05 (dynamic range) | 0.2286 | **PASS** |
| Gray N1 plateau | ~1.10 absolute | 0.8807 | **MISS (0.80x)** |
| Gray N1 half-rise | Gray's closed form 100.2 at 350 eV | 175.8 | **MISS (1.76x)** |
| Gray N2 (CF2 reduces yield) | monotone decrease | decreases to clog | **PASS** |
| Gray six measured yield points | Table 5-10, 20-2000 eV | reproduced within fit residual | **PASS** |

Half-rise note: the previously quoted target of 27 +/- 8 came from a digitized
*replot* (Kwon Fig. 3.4) that the harness itself flags `[VERIFY]`. Gray's own
printed parameters predict 96-127 across the plausible energy assignment, so
the primary source refutes the replot. Graded against the primary source, the
pair transplant moved petch from 0.038x to 1.76x Gray's own model — a 26x
improvement — and the residual 1.76x is recorded as an open MISS, not smoothed.

## (b) SF6/O2 silicon arm

| benchmark | target | measured | verdict |
|---|---|---|---|
| de Boer direct-validation gates | committed suite gates | all pass | **PASS** |
| provenance grade | — | L3 (profile-fitted, per Belen 2005's own abstract) | **declared** |

The silicon arm's constants are profile-fitted by their source, so the de Boer
comparison is a *transfer* test, not a blind prediction. This was corrected in
the beam-constants atlas and is restated here so the certification does not
inherit an overclaim.

## (f) Krueger trench — the dossier gates, final state

| benchmark | target | measured | verdict |
|---|---|---|---|
| mask remaining | 850 nm | 850.2 | **PASS (exact)** |
| aperture at 270 nm depth | 41.1 nm (his sim) | 41.8 | **PASS (1.7%)** |
| clog at O2 0.5 | sealed | 0.000 sealed | **PASS (exact)** |
| necking absent at O2 2.5 | open | open | **PASS** |
| closure/etch, t >= 8 s | 0.0310 | 0.0257 / 0.0153 / 0.0160 | **in band** |
| mask constriction | 45 nm | 50.9 equilibrium | **+13%, bounded** |
| constriction depth | 200 / 271 nm | 170 | **-15%** |
| closure/etch, t = 1-4 s | 0.0310 | 1.83x (pre-multiplicity) | **above band** |
| depth, reference energies | 825 +/- 5% | 811 / 791 / 852 | **PASS** |
| depth, feature keV (pre-Gray) | 825 +/- 5% | ~1066 (+29%) | superseded |
| depth, feature keV (final) | 825 +/- 5% | ~418 projected (-49%) | **MISS** |

The depth MISS changed sign when the inherited two-row energy-law anomaly was
replaced with Gray's measured laws. It is decomposed, not merely attributed:
the physical channel read 1.22x too strong and the chemical channel 2.8x too
weak against Gray's absolute yields, and the residual now sits on fluorine
delivery to the front. Two declared-open supply items remain (the Kwon/Sawin
adsorption element beyond the scalar; the unpublished CFx+ fraction, swept and
shown immaterial). E8 was built, gated, and measured immaterial at this
geometry (+0.02% at AR 200 over the whole physical band).

## (a) Eight-condition scorecard

Run 2026-08-06 on box 46982614, all six conditions launched in parallel from
the same commit. **The 60 s endpoints were not reached and are not claimed.**
The box's warp build exposed no CUDA device, so the transport ran on CPU at
~2-4 min/step; a 60 s endpoint is ~4 h per condition on that hardware. What
follows is graded at a **matched simulated time of t = 1.948 s** -- the largest
time all six conditions had reached -- and the reduced duration is declared
rather than papered over.

Grader: `scripts/grade_scorecard_matched_time.py`, which refuses to grade the
endpoint-only criteria.

| condition | depth at t = 1.948 s | aperture | reached |
|---|---|---|---|
| base (6 kW, O2 1.0) | 13.547 nm | 83.717 nm | 1.96 s |
| O2 0.5 | 12.033 | 79.613 | 4.08 s |
| O2 1.5 | 12.822 | 83.115 | 2.95 s |
| O2 2.5 | 12.995 | 84.991 | 2.93 s |
| 4 kW | 12.217 | 83.688 | 3.52 s |
| 8 kW | 13.776 | 83.711 | 1.95 s |

**Power transfer -- both ratios in band, and this is the headline.**

| criterion | band | scorecard-1 (ml9a) | this mechanism | verdict |
|---|---|---|---|---|
| r(4/6) | [0.84, 0.94] | 0.672 (MISS) | **0.902** | **in band** |
| r(8/6) | [0.97, 1.06] | 1.085 (near-miss) | **1.017** | **in band** |

Scorecard-1 missed both power ratios and traced the misses to the over-narrowed
mouth modulating feature flux. Under the final mechanism both land inside their
published bands at matched time. Caveat stated plainly: these are ratios at
t = 1.9 s, not at 60 s, and the published bands are endpoint quantities. They
are reported as a strong early indication, not as a certified endpoint pass.

**Oxygen ordering -- correct at the extremes, one sub-cell inversion.**

  O2 0.5 -> 79.61 | O2 1.0 -> 83.72 | O2 1.5 -> 83.11 | O2 2.5 -> 84.99

The narrowest aperture is the starved 0.5 case and the widest is 2.5, which is
the published direction (oxygen thins the lip film). The 1.0/1.5 pair is
inverted by 0.60 nm, which is sub-cell at dx = 10 nm and is recorded as an
inversion rather than dismissed as noise.

**Not gradeable at this time, and not claimed:** the clog verdict at O2 0.5,
necking-absent at O2 2.5, and absolute depth against 825 nm are all endpoint
states. The trench dossier's own runs cover the first two (clog exact, necking
absent) under this mechanism.

**Cost to close properly:** one GPU-visible box (warp with a working CUDA
device), six conditions at 60 s, ~2-3 box-hours in parallel, well under $1.
The only reason it is open is the CPU-only box drawn this pass.
