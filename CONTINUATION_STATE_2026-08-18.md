# Multiphysics continuation state — 2026-08-18

Branch: `codex/validation-first-multiphysics`.

This is the current authoritative handoff. All campaign workers are destroyed.
The unrelated untracked mixed-layer log and `mouth_equilibrium_probe_dx/`
directory predated this work and remain untouched.

## Bottom line

The deterministic reactor-to-feature stack is materially stronger, but neither
the Freddie/Oxford final SEM nor Krüger's `825 nm` depth is yet an absolute,
validated prediction.

- The Oxford reactor boundary now includes conserved daughter chemistry,
  electron kinetics, collisional/sheath sensitivity, radial wafer transport,
  and a species-resolved dose envelope.
- The complete blind conditional square-pillar board contains 56 profile
  endpoints across seven widths, two rate endpoints, two ion-energy cases, and
  two angular-tail cases. It is rendered in one common-scale SVG atlas.
- A preregistered `5 nm` refinement closes depth, top/middle CD, timestep,
  CPU/CUDA replay, symmetry, particle balance, and state-remap checks. The
  literal bottom-CD gate fails (`8.612 nm` change versus `5 nm` allowed), so
  bottom CD is not certified.
- That failure is localized to the sidewall junction with the still-unetched
  film. The body through `75%` relief changes by at most `2.506 nm`; the `5 nm`
  threshold is first exceeded at `85%`. The sentinel leaves `439.203 nm` of
  continuous TiO2 below the floor, so it has no independent pillar bottom.
- The wider lower section is not nonphysical growth. The conditional law
  permits only nonnegative TiO2 removal, pins Cr/fused silica, and disables
  deposition/redeposition. The exact pre-clear junction flare remains
  grid-unresolved and is not a target prediction.
- Independent TiO2 measurements now force the missing target surface topology:
  ion-energy-dependent removal, TiOxFy/TiFx surface state, neutral supply,
  physical passivation inventory, ion-assisted passivation removal, Cr-mask
  motion, and pattern-dependent transport.
- The common conservative surface kernel is now exposed as a material-labelled
  `ReducedFluorinatedOxideMechanism`. It is ready for an evidence-bearing TiO2
  deck; tests prevent a relabelled SiO2 deck from becoming predictive. An
  optional evidence-bearing passivation bulk density converts its conserved
  inventory into real outward/inward level-set velocity; absent that density,
  no thickness is invented and legacy behavior is unchanged.

No claim of atomic accuracy is supported. The present best grid rung is `5 nm`,
and the Oxford TiO2/Cr coefficients have not been measured for the supplied
condition.

## Freddie / Oxford NPG80 condition

Frozen supplied process:

- Oxford PlasmaPro NPG80 RIE
- `55/5/1 sccm CHF3/SF6/O2`
- `30 mTorr`, `150 W` forward RF, `20 C`, `1200 s`
- `700 nm` ALD TiO2, `45 nm` Cr, fused silica
- conditional geometry prior: square pillars, `400 nm` pitch, widths
  `80--320 nm` (not target GDS confirmation)

The versioned profile atlas is
`results/curated/zhu_npg80_conditional_profiles_v1/profile_atlas.svg`.
It shows all 56 endpoints without using a target SEM or target depth. Fifty
conditional endpoints clear; six wide/low-rate cases do not. Cleared panels
show the last pre-clear TiO2 geometry because the one-moving-material rung
does not invent a post-clear fused-silica profile.

The reactor predicts a central positive-ion flux near
`1.03--1.06e19 m^-2 s^-1` over the declared boundary ensemble. This is a
physics-constrained boundary estimate, not a direct Oxford measurement. The
achieved self-bias/electrode waveform remains unresolved by `150 W` forward
power alone.

### Surface-science advance

Three independent response boards now constrain any TiO2 mechanism:

1. Choi 2013: TiO2 rate rises `130.9 -> 197.2 nm/min` as DC-bias magnitude
   rises `50 -> 250 V`; XPS and the authors' mechanism require fluorination,
   bond breaking, and ion-assisted product desorption.
2. Ji 2024 RF board: interfeature gap narrows `95.96 -> 18.02 nm` as RF power
   rises `90 -> 210 W`. A removal-only solid cannot reproduce that response;
   positive retained/deposited surface volume and mask/passivation evolution
   are required.
3. Ji 2024 spacing board: the `100` and `70 nm` gap points separate from the
   `350--750 nm` cluster. The paper's verbal `100 nm` boundary conflicts with
   the already-shifted `100 nm` datum, so the empirical transition is only
   bracketed between `100` and `350 nm`.

These experiments establish topology and response signs. They change tool,
feed, film, mask, and/or boundary, so none supplies an Oxford coefficient.

### What still blocks a blind absolute Oxford SEM prediction

Required same-run data, in descending value:

1. achieved DC self-bias or the electrode voltage/current waveform,
2. blanket TiO2 loss after the exact `1200 s` recipe,
3. remaining Cr thickness (or Cr loss),
4. actual GDS dimensions, sample radius/orientation, and ALD film phase/density,
5. cross-section and top-down SEMs as the frozen answer key.

The SEM is not needed to run a blind prediction; it is needed to score it.
Items 1--4 determine whether the absolute prediction is unique and
mechanistically attributable.

## Krüger depth

Krüger remains an honest miss: published-boundary simulation is about
`346.833 nm` versus the reported `825 nm`. Do not revive the retracted
"blanket-anchor" story: `13.75 nm/s` is a feature-average rate, not a published
blanket measurement.

This continuation added four C4F6 boundary advances:

- source-backed direct-ion fragmentation from NIST C4F6 spectra plus CFx
  secondary ionization,
- Benck neutral-radical ratio and voltage-response boards,
- a pixel-audited light-ion inverse that rejects one common loss rate,
- a differential-loss closure whose source remains nonnegative over all 15
  Ar-mixture electron-temperature rows.

The Ar-mixture board supports selective CF3+ loss/common loss of roughly
`0.1025--0.6135`; pure C4F6 changes sign, so one fixed operator across both
regimes is still rejected. This narrows missing reactor physics but does not
supply Krüger's absolute species-resolved wafer boundary. Exact depth still
requires the authors' HPEM/PCMCM species flux/IED output or an independently
validated C4F6 reactor boundary.

## Validation and reproducibility

Relevant exact receipts:

- `results/curated/zhu_npg80_conditional_profiles_v1/`
- `results/curated/zhu_npg80_profile_convergence_refinement_v1/`
- `results/curated/zhu_npg80_tio2_surface_topology_v1/`
- `data/experimental/ji_2024_tio2_hierarchical/`
- `data/experimental/choi_2013_tio2_cf4/`
- `results/curated/benck_c4f6_differential_loss_v1/`

The literature library now indexes 170 sources with reverse
constant/law/decision links. Fetch-to-extract-to-entry and bibkey provenance
remain binding.

Final all-files suite result after physical passivation growth and legacy
provenance compatibility: **2101 passed, 1 skipped in 1209.06 s**. This was run
at commit `3ee3146` after the exact C4F6 topology-replay regression was repaired.
The regression was provenance-only: a disabled optional density had appeared as
a new `null` field in a historical receipt. The default now reproduces the old
receipt exactly, while an explicitly supplied density remains recorded and
drives physical growth.

## Correct next build order

1. Land a TiO2-specific parameter deck on the generic fluorinated-oxide
   contract. Keep every missing Oxford coefficient nonpredictive/fail-closed.
2. Constrain the now-implemented physical passivation-density conversion and
   its neutral-radical supply with target-relevant evidence; preserve analytic
   bounded state updates.
3. Evolve Cr as a second moving material and validate mass/state remap.
4. Grade the mechanism first against Choi energy response and Ji RF/spacing
   boards, without using Freddie's target SEM.
5. Re-run the frozen Oxford width board with same-condition bias/blanket/Cr
   measurements when available, then reveal and score the SEM.
6. For Krüger, obtain or reconstruct the missing C4F6 species-resolved wafer
   boundary; do not tune a surface yield to `825 nm`.

Standing rules: deterministic/differentiable operators; no Monte Carlo in the
production path; zero undeclared knobs; one physics change per validation;
source and uncertainty receipts for every coefficient; no target-selected
coefficient; destroy campaign workers after results are copied.
