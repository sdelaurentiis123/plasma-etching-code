# Resona design-partnership call — outcomes, challenges, prep plan

Date: 2026-07-20. **LOCAL ONLY — never push this file; it is a standing pre-push scrub item
(with `RESONA_PATTERN_TRANSFER_PARTNERSHIP_2026-07-17.md`).**

## Outcome

Design partnership accepted in spirit on the call. Arun (Resona Semiconductor) is sending **two
STL challenge files (~1 week)** with write-ups of expected outputs; results wanted ~2 weeks after,
format = 2-D/3-D cross-sections he can eyeball ("you cross the threshold for usability or you
don't"). He offered clean-room access (Columbia leverage, plus Stanford/Caltech through him), a
possible residency, and named three end states verbatim: "funding it," "being a collaborator," or
"you could do this inside of our company." Quote to keep: "if you start generating really novel
data here, I wanna be a part of whatever you end up doing."

## The two challenges

1. **Memory (aspect ratio).** Target: **200:1** — his number for 3-D DRAM (~400 stacked layers,
   ~30 um deep at ~20 nm features). Calibration points he gave: Samsung cryo leading edge ~100:1;
   40:1 = "a very good etch tool for 2 nm node at Intel." Bar: "if you can get to 200 in
   principle, then I believe you get it in practice, and you've got something real."
2. **CMOS (fidelity).** Very small feature sizes; the metric is **line-edge roughness and edge
   damage** — he says the industry lever is chemical (chemisorption pathways), his core interest.

## Immediate commitments

- **Today:** send him the lab-page link (LIVE: standelaurentiis.com/lab/petch-evidence.html) plus
  the 2-D cross-section images he asked for (necking + near-top roughness).
- **Wednesday:** Stan flies to SF (marathon, through the weekend, back Monday). Visit the new
  Resona office during business hours, meet the team. Coordinate by email.
- **~1 week:** STL files arrive. **~2 weeks:** results call.

## His technical tells (shape the roadmap around these)

- Wants to build a **"Resona etcher"** after the litho tool — simpler physics, not "layers of
  complexity on existing methodologies." Etch is "probably more disruptive than lithography."
- Actively designing a **metastable-atom (helium) plasma etch** concept: Lace Lithography's
  physics (metastable atoms; he thinks they fail as litho) repurposed as a precise, low-damage,
  FIB-like / ALE-adjacent directed-beam etch. He is parameterizing the metastable He source now.
- Believes tool companies win via "proprietary dataset → custom model → powers our tools";
  bullish on an **AI-native etch company** ("I think Lam sucks"); likes focused Periodic-style
  data plays (cites Radical for alloys), dislikes unfocused ones.
- On models: predicted mean-field/MC will miss hardware non-idealities — expects layering
  hardware-reality corrections on the MC core (this matches the calibrate-per-chamber contract).

## Fact-check of what we told him (correct in follow-up, casually)

- **Necking is NOT charging** in the Krüger case — it is mask-film growth vs ion sputter; our own
  causal audit found charging negligible there. Charging drives *notching* (demonstrated) and is
  implicated in twisting (not yet demonstrated by us). Correct this in the email.
- We score against the **experiment** (SEMs), not against Kushner's MCFPM — the stronger claim.
- de Boer Fig 9 is development data under the current program (re-earning); Krüger sealed reveal
  is the clean held-out story.
- AMR: built, measured, **no-go'd by its own gate** (cost is chemistry/transport-dominated;
  1.415x memory vs 3x bar) — say that, not "in flight."
- Speed framing: we built a separate engine that beats ViennaPS head-to-head (8-15x depth-matched,
  ~50-100x on specific cases); we did not CUDA-patch ViennaPS itself.

## Prep plan for the STL challenges (start before files arrive)

1. **STL -> level-set importer** (engine consumes level sets; voxelize + signed distance +
   `FeatureGeometry3D`). First thing that breaks otherwise. Small, do first.
2. **200:1 feasibility program**: free-molecular transport + charging at extreme AR through the
   common engine; expect charging + sub-degree IADF to dominate (the de Boer >20:1 frontier).
   Deliverable framing: "physics-faithful 200:1 exploration with stated uncertainties," never
   "validated at 200:1" (nobody is). Cryo SF6/O2 chemistry is the natural regime (de Boer kernel
   exists). Wall-clock check at 200:1 grids needed early.
3. **LER/roughness modality — the genuine capability gap.** We explicitly do not model physical
   LER today (profile jaggedness = numerics). Minimum credible: propagate a prescribed mask-edge
   roughness spectrum through 3-D etch (needs true 3-D lateral extent, not quasi-2D cells) and
   report post-etch LER transfer. Stochastic chemistry LER is research-tier; scope honestly.
4. **Deliverable format**: 2-D cross-sections + 3-D renders (movie pipeline already built);
   annotate necking/roughness features; publish as lab pages like the evidence pack.
5. Keep the **Kruger blind reveal** moving — a completed held-out validation upgrades every
   conversation above.

## Standing discipline

Resona specifics stay out of pushed branches and the shared docs repo. This file and the
2026-07-17 partnership doc are local-only. The public lab page is partner-neutral by design.
