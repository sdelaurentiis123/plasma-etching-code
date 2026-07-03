# Chemistry expansion: from one chemistry to the industry workhorses

Deep plan, literature-first (researched 2026-07-03). Motivation: the #1 documented honest limit
(what-petch-does.html, RECONCILIATION.md) is **"one wafer dataset, one chemistry"** — petch's
surface chemistry is SF6/O2-on-Si only (Belen/Ertl), validated on de Boer cryo. This plan adds the
four industry-workhorse chemistries, each with a PUBLISHED surface model and a digitizable
validation dataset, in dependency order. Working rules unchanged: every claim gets a reproducer
script + numeric gate vs published data; no unpublished fudge knobs; honest fails documented.

Relation to other plans: independent of CHARGING_PHYSICS_PLAN.md (that is transport/charging; this
is surface chemistry) — they interleave freely, and W3 here has an optional gate that CONSUMES the
charging capability (charging-induced etch stop in HAR dielectric — a measured phenomenon).

## W0 — Chemistry plugin interface + real multi-material (architecture prerequisite)

**What.** Today `chemistry.py` hardcodes Belen SF6/O2 (theta_F/theta_O fixed point) and `threed.py`
knows materials only as mask / substrate / etch_stop_z. Refactor so a chemistry = a small module
implementing `{coverages(fluxes) fixed point, yield triplet(E,theta), rate law}`, selected per
MATERIAL ID; the grid carries a material-id field (PR / poly-Si / SiO2 / Si3N4 / Si) generalizing
`etch_stop_z` (a stop is just material with rate 0); each mesh face inherits the material of the
solid cell it borders. Also: allow SIGNED extended velocity (deposition = outward front motion) —
the upwind-Godunov advection already handles sign; today's `V = max(V-..., 0)` clamps growth, so
deposition needs a dedicated path + a mask-pinning decision (mask does NOT stay pinned under
deposition). Add a deposition-only unit test (uniform growth of a flat surface at known rate).

**Gates (W0).** (1) SF6/O2 regression: de Boer evolving gate (RMSE ≤ 0.05) and the notch gates
reproduce EXACTLY (same seeds → same numbers) through the plugin path. (2) Deposition test: flat
front advances outward at the imposed rate within 2%. **Effort: 2–4 days.** This is pure
refactor + plumbing; no new physics constants.

## W1 — O2 plasma etch of polymers (resist trim / ash / descum)

**Why first.** Cheapest new chemistry; PMMA is the SAME material and SAME research group as the
charging SEE anchor (Athens/Demokritos), so constants cross-check; and resist processing (trim to
sub-print CDs, descum, strip) is the most-demanded companion step to any pattern-transfer etch.

**Model (published).** Gogolides et al., "Oxygen plasma etching of hydrocarbon-like polymers:
Part I Modeling," Plasma Process. Polym. (2018): ion-enhanced etching DOMINATES (synergy of
adsorbed O + ion energy; physical sputtering / pure chemical / UV channels negligible at
~100 eV, 300–400 K) — structurally the same Langmuir-coverage + sqrt(E)-yield form as our Belen
model, one adsorbing species instead of two. Constants from the paper; nothing fitted.

**Validation gates.**
1. Blanket-rate gate: "Part II experimental validation" (Bès et al., PPP 2018) — measured PMMA
   (+PS/PET) etch rates vs ion flux/energy. Digitize; gate RMSE on normalized rates ≤ 0.15
   (blanket-rate tolerance; state it).
2. Profile gate: resist-line trim sim (line of PR on Si, isotropic-ish O flux + ion vertical) —
   monotone CD reduction, no artifacts, trim rate consistent with the blanket gate.
3. Cross-consistency: same PMMA surface as charging-SEE work (Dapor yields) — one material card.
**Effort: 2–4 days** after W0.

## W2 — Cl2 / poly-Si (the gate-etch chemistry; makes the notch story chemistry-consistent)

**Why.** Our charging gates (Hwang–Giapis) and the notch experiment ARE Cl2 HDP data — but the
notch sim currently runs SF6/O2 surface chemistry under a Cl2-derived charging table. Adding Cl2
closes that inconsistency and gives petch the classic anisotropic gate-etch chemistry.

**Model (published).** Chang & Sawin, JVST A **15**, 610 (1997) — beam study of ion-enhanced
poly-Si etching: Langmuir Cl chemisorption (yield saturates with neutral-to-ion flux ratio),
ion-energy dependence **Y ∝ (√E − √E_th), E_th ≈ 10 eV**, angular dependence from the companion
literature. Negligible spontaneous etch at room temperature (undoped poly). Same structural form
as Belen — the W0 plugin makes this a parameter card + one-species coverage solve.

**Validation gates.**
1. Yield gate: digitized Chang–Sawin yield-vs-energy and yield-vs-flux-ratio curves (RMSE ≤ 0.15
   on normalized yield).
2. Profile gate: Chang, Arnold & Sawin, JVST B **18**, 172 (2000) — "Ion-assisted etching and
   profile development of silicon in molecular chlorine": reproduce the reported profile trends
   (anisotropy vs neutral-to-ion ratio; undercut onset). Qualitative-to-±30% first-wiring gate.
3. Consistency gates (re-runs, no new fits): microtrenching mechanism vs Hoekstra–Kushner
   (JVST B 16, 2102 — Cl2, already our ion-reflection anchor); **notch gate re-run with Cl2
   chemistry** — gates A/B/C of scripts/notching_depth_gate.py must still pass (the mechanism is
   charging-driven, so they should; measuring this is the point).
4. Stretch (optional): HBr/O2 passivation card for selectivity/taper — public model basis is thin;
   scope only if a published parameter set surfaces. Do NOT invent constants.
**Effort: 3–5 days** after W0.

## W3 — Fluorocarbon etch of SiO2/Si3N4 (the dielectric workhorse; biggest lift, biggest moat)

**Why.** Dielectric etch (contacts, vias, hard-mask opens, 3D-NAND/DRAM HAR holes) is the largest
plasma-etch segment and the one where charging matters most (insulating everything). No open
feature-scale code carries a gated FC chemistry + charging together.

**Model (published).** Two-layer surface state: a steady-state fluorocarbon film of thickness
h_FC per face (NEW per-face state variable) mediating etch vs deposition:
- Standaert et al., JVST A **15**, 1881 (1997) (CHF3 ICP): thin steady-state FC film controls the
  SiO2 etch rate; dep→etch transition with bias power; extreme sensitivity (~1 nm film ↔
  ~400 nm/min rate swing) — the reason h_FC must be an explicit state, not folded into a yield.
- Standaert et al., JVST A **22**, 53 (2004): same framework across Si, SiO2, Si3N4, SiC —
  ion-induced DEFLUORINATION of the film is a major channel (the film is a fluorine SOURCE, not
  just an inhibitor); gives the multi-material parameter set.
- Reference feature-scale implementation: Huard & Kushner et al., "Plasma etching of high aspect
  ratio features in SiO2 using Ar/C4F8/O2," JVST A **37**, 031304 (2019) — OPEN PDF (cpseg) —
  integrated modeling to AR 80; adopt their surface-reaction structure where Standaert
  under-determines it.
Balance per face: d(h_FC)/dt = polymerizing-neutral dep − ion-induced consumption(E, flux);
substrate etch rate = f(h_FC, ion flux, E) with the defluorination channel; oxygen chemistry
(from O2 feed or the oxide itself) consumes the film — that is what differentiates SiO2 from Si
selectivity. Deposition uses the W0 signed-velocity path when h_FC exceeds the thin-film regime.

**Validation gates.**
1. Blanket gate: Standaert CHF3 SiO2 etch-rate-vs-bias curve including the dep→etch transition
   bias (RMSE ≤ 0.15 normalized; transition bias ±20%).
2. Selectivity gate: SiO2:Si and SiO2:Si3N4 rate ratios vs the 2004 multi-material data (±30%).
3. HAR gate: ARDE in high-AR holes vs the published Ar/C4F8/O2 HAR data (Huard figures; tapered
   sidewall + slowing floor). First-wiring ±30% on the normalized depth-vs-AR curve.
4. STRETCH (consumes CHARGING_PHYSICS_PLAN): charging-induced etch STOP in sub-100 nm HAR holes —
   measured phenomenon (dual-frequency CCP studies report etch stop at 60 nm width that ion energy
   cannot recover but ion flux can — the electron-shading signature). Reproducing the stop
   qualitatively with charging ON vs OFF would be a headline second only to the notch.
**Effort: 1–2 weeks.** The h_FC state + deposition coupling is the real work; do NOT start before
W0's deposition test passes.

## W4 — Atomic layer etching (Cl2-dose / Ar-ion cyclic; precision regime)

**Why.** ALE is where etch goes at sub-50 nm CDs; it is ion-energy-WINDOW engineering (petch
resolves IEDFs and energy-dependent yields natively); self-limiting cycles make petch's
differentiability directly useful (recipe = sequence optimization, a demo no other open tool can
run end-to-end).

**Model.** Reuses W2 entirely: phase A = Cl2 dose with saturating chemisorption (Langmuir, no
ions); phase B = Ar-ion pulse removing the chlorinated layer, yields from the W2 card evaluated at
the pulse energy; a cycle scheduler alternates flux boundary conditions (pure parameter schedule —
no new solver). The self-limiting behavior comes out of the coverage saturation + finite
chlorinated-layer inventory (cap the removable depth per cycle at the chlorinated thickness).

**Validation gates (all published numbers).**
1. EPC gate: Kanarik et al. (Lam) Si Cl2/Ar ALE — EPC ≈ 1.2 Å/cycle at 25 eV; saturation curves
   (EPC flat past ~5 s Cl2 dose, ~10 s Ar pulse). Gate: EPC within ±30%, saturation shape monotone
   → flat.
2. Window gate: ALE window ~15–20 eV for normal-incidence Ar+ (MD + reduced-order-model study,
   J. Phys. Chem. B 2025; Athavale & Economou 1996 is the original) — below: no etch; above:
   sputtering breaks self-limitation. Gate: EPC(E) reproduces the window edges ±5 eV.
3. Synergy gate: Kanarik synergy metric (full-cycle EPC vs sum of individual half-steps) > 90%.
4. Profile-level reference: Kushner-group "ALE of 3D structures in Si" JVST A **35**, 031306
   (2017), OPEN PDF — qualitative comparison of profile fidelity vs continuous etch.
5. Demo (differentiability headline): gradient-optimize (dose time, pulse energy, n_cycles) to a
   target depth with minimum over-etch — end-to-end autodiff through the cycle schedule.
**Effort: 3–5 days** after W2.

## W5 — Sub-100 nm grid stability (cross-cutting prerequisite for leading-edge CDs)

Documented frontier (RECONCILIATION.md): sub-micron grids destabilize sidewalls even without
charging (seen at W = 0.5 µm, dx = 0.03). None of W1–W4's *gates* need it (they validate at
published-experiment scales), but applying any of them at modern CDs does. Scope: characterize the
instability vs dx (advection CFL? reinit frequency? narrow-band width? MC flux noise at small
faces?), isolate with a charge-off straight-wall etch at W = 100 nm, dx = 5 nm; fix the dominant
term. Gate: straight-wall anisotropic etch at W = 100 nm holds the wall to < 5% of W over a full
etch-through. Unbounded risk — timebox to 3 days of diagnosis before committing to a fix; if the
dominant term is MC noise, the deterministic transports (radiosity/knudsen) may already be the
answer — test those first.

## Order, interleaving, kill criteria

**W0 → W1 → W2 → W4, then W3; W5 timeboxed whenever blocked elsewhere.** W1 is the fast win and
validates the plugin interface on a second chemistry; W2 unlocks both the notch-consistency re-run
and W4; W3 is the biggest and lands last. CHARGING_PHYSICS_PLAN interleaves freely (different
files); its W1 (integrator) before this plan's W3 stretch gate.

**Kill criteria.**
- W0: if SF6/O2 regression cannot be made bit-exact through the plugin, keep the legacy path and
  make plugins additive-only (never break the validated baseline for architecture's sake).
- W1/W2: if digitized-data gates miss at > 2× tolerance with published constants, publish the gap
  (as with DDA) — do not tune. Check units/conventions first (the rate_scale minutes-vs-seconds
  slip cost a day once already).
- W3: if the h_FC balance cannot reproduce the dep→etch transition bias ±20% with Standaert's
  constants, bisect: film model vs flux inputs (their ICP fluxes are reported; use them, not ours).
  The HAR stretch gate is OPTIONAL — do not let it block landing the blanket+selectivity gates.
- W4: if the ALE window edges miss by > 5 eV, the likely cause is our IEDF width vs their beam
  monochromaticity — document, gate on beam-like narrow IEDF, keep the plasma-IEDF case as info.
- W5: 3-day diagnosis timebox before any fix work.

**What this buys.** Four chemistries × published gates turns "one wafer dataset, one chemistry"
into "the industry-workhorse set, each gated" — with two capabilities no open tool combines:
charging+dielectric HAR (W3 stretch) and differentiable ALE recipe optimization (W4 demo).

## Sources

- Gogolides et al., "O2 plasma etching of hydrocarbon-like polymers: Part I Modeling," Plasma
  Process. Polym. 15, 1800046 (2018); Bès et al., "Part II experimental validation," PPP (2018).
  https://onlinelibrary.wiley.com/doi/abs/10.1002/ppap.201800037
- Chang & Sawin, JVST A 15, 610 (1997) — Cl/Cl2/Cl+ beam study, Y ∝ (√E−√10 eV), Langmuir coverage.
  https://pubs.aip.org/avs/jva/article-abstract/15/3/610/961780
- Chang, Arnold & Sawin, JVST B 18, 172 (2000) — Cl2 profile development.
  https://pubs.aip.org/avs/jvb/article-abstract/18/1/172/470387
- Standaert et al., JVST A 15, 1881 (1997) — steady-state FC film, CHF3/SiO2.
  https://pubs.aip.org/avs/jva/article/15/4/1881/98733
- Standaert et al., JVST A 22, 53 (2004) — FC film role across Si/SiO2/Si3N4/SiC; defluorination.
  https://pubs.aip.org/avs/jva/article-abstract/22/1/53/242881
- Huard, Kushner et al., JVST A 37, 031304 (2019) — HAR SiO2 in Ar/C4F8/O2, integrated modeling
  (open PDF). https://cpseg.eecs.umich.edu/pub/articles/JVSTA_37_031304_2019.pdf
- HAR SiO2 dual-frequency CCP study (charging-induced etch stop at 60 nm; flux recovers it, energy
  does not). https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10222222/
- Kanarik et al. (Lam) — Si Cl2/Ar ALE: EPC 1.2 Å/cycle @ 25 eV, saturation, synergy metric
  (JVST A 33, 020802 (2015) review + follow-ups).
- Si–Cl2–Ar+ ALE window 15–20 eV: MD + reduced-order model, J. Phys. Chem. B (2025).
  https://pubs.acs.org/doi/abs/10.1021/acs.jpcb.5c01378 (OSTI open copy:
  https://www.osti.gov/servlets/purl/2586627)
- Kushner group, "ALE of 3D structures in silicon: self-limiting and nonideal reactions," JVST A
  35, 031306 (2017) (open PDF). https://cpseg.eecs.umich.edu/pub/articles/JVSTA_35_031306_2017.pdf
- Athavale & Economou, JVST B 14, 3702 (1996) — original ALE simulation.
