# Krüger thesis Appendix B — verbatim constants (fetched 2026-07-27, curl -k)

Source: Krueger_Florian_PhD_Thesis_2024.pdf, cpseg.eecs.umich.edu (TLS cert
invalid → agents' fetch blocked; direct download succeeded). Columns:
p0, Eth(eV), n(=q), E0(eV), angular-form. (s)=fresh, (xs)=crosslinked,
X+ = ion, X# = hot neutral. All previously [VERIFY] items now resolved.

## Reflection / hot-neutral model (thesis §2.2.1-2.2.2, verbatim semantics)
- "All positive ions neutralize upon their first collision with a surface and
  return as a hot neutral that ... behave identically to their ion counterpart."
- Sputter reactions EMIT the continuing hot neutral in the product line
  (e.g. CF(s)+Ar+ → EP + Ar#): sputter AND continue — ADDITIVE, verbatim.
- Hot neutrals have identical reaction rows (multi-bounce explicit).
- Eq. 2.41 scattering: E > E_ts → purely specular, "preserve ALL the energy";
  below E_c or angle < theta_c → partially/fully diffusive with energy loss.
  (Numeric E_ts/E_c/theta_c values not located in the appendix table; likely
  in the MCFPM parameter set — remaining single [VERIFY].)

## Polymer (film) energetic rows
- CF/CF2/CF3(s) + ion/# → EP + #:      0.9  20  0.5  500   (confirms audit lift)
- CF(s) + ion → AC(s) + F + #:          0.01 20  0.5  500   (carbonization branch)
- CF/CF2/CF3(xs) + ion → EP + #:        0.6  50  0.5  500   (PC sputter resistance — was [VERIFY])
- CF/CF2/CF3(xs) + ion → CF(s) + #:     0.3   8  0.5  500   (ion DE-crosslinking — new mechanism)

## Deposition (neutral radicals)
- on fresh polymer (s):    CF/CF2/CF3 0.1, C2F3 0.03   (matches ours)
- on crosslinked (xs):     0.02                         (matches ours)
- on AC mask:              0.2   (appendix-converged; paper Table V optimized 0.0842 — reconcile)
- on bare SiO2 (complex):  CF/CF2 0.278, CF3/C2F3 0.2, heavies 0.001

## Oxygen
- polymer (s/xs) + O → EP: 0.0423  (appendix; paper optimized 0.0628 — reconcile)
- AC + O → CO:             1.0e-5   (mask essentially O-inert)

## AC mask energetic
- AC(s/xs) + ion/# → C + #: 0.001  200  0.4  250  (mask top nearly sputter-immune —
  explains experimental hm≈850; our mask-top over-erosion used far stronger removal)

## SiO2
- bare sputter: SiO2(s)+ion → SiO2 + #: 0.0852  70  1  140
- ACTIVATION: SiO2(s)+ion → SiO2*(s) + #: 0.9   (activated oxide state)
- complex formation on ACTIVATED oxide: CF 0.8, CF2 0.85, CF3 0.9
  (vs 0.2-0.278 unactivated → ion-activated chemisorption is the DOMINANT
  complex channel — mechanism we lack)

## Implications for petch mixed layer (priority)
1. Additive reflection + multi-bounce + full specular energy retention above
   threshold (we have single-bounce @0.90 retention; grazing-only trigger).
2. Mask-top fix: AC sputter 0.001@200 — replaces our over-strong mask carbon
   removal; with 1e-5 O-inertness this pins hm≈850.
3. PC dynamics: sputter 0.6@50 + de-crosslink 0.3@8 (both directions now exact).
4. Ion-activated SiO2 chemisorption (0.9 activation; 0.8-0.9 on activated) —
   likely raises F delivery and may re-balance mouth vs floor.
5. Reconcile appendix-converged vs paper-optimized values (0.2 vs 0.0842 AC
   deposition; 0.0423 vs 0.0628 O-etch; 0.0852 vs 0.0909 bare) — the fig-7
   base case likely uses the OPTIMIZED set; declare which set petch targets.

## Huang thesis (local, 372 pp) — MCFPM scattering semantics VERBATIM (2026-07-28)

Selection rule (§2.3, verbatim): per collision the process is randomly
selected from a probability array; "if the cumulative yield of non-reflective
processes is less than unity, then the elastic scattering yield is increased
so that the cumulative yield is unity. If the cumulative yield is greater
than unity, the elastic yield is reduced so that the cumulative yield is
unity" (then all scaled if still >1). => reflection weight = max(0, 1 - sum
of reaction probabilities) per event; particle consumed when reactions
saturate.

Eq. 2.34 energy retention with **E_ts = 100 eV, E_c = 10 eV, theta_c = 70
deg**: E_i > 100 eV AND theta > 70deg (from normal) => PURE SPECULAR, retain
ALL energy; theta < 70deg or E_i < 10 eV => diffusive scatter (random
direction, interpolated loss per Eq. 2.34); between => interpolation.
Multi-bounce until consumed. ("specularly reflected ions to as low as a few
eV" survive after charging deceleration; "majority of the energy is
maintained" on sidewalls.)

Open interpretation question: the appendix angular-form column ("1") for
polymer rows vs the paper's Kress attribution — whether the SELECTION
probability carries the angular yield factor decides wall transmission at
78-87 deg. DISCRIMINATOR (no feature run needed): Huang reports hot-neutral
vs ion flux fractions at the etch front vs AR (thesis chs. 4-5 text around
lines 5089-5445 of extraction) — implement the 2-3 candidate interpretations
as a 1-D wall-cascade model and pick the one reproducing HIS reported
funneling efficiency; that is a published-observable calibration of the
semantics, not a knob.

Implementation (exact, replaces the three tested variants):
1. reflection_weight = max(0, 1 - sum(reaction probabilities at E, theta))
2. specular full-E when E>100 & theta>70; Eq-2.34 interpolated loss else;
   diffusive random direction below cutoffs
3. multi-bounce (iterate splitter on secondaries until weight < epsilon)
4. y-symmetrize secondary launches (extrusion-guard fix, ml12 trip)
