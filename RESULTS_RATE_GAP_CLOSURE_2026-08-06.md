# The rate gap: closed as a question, not as a number (2026-08-06)

The 60 s endpoints (`RESULTS_SCORECARD_ENDPOINT_2026-08-06.md`) left one
defect wearing four costumes: the etch runs **2.38x too slow** relative to the
published 825 nm, and depth, late mouth over-closure, clog timing and the O2
ordering all follow from it. The mandate was to mine the un-searched space of
measured fluorocarbon-oxide yield data for a magnitude that closes it, land it,
and re-run the board.

**Outcome: the magnitude is not the defect.** An independent beam experiment
corroborates petch's channel magnitudes to 4.7%, and the channel a calibration
would act on is supply-bounded at the feature floor, so no value of it can
close the gap. No calibration was landed and no box time was spent re-running a
board that nothing changed.

## 1. The measurement that settles the magnitudes

Karahashi, *Hyomen Kagaku* **28**, 60 (2007) — the open-access full-text review
of Karahashi *et al.*, *JVST A* **22**, 1166 (2004), DOI 10.1116/1.1761119.
Mass-analyzed single-species ion beam, UHV, SiO2, 250–2000 eV, **no radical
flux**. Archived full text (not an abstract):
`research_sources/thesis_extracts/karahashi_2007_sio2_cfx_ionbeam.txt`,
PDF at `research_sources/karahashi_2007_hyomen_kagaku_28_60.pdf`.

Verbatim, L118-127 (Fig. 3, 1000 eV):

> F＋イオンに関しては 0.3 molecules／ion と低く，質量数の近い Ne＋イオン照射と
> 同程度である。したがって物理的スパッタによるエッチングと考えられる。
> …CF3＋イオンの場合エッチングイールドの値が 1.5 molecules／ion となる。

F+ gives **0.3 molecules/ion** and is physical-sputter-like (comparable to Ne+);
CF3+, the most reactive of the series, gives **1.5 molecules/ion**.

English abstract, L30-35: *"Etching yields of CFx+ increased with increasing ion
energy and with increasing number of fluorine atoms in the ions. **Above 1000
eV, etching yields is gradually saturated.** Below 500 eV, etching yields
abruptly dropped with decreasing ion energy, and fluorocarbon films grew on the
surfaces."*

### Cross-experiment check against petch (gated)

petch's two SiO2 channels carry Gray's beam-measured laws (MIT thesis 1993,
Ar+/F on SiO2, QCM) — a *different* apparatus, a different laboratory, eleven
years earlier, never cross-fitted to Karahashi:

| channel at 1000 eV | petch | Karahashi measured | ratio |
|---|---|---|---|
| physical sputter | 0.381 | 0.3 (F+) | 1.27x |
| chemically enhanced | **1.570** | **1.5 (CF3+)** | **1.047x** |

Gated in `tests/test_rate_gap_supply_bound.py`. The chemically enhanced
channel's **absolute per-ion magnitude** — the quantity a magnitude calibration
would move — is corroborated to **4.7%** by an independent measurement.

This does not contradict the earlier "2.8x too weak" reading recorded in
`BENCHMARK_CERTIFICATION_2026-08-06.md`: that comparison was against Gray's
*saturated plateau* at 350 eV, a coverage-curve quantity, not an absolute
per-ion yield. Both can hold. What the Karahashi check establishes is narrower
and is exactly what the mandate needed: the absolute magnitude is not the
defect, so raising it is not the fix.

Two further corroborations of landings this campaign already made, from the
same source:

- **√E scaling confirmed.** L200-215: the per-fluorine-atom yield is
  proportional to the 0.5 power of the energy allocated to that atom
  (*"エッチングイールドがイオンエネルギーの 0.5 乗に比例している"*), the
  Steinbrüchel extension of Sigmund. This independently refutes Krüger's
  anomalous `q = 1` on the two SiO2 rows, which commit `2f1e218` replaced.
- **Angular class confirmed.** L193-210: yields rise with incidence angle,
  peak near 60°, then fall through reflection; the 60°/0° ratio is ~2.0 for CF+
  and **1.3 for CF3+**. petch's oxide class-1 bound (B = 1.7, peak 1.31,
  commit `830e5c5`) sits on the measured CF3+ value.

## 2. Why a magnitude calibration cannot work

At HAR-floor delivery (neutral 0.10x, ion 0.70x, 3406 eV front) the complex
channel is carbon/fluorine-supply limited. Scaling its magnitude:

| complex-yield scale | floor rate | vs base |
|---|---|---|
| 1x | 4.100 nm/s | 1.00 |
| 2x | 4.390 nm/s | 1.07 |
| 4x | 4.518 nm/s | 1.10 |
| **8x** | **4.518 nm/s** | **1.10** |

Saturated by 4x. An 8x magnitude buys **10%** where the board needs **238%**.
At blanket delivery it saturates at 1.19x. Gated as
`test_complex_magnitude_cannot_close_the_floor_rate`.

This is the quantitative reason the exhaustive model-space solve
(`953f01c`/`f322cb4`) found an empty survivor set: the parameter the search was
ranging over has no authority at the operating point that matters.

## 3. What the target actually requires

Krüger's paper publishes the ion flux the thesis omits — *JVST A* **42**,
043008 (2024) Table I, `research_sources/thesis_extracts/krueger-2024.txt`
L298: **`Ions 1.2 × 10^16`** cm⁻²s⁻¹ = 1.2e20 m⁻²s⁻¹, which is the value the
bundled deck carries (`data/experimental/krueger_2024/base_case_boundary_fluxes.csv`).
Note the thesis's Table 6.1 (`krueger_thesis.txt` L4342-4352) has **no ion row
at all** — earlier 0-D probes in `scripts/` inferred 9.6e19 from Krüger's stated
band; the feature runs always used the published 1.2e20.

- 825 nm / 60 s over SiO2 at 2.2e28 m⁻³ = 3.025e20 units/m²/s.
- Against the published 1.2e20 ion flux: **2.52 units per incident ion**,
  sustained at the floor of an AR≈21 feature.
- Measured ceiling for the most reactive single ion at ≥1000 eV, saturating:
  **1.5 molecules/ion**.

The target is **1.7x the measured pure-ion ceiling**, required not at a blanket
surface (where radical synergy legitimately exceeds beam yields — petch reads
1.81 units/ion there) but at a floor whose radical delivery is attenuated ~10x.

### How the source reaches it

`krueger_thesis_2024.txt` L4570-4583, verbatim — Table 6.3 optimizer bounds and
Table 6.4 targets:

> ps,SiO2 Sputter probability of SiO2 **0.0 0.3**
> …
> Table 6.4:Target metrics for gradient descent optimization
> … hf Etch depth **825 nm**

The 825 nm depth is a **fitted optimizer target**, and the SiO2 sputter
probability was a free parameter ranged 0.0–0.3. Converged to 0.0852 under his
`q = 1` linear energy law it yields 4.06 units/ion at the front energy — a
supply-free channel large enough to hit the target. Under the measured √E law
the same row yields **0.752**. That 5.4x is the whole rate gap, and it lives in
an energy exponent that two independent beam experiments contradict.

## 4. Supply-side levers, measured and refuted

| lever | result | verdict |
|---|---|---|
| atomic F at the sourced band (2e20–1e21 m⁻²s⁻¹) | floor 4.100 → 4.112 nm/s | **+0.3%**, refutes the 5–25% estimate in `RESEARCH_F_SUPPLY_BAND` |
| activated-site fraction | θ_act self-limits at 0.16–0.20 at the floor | formation and consumption balance; not a bug |
| thermalized-ion return (E8) | +0.02% at 200:1, 8.9% floor composition | already recorded immaterial (`f83ede8`) |
| ballistic transport | AR-50 delivery 0.025287 vs Lam's published 2.5% | exonerated, 4-digit agreement |

The floor is **carbon/complex-site limited**: reactive FC arrival at 0.10x
delivery, chemisorbed at 0.278, bounds removals at ~0.71 per ion, and the
measured floor total is 1.07–1.18. Adding fluorine does not create complex
sites, which is why the F band moves nothing.

## 5. Candidate magnitude table

| candidate | value | provenance class | closes gap? |
|---|---|---|---|
| Gray β_e (landed) | 0.053(√E−√4) | **L1 beam-measured**, MIT 1993 | corroborated by Karahashi to 4.7%; supply-bounded at floor |
| Karahashi CF3+ | 1.5 molecules/ion @1000 eV, saturating | **L1 beam-measured**, Osaka 2004 | it is a *ceiling*, below the 2.52 required |
| Gogolides *et al.*, JAP **88**, 5570 (2000) via ViennaPS (Rodrigues 2023 Table 1, `tuwien_rodrigues_2023_fc_silica.txt` L264-288): A_n 0.0361, A_n/p 0.1444, A_s 0.0139, B 9.3, k_n 2 | A_n·k_n = 0.0722(√E−√4) | **model-fitted** (a surface-model paper, not a beam measurement); A_s and E_th match Gray exactly | 1.36x on the law only; still supply-bounded, and exceeds the measured ceiling 2.7x |
| Krüger `q=1` linear row | 4.06 units/ion @3406 eV | optimizer-fitted, bounds 0.0–0.3 | would close it, contradicts two beam experiments |

**Quarantine (full text not obtained, cannot support a landed number):** Kwon
*JVST A* **24**, 1906/1914/1920 (2006); Yin *JVST A* **26**, 161 (2008); Chae/
Vitale/Sawin (2003) 10.1116/1.1539085; Joubert *JVST A* **12**, 665 (1994).
Nothing in this document rests on them. The Karahashi measurement supersedes
their role for the magnitude question: it is the same class of measurement, on
the same system, over the energy range at issue, and we hold its full text.

## 6. Verdict

petch's SiO2 surface chemistry is validated against **two independent ion-beam
experiments** — Gray 1993 (its source) and Karahashi 2004/2007 (never fitted
to) — at 4.7% on the chemically enhanced channel and 27% on physical sputter,
with the energy exponent and the angular class independently confirmed.

The residual rate gap is not in that chemistry. Reaching 825 nm at the published
1.2e20 ion flux requires 2.52 units per ion at the floor of a 21:1 feature
against a measured ceiling of 1.5, and the source itself reaches it only by
fitting a sputter probability inside an energy law that measurement rejects.

The declared-calibration fallback was not landed because it is refuted rather
than merely unattractive: the constant it would move buys 10% where 238% is
needed. Landing it would have produced a number that agrees with Krüger's
figure and disagrees with every measurement of the underlying surface.

**The board stands as measured.** What would move it is boundary data that does
not exist in the source — a measured ion flux and a species-resolved IEAD for
this reactor, or an in-feature flux measurement at AR≈20 — which is the
partner-measurable path the dossier already names.
