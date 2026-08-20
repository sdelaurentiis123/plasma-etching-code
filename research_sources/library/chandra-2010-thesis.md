# chandra-2010-thesis

**Chandra, "Capillary Force in High Aspect-Ratio Micropillar Arrays" —
the published PILLAR variant of the collapse criterion**

- **Citation:** D. Chandra, *Capillary Force in High Aspect-Ratio Micropillar
  Arrays*, PhD dissertation, Department of Materials Science and Engineering,
  University of Pennsylvania (Publicly Accessible Penn Dissertations, 2009/2010;
  advisor Shu Yang). Handle `20.500.14332/31394`.
- **Primary record:**
  `https://repository.upenn.edu/entities/publication/02a2a79a-790c-4fc0-a46f-0122533535cb`
- **Retrieval route that worked:** the DSpace-7 API content endpoint
  `https://repository.upenn.edu/server/api/core/bitstreams/16dc42c7-5336-4f5a-8d0f-369a8a0cd218/content`
  (the `/bitstreams/<uuid>/download` path returns an HTML shim, not the PDF).
- **Status:** PRIMARY FULL TEXT ARCHIVED —
  `research_sources/thesis_extracts/chandra_2010_upenn_thesis_capillary_pillars.txt`;
  **Eqs. 4.1, 4.2, 4.3, 5.7, 5.8, 5.9, 5.10a, 5.10b, 5.11 and Fig. 5.4b
  VISUALLY VERIFIED at 170 dpi from the PDF** (pdftotext mangles the radicals
  and the √2 in `w = √2p − d`).
- **Corresponding journal articles (paywalled, [unquoted]):** D. Chandra and
  S. Yang, "Capillary-Force-Induced Clustering of Micropillar Arrays: Is It
  Caused by Isolated Capillary Bridges or by the Lateral Capillary Meniscus
  Interaction Force?", *Langmuir* **25**, 10430–10434 (2009), DOI
  `10.1021/la901722g` (= thesis Ch. 5); D. Chandra, J. A. Taylor and S. Yang,
  "Replica molding of high-aspect-ratio (sub-)micron hydrogel pillar arrays and
  their stability in air and solvents", *Soft Matter* **4**, 979–984 (2008),
  DOI `10.1039/b800257a` (= thesis Ch. 4).
- **Topic:** capillary and adhesive collapse of 2-D pillar arrays; the only
  published closed-form pillar analogue of Tanaka's line criterion we located

## Claims table

| # | verified source claim | use and boundary |
|---|---|---|
| Q1 | Eq. 5.11 (p. 77), verbatim: `E'crit = 128γh³( 3h cosθ + w sinθ + √(9h² cos²θ + 3hw sin(2θ)) ) / (3π d³ w²)`, "the critical elastic modulus, E'crit, in case of Laplace pressure difference due to isolated capillary bridge between the four pillars is given by^{8,15}", with "Here, `w = √2 p − d` is the spacing between the diagonally opposite micropillars." | **The pillar criterion.** Circular cross-section, square lattice, 4-pillar diagonal cluster, isolated bridge spanning the full pillar height. Its bracket is algebraically identical to Tanaka's (Mack Eq. 4) with As = h/w. Cited to refs 8 (Chandra/Taylor/Yang Soft Matter 2008) and 15 (Tanaka 1993). |
| Q2 | Worked example (p. 77), verbatim: "in case of circular micropillar arrays of Chapter 4 (height 9µm, diameter 0.75µm, pitch 1.5µm), and assuming a contact angle θ = 60⁰… From Eq. 5.11, E'crit is calculated as 27 GPa". | Numerical self-check: our transcription reproduces **26.6 GPa** with γ = 72.7 mN/m — transcription verified to <2 %. |
| Q3 | Eq. 5.9/5.10a/5.10b (p. 74–75), verbatim: `Ecrit = 32√2 γ cos²θ h³ f(r)/(3d⁴)`; `f(r) = 1/(r−k) ( √(2/(k²−1)) + √(1/(2k²−1)) )`; `r = (1/k)[ (√2(k²−1)^{−1/2} + (2k²−1)^{−1/2}) / (√2(k²−1)^{−3/2} + 2(2k²−1)^{−3/2}) ] + k`; r = p/d. Fig. 5.4b: f(2) ≈ 3.2. | The **lateral-meniscus-interaction** branch (continuous liquid body, no isolated bridge). Reproduces the paper's 2 GPa for the same geometry. NOTE: as printed, Eq. 5.3 (N/m) and Eq. 5.7 (N) are dimensionally inconsistent with each other — our observation, flagged for verification against Langmuir 2009. |
| Q4 | Eq. 5.6 (p. 73) and text, verbatim: "for a typical case of θ = 60⁰ and aspect ratio of 10, the torque in case of isolated capillary bridge is at least 12 times greater than that from lateral capillary meniscus interaction." | Quantifies how much more severe the isolated-bridge branch is. |
| Q5 | Ch. 6 (p. 99), verbatim: "the micropillar arrays in both the geometries A and B are stable at elastic modulus of 1.2 GPa (75% MMA) and unstable for elastic modulus of 745 MPa (67% MMA)… by Eq. 5.11 (based on Laplace pressure due to isolated capillary bridge)… the critical elastic moduli are estimated as 16.7 GPa (geometry A) and 4.1 GPa (geometry B) which are much larger than 1.2 GPa for which the micropillars in both the geometries were found stable experimentally." | **The one quantitative experimental calibration of the pillar bridge criterion we have.** In a *continuous-liquid* drying regime it over-predicts the required modulus by ≈14× (A) and ≈3.4× (B) — i.e. it is conservative there. Does not license extrapolation to a true dry-out where isolated bridges do form. |
| Q6 | Ch. 5 (p. 71), verbatim: "It should be noted, however, that for 1D array of line patterns, Laplace pressure argument is applicable because in that geometry isolated liquid between the lines could exist, resulting in different Laplace pressures." | Explicit statement of the regime boundary between the two branches. |
| Q7 | Ch. 5 (p. 78), verbatim: "in later stages of liquid evaporation, there may no longer be a continuous liquid body surrounding the micropillars and isolated capillary bridges may likely form. However, such isolated capillary bridges will be near the base of the micropillars and thus will exert much less torque as compared to a capillary bridge spanning the whole micropillar length." | The physical reason Eq. 5.11 is an upper bound on the load. |
| Q8 | Eq. 4.2/4.3 (p. 56), verbatim, **adhesive** (dry) lateral collapse, attributed to ref. 8 (Glassmaker et al., J. R. Soc. Interface 1, 23 (2004)): `E*L = 8h⁴γs/(3w²d³)` (square pillars); `E*L = 2^{1/4}·32 h³ γs (1−ν²)^{1/4} / (3^{3/4} π d^{5/2} w^{3/2})` (circular pillars); `w` = spacing between two adjacent pillars, `γs` = surface energy of the pillar material. | **The published square-pillar formula** — but for *adhesion*, not capillarity. Useful as the post-dry stiction check. |
| Q9 | Eq. 4.1 (p. 56), verbatim, ground collapse: `E*g = 2^{11/2} 3^{3/4} (1−ν²)^{1/4} h^{3/2} W / (πd)^{5/2}`, W = work of adhesion. | Base/ground-collapse branch (Roca-Cusachs lineage). |
| Q10 | Ch. 4 (p. 53), verbatim geometry definitions: "h is the pillar height, d is pillar width or diameter, and w is the spacing between two adjacent pillars". | Symbol convention: in Ch. 4 `w` = adjacent spacing; in Ch. 5 Eq. 5.11 `w` = **diagonal** spacing √2p − d. Do not conflate. |

## Use in petch

Sole published pillar-geometry collapse criterion consumed by
`research_sources/RESEARCH_PATTERN_COLLAPSE_CRITERION_2026-08-20.md`; also the
cross-check that our Tanaka/Mack line transcription is faithful (identical
bracket). Q5 is the only experimental anchor available for the pillar branch.
