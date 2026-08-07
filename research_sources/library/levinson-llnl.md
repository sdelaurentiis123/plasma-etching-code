# levinson-llnl

**Levinson, Shaqfeh, Balooch & Hamza controlled chlorine-beam corpus**

- **Primary article:** J. A. Levinson, E. S. G. Shaqfeh, M. Balooch, and
  A. V. Hamza, “Ion-assisted etching and profile development of silicon in
  molecular chlorine,” *Journal of Vacuum Science & Technology A* **15**,
  1902–1912 (1997).
- **DOI:** `10.1116/1.580658`
- **Primary full-text route:** author-uploaded article transcription and
  figures at
  `https://www.researchgate.net/publication/224449171_Ion-assisted_etching_and_profile_development_of_silicon_in_molecular_chlorine`.
- **Follow-up:** J. A. Levinson et al., “Ion-assisted etching and profile
  development of silicon in molecular and atomic chlorine,” *Journal of
  Vacuum Science & Technology B* **18**, 172–190 (2000),
  DOI `10.1116/1.591170`; official abstract at OSTI 20215307.
- **Status:** PRIMARY FULL TEXT ONLINE for the 1997 article; article pixels
  and thesis not archived. The article is sufficient for a controlled-beam
  surface/transport board but not an absolute feature-depth prediction
  because it omits case-specific ion flux and exposure time.
- **Topic:** Ar+/Cl2/Si beam yields, coverage-dependent sticking, Knudsen
  neutral transport, and controlled-beam feature profiles.

## Claims table

| # | verified source claim | use and boundary |
|---|---|---|
| Q1 | The MIBE apparatus sends a nearly monoenergetic Ar+ beam normal to the sample while an effusive doser supplies an isotropic room-temperature Cl2 background. A Faraday cup at the sample position measures ion intensity. | This removes reactor/sheath ambiguity and is a strong independent transport board. |
| Q2 | The source is designed for 10–1200 eV and about 0.45 mA/cm2 at 100 eV; at least 70% of emitted ions lie within ±10 eV of the anode setting and the distribution is near Gaussian. | These are apparatus capabilities, not the case-specific Figure 11 ion flux or fluence. They must not be substituted for the unreported run value. |
| Q3 | The paper models coverage as `Q = 1 / (1 + A Ychem GI / (2 S0 GN))` and total yield as `Ysput(1-Q) + Ychem Q`. It regresses `S0=0.75` and `A=2.0` to its planar yield data. | This is a mechanistic beam regression, not an ab initio or fit-free surface law. It is independent of the feature SEMs and can be transferred retrospectively with that label. |
| Q4 | Its feature solver uses line-of-sight molecular flow and coverage-dependent diffuse re-emission. The source explicitly treats the coupled coverage/neutral-flux problem as nonlinear and advances the interface with fourth-order Runge–Kutta. | Directly tests petch’s deterministic neutral-exchange and evolving-interface core. |
| Q5 | Figure 11 reports three 100 eV feature cases: 2.67 µm width and 0.38 µm depth at ion/neutral ratio 0.004 with 0.97 µm initial mask; 5 µm width and 1.2 µm depth at ratio 0.008 with 0.8 µm initial mask; 1.9 µm width and 1.18 µm center depth at ratio 0.008 with 0.75 µm initial mask. | Exact caption values can seed a geometry/shape board. They cannot determine absolute time or fluence. |
| Q6 | The 1997 article gives neither the exposure time nor the measured case-specific ion current for those Figure 11 runs. It only says the current was measured and gives a typical saturated etch rate of about 15–20 Å/min. | Figure 11 cannot be an absolute-depth gate from the article alone. Choosing time so simulated depth matches the SEM would be target fitting. |
| Q7 | The measured features show minimal lag at low ion/neutral ratio, plus microtrenching and sidewall slope that the source simulation misses; the paper names ion reflection, oxide-mask redeposition, or mask charging as candidates. | A transport-only model may gate normalized ARDE and main profile shape, but must fail honestly on microtrenching unless those mechanisms are enabled. |
| Q8 | The follow-up adds atomic/molecular chlorine composition, surface recombination, ion-induced desorption, neutral-limited and ion-limited features, and direct-versus-isotropic neutral delivery. | This is the preferred extension after the 1997 molecular-Cl2 board, but only after full text and pixels are archived. |

## Numerical-use gate

The article-level audit is replayed by
`scripts/audit_levinson_1997_feature_identifiability.py`. Before any Figure 11
profile score is claimed, archive the original PDF or thesis, hash it, render
the SEM page at native/high resolution, verify it with PIL, digitize the
experimental contour separately from the authors’ light simulation overlay,
and save an overlay. Before any absolute-depth score is claimed, obtain the
case-specific measured ion fluence or both measured ion current and exposure
time. Neither the source’s apparatus maximum nor the measured target depth may
fill that missing boundary.
