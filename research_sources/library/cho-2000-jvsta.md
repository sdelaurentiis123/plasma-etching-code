# cho-2000-jvsta

**Cho et al. (2000)**

- **DOI/URL:** 10.1116/1.1318193
- **Retrieval route:** publisher
- **Status:** ABSTRACT-ONLY
- **Topic:** angular-laws — Angular yield laws (class-1 physical, class-2 chemical)

## Claims table

Rows relocated verbatim from the repo's research/results docs. `consumed by` = the
doc that recorded it (that doc states which constant/decision it fed).

| # | recorded claim (as written in our docs) | consumed by |
|---|---|---|
| Q1 | - **Cho, Hwang, Lee, Moon, *JVST A* 18, 2791 (2000)**, DOI 10.1116/1.1318193, CF₄, 5 mTorr, bias −100 to −800 V, Faraday cage: *[abstract]* "The normalized etch-yield curves showed virtually the same angular dependence regardless of the ion incident energy. The curve shape was similar to that of physical sputtering **except that the ratio of the maximum yield to that at 0° was as low as about 1.3.** … partly attributed to the fluorocarbon polymer film, which existed **as a few monolayers-thick film on the substrate surface at low angles near 0° but as a submonolayer at high angles between 45° and 75°.**" - **Lee, Hwang, Min, Moon, *JVST A* 20, 1808 (2002)**, DOI 10.1116/1.1503786, CHF₃, −20 to −600 V: *[abstract]* "**When the absolute value of the bias voltage was smaller than 200 V, the normalized etch rate … changed following a cosine curve** … When the magnitude of the bias voltage wa | `RESEARCH_LIP_CERTAINTY_2026-08-04.md`:230 |
| Q2 | petch's B=9.3 curve peaks at **52.6°** and is zero at 90° — the right shape. - **Provenance of the *value 9.3*: none found, and the source paper is off-domain.** Kress, Hanson, Voter, C. L. Liu, X.-Y. Liu, Coronell, *JVST A* **17**, 2819 (1999) is titled *[abstract, OpenAlex]* **"Molecular dynamics simulation of Cu and Ar ion sputtering of Cu(111) surfaces"**, for "ionized physical vapor deposition, used in **Cu interconnect technology**", over "**10–100 eV for Cu ions and 50–250 eV for Ar ions**". It is an MD study of a **metal** at **≤250 eV**. petch applies its angular parameter to a **fluorocarbon polymer** at **~1.5 keV**. `[VERIFY]` — I did not obtain the Kress body and found no source printing 9.3. - **Magnitude check against the in-chemistry measurements.** petch's `kress` is the yield per incident ion (petch's `atom_flux` is the *areal* flux onto the face and already carries cos | `RESEARCH_LIP_CERTAINTY_2026-08-04.md`:278 |
| U3 | [unquoted — verify on next use] \| oxide/mask angular class \| `B = 1.7` (peak 1.31) \| bounded by Cho 2000 / Schaepkens 1998 \| | `BENCHMARK_CERTIFICATION_2026-08-06.md`:17 |
| U4 | [unquoted — verify on next use] The only angular sputter measurements on SiO2 in fluorocarbon bound peak/normal at 1.30 (Cho 2000, JVST A 18, 2705) and 1.33 (Schaepkens 1998, JVST A 16, 3281). The oxide/mask rows now use B = 1.7 (peak 1.31, inside the measured band); the polymer row keeps Krueger's cited B = 9.3, so every validated lip/mouth result is untouched. `f(0) = 1` for any B, so all normal-incidence and blanket results are bitwise unchanged. Off-normal amplification drops 2.54x -> 1.16x, which should restore the timestep. | `RESULTS_LIMITING_REGIME_2026-08-05.md`:103 |
| U5 | [unquoted — verify on next use] Our `B = 9.3` gives peak/normal = 4.17 at 52.6 deg. The in-chemistry measurements bound the peak near 1.3 (Cho 2000) and 1.33 (Schaepkens 1998). Adopting a measured-bounded shape therefore *reduces* removal everywhere off normal and closes the top **faster** -- the opposite of what is needed. (Noted already in `RESEARCH_LIP_CERTAINTY_2026-08-04.md`; confirmed here.) | `RESULTS_LIP_REMOVAL_AUDIT_2026-08-04.md`:66 |

_Harvested 2 quoted + 3 unquoted mentions across the repo's docs._
