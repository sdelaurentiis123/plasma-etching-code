# ji-2024-tio2-hierarchical

**Ji et al., same-gas TiO2/Cr two-hierarchical profile evolution**

- **Citation:** X. Ji et al., “One-Step Dry-Etching Fabrication of Tunable
  Two-Hierarchical Nanostructures,” *Micromachines* **15**, 1160 (2024).
- **DOI:** `10.3390/mi15091160`
- **Primary open full text:**
  `https://in.iphy.ac.cn/upload/s22/m/02669_2024121818070514127.pdf`
- **Retrieved PDF SHA256:**
  `92a57c600e93e4113e574ef286aae4a248d787c3445f8d21e2ff151c73647e81`
- **Local extraction:**
  `research_sources/thesis_extracts/ji_2024_tio2_hierarchical.txt`
- **Figure-3 page render:** PDF page 6, 600 dpi, 4961 by 7016 pixels,
  SHA-256
  `976db637feed609a182301266654e675a42d1422f6e7c08cff17ea3744d6c148`
- **Status:** PRIMARY CC-BY FULL TEXT READ + FIGURE 3 QUANTITATIVELY
  PIL-AUDITED; PROPOSED SURFACE REACTION SCHEME NOT A QUANTITATIVE KINETIC LAW
- **Topic:** TiO2/Cr feature etching, CHF3/SF6/O2, mask shrink, lower-feature
  passivation, and physical wide-base profiles

## Claims table

| # | verified source claim | use and boundary |
|---|---|---|
| Q1 | The TiO2 run used an ICP tool at `350 W` source power, `120 W` RF power, `40/10/5 sccm SF6/CHF3/O2`, `10 mTorr`, and `40 C`, on `800 nm` electron-beam-deposited TiO2 with `60 nm` Cr. | This is a close constituent/material/mask profile witness, but it is an ICP, SF6-rich, lower-pressure process on a different film. No flux, rate, bias, or coefficient transfers to the Oxford CCP condition. |
| Q2 | A `200 nm` CD, `300 nm` pitch TiO2 grating evolves from a vertical profile to a trapezoidal upper region over a rectangular lower region and finally to an upper triangle over a lower rectangle. | Direct experimental evidence that a positive bottom-minus-top width and a two-zone profile can be physical in TiO2/Cr under CHF3/SF6/O2. It does not validate the magnitude or mechanism in the conditional Oxford profile. |
| Q3 | The authors report Cr linewidth loss and lower-feature passivation occurring together; lateral Cr shrink is visible once TiO2 depth reaches about `300 nm`. | Requires explicit Cr geometry evolution and a passivating/growth state before the reported shape may be mechanistically reproduced. A pinned-mask, removal-only model cannot claim this pathway. |
| Q4 | Increasing ICP power changes the profile from upper-triangle/lower-rectangle through near vertical at `500 W` to upper-rectangle/lower-trapezoid at higher source power; changing RF power also changes the upper height, tip radius, and gap. | Shows the sign of the profile response is plasma- and surface-state-dependent. It rules out treating a generic flare as a universal geometry-only law. |
| Q5 | At fixed `200 nm` CD, designed gaps at or above `100 nm` leave the reported morphology metrics comparatively stable, while the `70 nm` gap changes the upper/lower heights and angle. | Qualitative pattern-loading/transport validation target. It cannot identify radial reactor nonuniformity. |
| Q6 | The paper proposes TiF4 passivation, sputter removal, and bottom accumulation as its mechanism. It does not report species-resolved fluxes, energy/angular distributions, reaction probabilities, TiF4 coverage/thickness metrology, or sputter yields. | Mechanism topology only. The printed CF3-negative-ion sputter channel is not installed: a powered negative wafer does not receive an unconstrained bulk negative-ion beam, and no target flux is measured. |
| Q7 | The checksum-bound Figure-3 RF board gives upper-triangle height `143.7, 206.9, 269.8, 347.8, 336.1 nm`, tip radius `99.9, 84.9, 15.0, 9.8, 9.8 nm`, and gap `96.0, 80.0, 69.9, 42.1, 18.0 nm` at `90, 120, 150, 180, 210 W`. | `JI-2024-TIO2-FIG3-RF-MORPHOLOGY-R1`; target-free response gate for a TiO2 mechanism. Strict gap narrowing requires retained/deposited surface volume and rejects a removal-only explanation for this experiment. The unspecified etch time and different ICP/material condition prohibit Oxford coefficient transfer. |

## Use in petch

This paper is the strongest direct reason not to delete a converged wide-base
TiO2 shape merely because it looks unusual. It is equally strong evidence that
the current rate-normalized Oxford sentinel cannot attribute such a shape to
the real target: that sentinel disables mask erosion, passivation, growth, and
redeposition. The source therefore supports a future TiO2 mechanism topology
with a fluorinated/passivated surface inventory and evolving Cr mask, while
leaving every target rate coefficient explicitly unresolved.

Figure 3(b1-b3) was inspected at the original 600 dpi render. Fifteen red
marker centers and three ordinate calibrations replay through
`scripts/digitize_ji_2024_figure3.py`; the color overlay was reconciled at
original resolution. The committed digitization bounds are pixel-localization
bounds, not experimental error bars, because the source reports none.
