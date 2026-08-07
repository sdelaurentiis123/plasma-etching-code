# benck-2003-c4f6

**Benck, Goyette & Wang, absolute mass-resolved C4F6/Ar ion flux and IEDs**

- **Citation:** E. C. Benck, A. Goyette, and Y. Wang, “Ion energy
  distribution and optical measurements in high-density, inductively coupled
  C4F6 discharges,” *Journal of Applied Physics* **94**, 1382–1389 (2003).
- **DOI:** `10.1063/1.1586978`
- **Primary metadata:** NIST publication record,
  `https://www.nist.gov/publications/ion-energy-distribution-and-optical-measurements-high-density-inductively-coupled-0`
- **Full-text retrieval:** author/NIST-origin AIP copy located through search;
  inspected locally but not redistributed because the PDF carries an AIP
  redistribution notice.
- **PDF SHA-256:**
  `96ef064afb0f804d3d853adcc275a63cd52507b0ea7e39f921d6cf6209e86d08`
- **Figure-9 page render:** PDF page 6, 600 dpi, 5209 by 6760 pixels,
  SHA-256
  `f423f6992e049ba40f47a540eb56450e9c080dddb4b03ad45ee8883b109a71fd`
- **Status:** PRIMARY FULL TEXT READ + FIGURE 9 QUANTITATIVELY PIL-AUDITED
- **Topic:** quantitative C4F6/Ar reactor-boundary validation for total ion
  current, ion composition, and mass-resolved IEDs

## Claims table

| # | verified source claim | use and boundary |
|---|---|---|
| Q1 | The experiment used a modified inductively coupled GEC reference cell with a grounded, water-cooled lower electrode and no lower-electrode RF bias. | Defines an ICP/grounded-surface validation case. It is not Krüger's dual-frequency, high-bias CCP. |
| Q2 | Ions were sampled through a 10 µm grounded side orifice 9.5 mm above the steel plate and 46 mm from the radial center. A Faraday cup at the same radius and height normalized the mass-spectrometer counts to absolute total ion current density. | This is a calibrated reactor-edge ion-flux measurement, not an inferred density or arbitrary mass-spectrum count. It is not a wafer-center or in-feature measurement. |
| Q3 | The energy resolution was 1 eV with an estimated ±1 eV scale uncertainty. Corrected ion transmission was estimated uniform to 20% over the measured mass/energy range. | Mandatory measurement uncertainty for a future reactor-model gate. |
| Q4 | Figure 9 reports absolute total and component ion fluxes for 25%, 50%, 75%, and 100% C4F6 in C4F6/Ar at 1.33 Pa (10 mTorr), 200 W. | A quantitative same-pressure C4F6 ion-composition board. O2 is absent and the power/source geometry differ, so the values cannot be transplanted into Krüger. |
| Q5 | In every C4F6/Ar mixture Ar+ was dominant. CF+ was the largest fluorocarbon ion, followed by CF2+; SiFx+/COFx+ etch-byproduct ions were also substantial. | Directly invalidates a single-projectile aggregate-ion closure and gives a target for a C4F6 reactor model. The byproduct-ion signal also proves two-way plasma/surface coupling in this apparatus. |
| Q6 | CF+ and CF2+ rose nearly proportionally from 25% to 75% C4F6 and then leveled, whereas CF3+ continued to rise; the total ion current decreased slowly with C4F6 fraction. | Ion composition is nonlinear in feed fraction. A feed-ratio interpolation is not a species-resolved boundary model. |
| Q7 | The measured IEDs were double peaked and mass dependent. C4F6 concentration increased the inferred plasma-potential oscillation amplitude. | Species identity, energy distribution, electronegativity, and sheath dynamics are coupled; one species-independent IEAD is not fundamental. |
| Q8 | The submillimeter diagnostic reports CF2/CF ratios, not absolute stable-C4F6 parent flux. | The paper advances the ion-side closure but does not supply the stable-parent wafer flux or C4F6/ion surface co-incidence law needed for Krüger depth. |
| Q9 | The checksum-bound Figure-9 digitization gives total positive-ion current `0.2644, 0.2164, 0.2012, 0.1557 mA/cm2` at `25, 50, 75, 100%` C4F6. CF+ rises `0.01319 -> 0.03153`, CF2+ rises `0.006133 -> 0.01752`, and CF3+ rises `0.002390 -> 0.008994 mA/cm2`; Ar+ falls `0.2067 -> 0.06973 mA/cm2` through 75% and is not plotted at 100%. | `BENCK-2003-C4F6-FIG9-ION-CURRENT-R1`; quantitative held-out target for a C4F6 reactor provider, retaining 10.1% digitization and 20% transmission uncertainties separately. |

## Vision audit and use

Figure 9 was inspected at the original 600 dpi render. The log ordinate is
positive-ion current density in `mA/cm2`; four mixture markers are shown for
the total and component ions. A checksum-bound 19-point table now retains the
unambiguous total, Ar+, CF+, CF2+, and CF3+ glyph centers. The three major
log-axis ticks and all full-page marker pixels replay through
`scripts/digitize_benck_2003_figure9.py`. The original-resolution color
overlay was visually reconciled. Overlapping C/F/Cx and SiFx/COFx glyphs are
excluded rather than guessed, and source pixels are not redistributed.

This source is now the preferred quantitative validation target for a C4F6
reactor model's ion side. It cannot normalize Krüger directly and it does not
close the stable-parent channel. Its main model-design consequence is that
porting a C4F8 global mechanism and calling C4F6 a “small delta” is not
authorized: C4F6-specific fragmentation, electronegativity, byproduct return,
and mass-resolved sheath transport must be graded against these measurements.
