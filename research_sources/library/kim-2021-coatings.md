# kim-2021-coatings

**Kim et al., measured C4F6/Ar neutral and positive-ion mass/energy spectra**

- **Citation:** Y.-H. Kim, J.-S. Kim, D.-C. Kim, Y.-W. Kim, J.-B. Park,
  D.-S. Han, and M.-Y. Song, “Ion and Radical Characteristics (Mass/Energy
  Distribution) of a Capacitively Coupled Plasma Source Using Plasma Process
  Gases (CxFy),” *Coatings* **11**, 993 (2021).
- **DOI:** `10.3390/coatings11080993`
- **Retrieval route:** open publisher PDF,
  `https://mdpi-res.com/d_attachment/coatings/coatings-11-00993/article_deploy/coatings-11-00993.pdf`
- **Status:** FULL TEXT:
  `research_sources/thesis_extracts/kim_2021_c4f6_ion_radical_characteristics.txt`;
  archived PDF:
  `research_sources/kim_2021_c4f6_ion_radical_characteristics.pdf`
- **PDF SHA-256:**
  `9641195e267c298ae0334633540880845f523d40e8b2a08656cb1ec8df16fc22`
- **Text SHA-256:**
  `298110ad195f66ebc628b86b0013e9888753c2121584866173a71341c20293af`
- **Topic:** measured C4F6 parent signal, fluorocarbon positive-ion mixture,
  and mass-dependent ion-energy distribution at a wafer-facing electrode

## Claims table

| # | verbatim source claim | use and boundary |
|---|---|---|
| Q1 | “The mass/energy analyzer (EQP1000, Hiden, Warrington, UK) was placed at the center of the powered electrode” and a “100-µm aperture plate was installed 1 mm below the bottom electrode.” | Establishes a wafer-facing, mass/energy-resolved diagnostic plane. It is not an in-feature measurement. |
| Q2 | For the reported C4F6 spectra, “the gas ratio was fixed at 1:2, the total pressure was maintained at 20 mTorr, and the RF power was set to 300 W.” | Defines the measured reactor support: C4F6/Ar only, 40/80 sccm in the figure. It does not match Krüger's C4F6/Ar/O2, 10 mTorr, multifrequency high-power reactor. |
| Q3 | “In the case of the C4F6/Ar plasma, the radical density was measured in the order of Ar, C3F3, C4F6, CF, C3F4, and CF3.” | Direct evidence that an undissociated C4F6 mass signal can survive to the electrode diagnostic in a C4F6 plasma. The plotted SEM count rate is not an absolute C4F6 wafer flux and is not imported as one. |
| Q4 | “In the case of C4F6/Ar plasma, the order was CF+, C3F3+, CF3+, CF2+, and C3F5+.” | Directly rejects treating an aggregate `Ions` boundary as a single Ar-like projectile. The ordering is reactor-specific and cannot identify Krüger's mixture. |
| Q5 | “smaller ions have an energy distribution that includes high energy” while “as the ion mass increases, the energy distribution mainly converges to approximately 270 eV.” | Ion identity and IEAD are coupled; one aggregate IEAD plus one effective projectile mass cannot reproduce the measured species structure. The authors explicitly defer the cause to future verification, so no causal law is imported. |

## Consumed evidence and prohibition

Figures 4 and 5 were rendered and visually inspected at 300 dpi. Figure 4
shows a distinct C4F6 peak in the neutral/radical mass spectrum; Figure 5
shows different energy spectra for multiple light and heavy fluorocarbon
ions. No curve is digitized into a petch default because the ordinate is
secondary-electron-multiplier count rate, the ionization/fragmentation
response is species dependent, and the source does not provide the
calibration needed to turn those counts into absolute incident fluxes.

This paper strengthens the depth-identifiability diagnosis: the two omitted
boundary classes are directly observed in a C4F6 plasma. It does not close
Krüger's boundary because the chamber, gas mixture, power, pressure, and
diagnostic calibration differ.
