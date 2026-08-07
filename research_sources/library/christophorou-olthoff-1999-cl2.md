# christophorou-olthoff-1999-cl2

**NIST critical evaluation of electron interactions with Cl2 and its
byproducts**

- **Citation:** L. G. Christophorou and J. K. Olthoff, “Electron Interactions
  With Cl2,” *Journal of Physical and Chemical Reference Data* **28**,
  131--169 (1999).
- **DOI:** `10.1063/1.556036`
- **Official NIST full text:**
  `https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=8765`
- **PDF SHA-256:**
  `a70bbed40bae014a551dd91fc96e322a258d959674ea53fe50ea2c6f4e020a6b`
- **Status:** PRIMARY NIST FULL TEXT READ; TABLE 25 PIXEL-AUDITED AT 300 DPI
- **Topic:** evaluated molecular-chlorine collision data, atomic-chlorine
  ionization, evidence gaps, and predictive reactor-deck boundaries

## Claims table

| # | verified source claim | use and boundary |
|---|---|---|
| Q1 | The review suggests cross sections for Cl2 total scattering, elastic scattering, total ionization, neutral dissociation, attachment, and ion-pair formation. | These evaluated tables are the preferred starting point for a predictive molecular deck; each channel still needs its own physical threshold and energy ledger. |
| Q2 | The authors identify missing measurements for Cl2 momentum transfer, vibrational excitation, electronic excitation, and dissociative ionization. | A “fundamental” deck cannot silently turn the legacy fitted electronic-excitation rows into measurements. These channels require declared uncertainty or a newer primary source. |
| Q3 | Three Cl2 electronic-excitation calculations disagree poorly, and the authors call for more experimental and computational work. | Lee--Lieberman Table 4 excitation rates cannot be elevated to truth merely because they close a power balance. |
| Q4 | Hayes et al. measured atomic-Cl single-ionization cross sections from threshold to 200 eV with ±14% absolute uncertainty; Table 25 lists selected values as the NIST suggested data. | Landed as a no-fit, tabulated Maxwellian rate provider. The quoted scale uncertainty propagates directly to the rate coefficient. |
| Q5 | Atomic-Cl excitation data discussed in the review are calculations for the 4s, 5s, 6s, 4p, 5p, 3d, 4d, and 5d states, not direct measurements. | Confirms that the first Lee Table 5 product denotes an excited `3d` state. Lee's printed extra electron violates charge conservation and is quarantined as a source defect. |
| Q6 | The review gives `D0(Cl2)=2.4793 eV` and a vibrational quantum of `0.0694 eV`; its threshold table gives 11.48 eV for Cl2 ionization, 15.45 eV for dissociative ionization, and 11.9±0.2 eV for ion-pair formation. | Physical energy-loss thresholds for later molecular-channel mapping; these are not interchangeable with Arrhenius fit exponents. |

## Executable decision

`nist_hayes_atomic_chlorine_ionization_rate()` integrates the audited Table 25
measurements analytically over a Maxwellian EEDF, enforces the NIST ASD
ground-state threshold, carries the ±14% scale uncertainty, and rejects
temperatures with a material unmeasured high-energy tail.

This evaluated provider does not overwrite the Lee--Lieberman reproduction
deck. The two remain side by side so source replay and predictive evidence
cannot be confused.
