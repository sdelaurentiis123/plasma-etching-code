# peko-2002-chf3-ion-molecule

**Measured reactive ion--molecule cross sections in CHF3**

- **Citation:** B. L. Peko, R. L. Champion, M. V. V. S. Rao, and J. K.
  Olthoff, "Measured cross sections and ion energies for a CHF3 discharge,"
  *Journal of Applied Physics* **92**, 1657--1662 (2002).
- **DOI:** `10.1063/1.1491276`
- **Official NIST record:**
  `https://www.nist.gov/publications/measured-cross-sections-and-ion-energies-chf3-discharge`
- **Official NIST PDF:**
  `https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=9195`
- **PDF SHA256:**
  `c0ab4fac39a611e364efc29c12632590d37e319cd8274382ae13be5a2d33d99c`
- **Local text extraction:**
  `research_sources/thesis_extracts/peko-2002-chf3-ion-molecule.txt`
- **Status:** PRIMARY FULL TEXT; FIGURE 2 SUM CURVE DIGITIZED AT 400 DPI,
  CHECKSUM-REPLAYED, AND ORIGINAL-RESOLUTION PIL OVERLAY AUDITED

## Claims table

| # | verified source claim | use and boundary |
|---|---|---|
| P1 | The cross-section abscissa is relative collision energy, not projectile laboratory energy. | Runtime converts a stationary-target CF3+ sheath energy by `E_rel = 70/(69+70) E_lab`; direct insertion of laboratory energy is forbidden. |
| P2 | Figure 2 reports the summed DCT cross section for CF3+ + CHF3; the reconstructed 21-point curve rises from about `2.63 A2` at `20.52 eV` to about `8.01 A2` at `195.85 eV`. | A C1 PCHIP in log-energy is executable only over the digitized support; the source's `+/-25%` DCT uncertainty and pixel-center allowance remain separate. |
| P3 | The summed CF3+ collision-induced-dissociation cross section is about `18 A2`; the three measured CID product groups are each about `6 A2` and approximately energy-independent for `20 < E_rel < 225 eV`. | Adds a measured reactive-destruction channel with the source's `+/-15%` CID uncertainty. It is not an elastic cross section. |
| P4 | At high measured energy, DCT plus CID makes CF3+ destruction substantial; the paper uses these reactions to explain discharge ion populations qualitatively. | Establishes that a collisionless CHF3 sheath is not a defensible final model at the Zhu target pressure. It does not supply a target-tool flux. |
| P5 | The article reports only selected reactive channels; it does not provide the elastic differential or total momentum-transfer cross section needed for an IEAD. | The executable provider is intentionally a reactive-removal floor. Complete molecular transport, target IEAD, and depth certification remain false. |

## Digitization receipt

- Script: `scripts/digitize_peko_2002_chf3_dct.py`
- Data: `data/experimental/peko_2002_chf3/figure2_cf3_chf3_dct_sum.csv`
- Manifest: `data/experimental/peko_2002_chf3/digitization_manifest.json`
- Committed CSV SHA256:
  `e10cc18dfa63e8fb79ccd386fba86bacc019b37f3714563b044b6b10c05105bb`
- Source pixels and the publisher PDF are not committed.

## Executable decision

`Peko2002CF3CHF3ReactiveCollisionModel` installs the measured DCT curve and
text-reported CID sum as a target-independent, differentiable destruction
kernel. It fails closed outside measured support and cannot be relabeled as a
complete scattering law or a TiO2 depth fit.
