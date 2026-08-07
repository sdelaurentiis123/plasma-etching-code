# Hamilton et al. 2018 Cl2 dataset extraction manifest

## Source identity

- Primary source: J. R. Hamilton, J. Tennyson, J.-P. Booth, T. Gans, and
  A. R. Gibson, “Calculated electron impact dissociation cross sections for
  molecular chlorine (Cl2),” *Plasma Sources Science and Technology* **27**,
  095008 (2018).
- Article DOI: `10.1088/1361-6595/aada32`
- Data DOI: `10.15124/b11c65cf-2913-4c63-a522-2f57006cfb8a`
- Official York dataset filename: `Cl2_PSST_2018_Dataset.opj`
- Dataset license: CC BY.
- Downloaded OPJ SHA-256:
  `1e53fd091c1685f38326a9da5a1e78ba97387bc9c375783e8bb6c2a1d7fa0272`
- Published article PDF SHA-256:
  `7073912aeeeacedb82c75bcb3018e721f1b6ea7d596c721e2fd81d1808d57add`
- Grep-ready published-text SHA-256:
  `803cfbaa533888d5bcd58426a59f13bc6554ab26457bbc0c83e3a3ae15a0e175`

The OPJ and PDF are not redistributed. The exact numerical extracts and
bounded text extraction are retained under the source's CC BY terms.

## Exact OPJ extraction

- Reader: open-source `liborigin` / `opj2dat` 3.0.4.
- Parsed Origin version: `9.42`.
- Parsed objects: 45 datasets, 3 spreadsheets, 1 two-sheet Excel workbook.
- `Book3`, label `Fig2`: common 0.02--1000 eV energy grid plus the four
  dissociative Pi-state cross sections.
- `Book4`, label `Fig3`: common grid plus the four retained dissociative
  Delta/Sigma-state cross sections.
- `Book5`, label `Fig4`: authors' summed “New” dissociation cross section.
- `Book10/Sheet1`, label `Fig5`: authors' effective-temperature grid and
  Rescigno/Hamilton total rates for `x=1`, `x=2`, and `x=0.5` EEDFs.

Exact extracted artifacts:

| artifact | rows | SHA-256 |
|---|---:|---|
| `hamilton_2018_cl2_state_cross_sections.csv` | 50,000 | `7328d289542e23f2d12b4b172a19271120a3c5b62dc0dcd22a831569365dd288` |
| `hamilton_2018_cl2_total_dissociation_cross_section.csv` | 50,000 | `4645b1a97348103ff8a49a592197eb9f86f9f7c23728f1a8d703b0149785e93d` |
| `hamilton_2018_cl2_reference_rate_coefficients.csv` | 1,000 | `768ebc20052dc204d6dfbc733e407ab4229fd2d8da84175693d188442680a0de` |

The sum of the eight exported state cross sections agrees with the authors'
separate total column to a maximum relative difference of `0.2595%` at
positive-cross-section nodes. This is a source-internal rounding check, not
an experimental grade.

## Published-pixel audit

- Renderer: Poppler `pdftoppm`, 600 dpi.
- Figure-4 full page (PDF page 6) SHA-256:
  `5471b29261efd406ca56a5c8413d24f7bcb64a2a02cd3eeacb5731a9387a9c31`.
- Figure-5 full page (PDF page 7) SHA-256:
  `45392fddd4173098b7218aec30baccbbc9126a29c4d85769efe1c73c49c53224`.
- Visual inspection confirmed the “This work” total-cross-section and
  `x=1` Maxwellian-rate identities used in the OPJ mapping, their axes and SI
  units, and the 0--10 eV plotted rate domain.

## Executable reduction and gate

`scripts/build_hamilton_2018_cl2_state_rates.py` analytically integrates every
one of the 50,000 points for all eight states over the authors' own Figure-5
temperature grid. The sum reproduces their separately supplied `x=1` total
rate at 235 nodes from 0.3--5 eV:

- maximum absolute relative difference: `0.4443%`;
- mean absolute relative difference: `0.1855%`;
- numerical reproduction gate: `1%`;
- verdict: PASS.

The runtime table has 237 nodes including exact 0.3 and 5 eV endpoints.
Its SHA-256 is
`dce2b18a3a16dbdd8836d902d1628d2b427576d4d1646e35a2480cf7b035bb51`.
It uses bounded positive `log(k)` interpolation against `1/Te`, carries each
state's Table-2 vertical excitation energy into the power ledger, and refuses
extrapolation.

## Physics boundary

The eight retained state cross sections are fixed-nuclei R-matrix
calculations for ground-vibrational-state `Cl2(v=0)`. Above the ionization
potential, they are extended by transition-specific high-energy scaling.
Hamilton et al. state that Cosby's experiment sampled a distribution of
vibrational states, so it is not directly comparable. The paper publishes no
single uncertainty for the calculated cross sections. These rates are
therefore `semi_empirical`, Maxwellian-only within 0.3--5 eV, and are not
promoted to measured evidence.
