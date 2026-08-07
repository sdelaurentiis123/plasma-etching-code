# tinacba-2021-jvstb

**Tinacba et al., DFT-informed SF5+ MD against mass-selected beam depths**

- **Citation:** E. J. C. Tinacba, T. Ito, K. Karahashi, M. Isobe, and
  S. Hamaguchi, “Molecular dynamics simulation for reactive ion etching of Si
  and SiO2 by SF5+ ions,” *Journal of Vacuum Science & Technology B* **39**,
  043203 (2021).
- **DOI:** `10.1116/6.0001230`
- **Retrieval route:** author-hosted full text,
  `https://researchmap.jp/satoshi_hamaguchi/published_papers/37820817/attachment_file.pdf`
- **Status:** FULL TEXT READ; audited extract:
  `research_sources/thesis_extracts/tinacba_2021_sf5_audit.txt`; copyrighted
  PDF checksum-pinned but not redistributed.
- **PDF SHA-256:**
  `c0be3b475aa17b396c1f788baee14ba37b9026b264bb870dc2553055f27b31ad`
- **Topic:** atomistic surface provider, mass-selected SF5+ beam, absolute
  depth per ion dose, Si and SiO2

## Claims table

| # | source evidence | use and boundary |
|---|---|---|
| Q1 | “The mass and energy of incident ions were measured by an energy-mass analyzer” (p. 043203-2). | The surface boundary is species/energy identified rather than inferred from reactor knobs. |
| Q2 | The Faraday cup occupied the sample location for current measurements; the sample depth was measured by contact profilometry; no radical beam was used (p. 043203-2). | Ion dose and depth are independent observables. The source equation `Y=dN/D` makes yield error exactly depth-per-dose error at common `N,D`. |
| Q3 | The S-F carrier potential was conditioned on B3LYP/6-311G DFT and a published S-F bond energy/distance; Si/O/F potentials predate this beam comparison (pp. 043203-2--6). | The MD markers were not fitted to the beam depth/yield and can be graded against the experimental markers. |
| Q4 | The evolving MD surface receives typically 4000 consecutive impacts before steady yield is evaluated (pp. 043203-3, 043203-7). | This is a finite-fluence atomistic provider, not a single-impact sputter coefficient. |
| Q5 | Figure 8 compares MD and experiment for Si and SiO2 at 150 and 2000 eV. The committed vision audit gives 5.88% mean and 15.04% maximum absolute depth-per-dose error across four overlaps. | Independent retrospective validation; no post-hoc PASS threshold is assigned. |
| Q6 | Figure 10 prints 2000 eV MD depth slopes of 12.5 nm for Si and 13.6 nm for SiO2 per `1e16 cm^-2`. | Independent nanometer conversion cross-check for the Figure-8 yield digitization. |
| Q7 | S-S, S-Si, and S-O reactions are deliberately disabled while S mass/radius and S-F carrier bonding remain (pp. 043203-3, 043203-11--12). | Forbids calling the provider a general reactive Si/O/F/S potential or transferring it into low-energy sulfur chemistry. |

## Consumed data

- `data/experimental/tinacba_2021/figure8_sf5_md_experiment.csv` retains all
  16 SF5+ MD/experiment markers, source pixels, setpoints, and digitization
  bounds.
- `scripts/digitize_tinacba_2021_sf5.py` verifies the 300-dpi source render and
  draws a full-resolution PIL overlay.
- `results/curated/tinacba_2021_sf5_depth/audit.json` grades the four exact
  MD/beam overlaps as absolute depth per fixed dose.
- `load_tinacba_2021_sf5_tables` and
  `TabulatedNormalIonRemovalMechanism` checksum-load the MD nodes into the
  common material router with exact `Si_atom`/`SiO2_formula` removal ledgers.
  The executable adapter refuses every source-absent dimension and leaves
  product routing unresolved.
- Nothing in this source supplies the species-resolved flux mixture or IEAD
  needed by the Yoshie or Krüger reactors.
