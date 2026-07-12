# Kounis-Melas et al. Si-Cl-Ar DeepMD interaction data

Source dataset: A. Kounis-Melas, J. R. Vella, A. Z. Panagiotopoulos, and D. B. Graves,
“Data from *Deep Potential Molecular Dynamics Simulations of Low-Temperature Plasma-Surface
Interactions*,” DOI [10.34770/rjv6-2w31](https://doi.org/10.34770/rjv6-2w31), OSTI 2589032.
The associated paper is DOI [10.1116/6.0004027](https://doi.org/10.1116/6.0004027), OSTI 2514378.
Both the dataset and these copied result tables are CC BY 4.0.

This is **molecular-dynamics evidence, not experiment**. It describes Si-Cl2-Ar+ and must not be used
as a SiO2/fluorocarbon parameter source. Its purpose in petch is to provide a sourced second-chemistry
interaction contract through the same general transport/surface-table architecture.

## Reproducible acquisition

Princeton catalog JSON:
`https://datacommons.princeton.edu/discovery/catalog/doi-10-34770-rjv6-2w31.json`.

Archive:
`https://g-ef94ef.f0ad1.36fe.data.globus.org/10.34770/rjv6-2w31/480/DeepMDData.tar.gz`

- Published size: 595,544,258 bytes.
- Downloaded archive SHA-256:
  `4c9fa0b9268ac314da77b1012906dff4e45c5af79afd7ea674b26ace48e0f269`.
- Dataset README SHA-256:
  `780d413e66fab00f5e51e951b29ffc07f9a3bffeff986bf46e00d0a71dc3df7f`.

The three CSV files here are byte-for-byte copies extracted from that archive:

| Path in archive | SHA-256 | Meaning |
|---|---|---|
| `DeepMDData/Sputtering/Sputtering.csv` | `80ae627c1cec67258496ee7d22bd130817b678c1fd3288d5141436fcf374ee3c` | Normal-incidence Ar+ physical sputter yield and amorphous-layer thickness versus energy |
| `DeepMDData/RIE/RIE.csv` | `7cc634ae1218ba12d1e30ba7e6b4aefc0f4f0cc6de04ced8120115a60786cc77` | Si etch yield versus Cl2:Ar flux ratio at 100 eV normal-incidence Ar+ and 298 K |
| `DeepMDData/ALE/Products.csv` | `79a7cd3a2618a3fc3d65946d2db5247870d428b58270f78b0ffe46b5116bd9bf` | 80 eV ALE product yields versus Ar+ dose |

The archive also contains the 13,020-row `ALE/ALE.csv` trajectory, the final DeepMD model, LAMMPS
inputs, and the training corpus. They are not vendored because the pinned archive is already public and
the complete package is 596 MB. The small result tables needed by the interaction contract are retained
locally and checksum-gated.
