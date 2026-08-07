# Krüger aggregate-ion identity envelope

Status: **the missing ion identity spans the required absolute rate; the
reactor composition remains unidentified**.

Krüger publishes one aggregate positive-ion flux and one combined IEAD, not
species-resolved ion fluxes.  This audit sends three declared endpoints
through the same production-order 10 nm feature operator for 0.25 s:

1. nominal non-incorporating aggregate ions;
2. every published aggregate ion counted as CF₂⁺;
3. every published aggregate ion counted as CF₃⁺.

These are composition-envelope endpoints.  They are not proposed reactor
mixtures, and no endpoint or mixture fraction was fitted to feature depth.
The aggregate ion flux remains exactly the published `1.2e16 cm^-2 s^-1`.

Implementation revision: `8e1fb65`.

## Result

| declared aggregate-ion endpoint | depth at 0.25 s (nm) | mean prefix rate (nm/s) | difference from 825/60 rate |
|---|---:|---:|---:|
| non-incorporating / unresolved | 2.34807141 | 9.39228565 | -31.692% |
| all CF₂⁺ | 3.39053968 | 13.5621587 | -1.366% |
| all CF₃⁺ | 3.91113123 | 15.6445249 | +13.778% |

Krüger's reported 825 nm after 60 s corresponds to 13.75 nm/s, or
3.4375 nm over 0.25 s if one uses the run-average rate only as a scale
reference.  The CF₂⁺ and CF₃⁺ endpoints bracket that scale.  The result is not
a prediction of the 60 s endpoint because the experimental early-time rate is
not published and the feature rate evolves.

The mask opening is 86.832 nm for all three endpoints at this prefix, as
expected: the sensitivity changes the Guo oxide ion-incorporation ledger while
leaving the separate mask mechanism and transport boundary unchanged.

## Why ion identity changes removal

Guo's translating layer atom-counts incident ions.  An inert endpoint supplies
energy but no C/F atoms.  CF₂⁺ and CF₃⁺ supply energetic fluorine and carbon
coincident with the impact, changing the steady Si/O/C/F/vacancy state and the
ion-enhanced volatile-product channels.  This is not equivalent to multiplying
a yield: each face resolves a new atom-balanced state from its local
energy–angle measure and neutral/ion ratios.

The planar endpoint board independently gives 2.6131, 3.2986 and
3.6612 SiO₂/wafer-ion for non-incorporating, CF₂⁺ and CF₃⁺ respectively.
The feature calculation above shows that the same ordering survives the
850 nm mask's angular filtering.

## Numerical and provenance controls

All endpoints use:

- 10 nm grid and conservative common-refinement remap;
- 16 source positions, three face quadrature points;
- 8×16 neutral directions and 16 ion azimuths;
- eight radiosity rays/face and 1e-12 balance tolerance;
- identical digitized IEAD, neutral fluxes, mask mechanism and random seed;
- no ion-flux normalization, oxide-yield scale, or depth parameter.

Maximum radiosity relative balance residual: `6.96e-14`.
Material ledger residual: `0 units/m²`.

Artifact receipts:

- CF₂ config:
  `934a1180dd2678f2ab8ece1538a490e51a2179a44fb96e368f77726e96083d9e`
- CF₂ audit SHA-256:
  `7704c947430361685a3c66a25af1395fb1212beae7ed5be3f2771460c81ddfa5`
- CF₃ config:
  `00b49ea09f830be4a1cce61f84b143c615443706495af1233a844cd996c5e58c`
- CF₃ audit SHA-256:
  `135b9beaf4e4709df05c7e1c9972ca316bd037eef4b20c299e9271204f9cc4ef`

Full local artifacts:
`/private/tmp/krueger_guo_ion_cf2` and
`/private/tmp/krueger_guo_ion_cf3`.

## Scientific verdict

The answer to “can the model match the absolute depth scale?” is **yes within
the physically admissible but unpublished ion-composition envelope**.  The
answer to “has the 825 nm depth been predicted from Krüger's published
boundary?” remains **no**.

The discriminating measurement is species-resolved positive-ion flux plus
species-resolved IEAD at the wafer (or an independently validated reactor
model producing them).  Fitting an ion-mixture fraction to 825 nm would turn
this identified evidence gap into a target fit and is explicitly prohibited.
