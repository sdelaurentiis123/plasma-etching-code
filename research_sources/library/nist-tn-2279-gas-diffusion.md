# nist-tn-2279-gas-diffusion

**NIST TN 2279, evaluated self- and binary-diffusion coefficients for gases**

- **Citation:** D. R. Burgess Jr., *Self-Diffusion and Binary-Diffusion
  Coefficients in Gases*, NIST Technical Note 2279 (2024).
- **DOI:** `10.6028/NIST.TN.2279`
- **Primary record:** `https://doi.org/10.6028/NIST.TN.2279`
- **Status:** AUTHORITATIVE FULL PDF READ; TABLE 1a AND CORRELATION
  DEFINITIONS AUDITED
- **Topic:** neutral-gas diffusion, Ar-in-Ar self diffusion

## Claims table

| # | verified source claim | use and boundary |
|---|---|---|
| Q1 | The recommended temperature form is `ln(D / (cm2 s-1)) = A + B/(T/K) + C ln(T/K)` at 101.325 kPa. | Implemented directly, with exact unit and pressure conversions. |
| Q2 | For Ar in Ar, Table 1a gives `D298 = 0.182 cm2/s`, `A=-11.097`, `B=-45.486`, `C=1.676`, over 235–418 K. | Supplies Lee’s bulk metastable diffusion term without fitting a plasma observable. |
| Q3 | The report evaluates and recommends diffusion values at standard room temperature and standard atmospheric pressure. | Pressure scaling is explicitly the dilute-gas `1/p` law; finite-density corrections are outside this closure. |

## Use decision

Lee and Lieberman assume 600 K, outside the NIST Ar-in-Ar fit's stated
235–418 K range. The implementation records that fact as an extrapolation and
therefore labels the composed transport state `published_model`, not
`validated_model`. The Lee Knudsen term is harmonically combined with this
bulk value exactly as printed in the source model.
