# open-adas-cl0-vuv

**Fine-structure Cl I collision strengths, radiative branching, and observed levels**

- **Provider:** OPEN-ADAS, ADF04 type-3 files `cophps][cl/dw/ic][cl0.dat`
  and `nist][17/ic][cl0.dat` (downloaded 2026-08-08).
- **Atomic calculation:** AUTOSTRUCTURE distorted-wave Cl I calculation,
  producer H. P. Summers (2012); separate NIST observed-level conversion.
- **Official records:**
  `https://open.adas.ac.uk/detail/adf04/cophps][cl/dw/ic][cl0.dat`
  and `https://open.adas.ac.uk/detail/adf04/nist][17/ic][cl0.dat`
- **ADF04 manual:** `https://www.adas.ac.uk/man/appxa-04.pdf`
- **Pinned normalized physical-record SHA-256:** collision
  `f215093e9ab5ac36a202cdb353a0e3a3b9651982158ba105531fb44cab74c4e7`;
  NIST levels
  `12905140158c763f9a1cc6efff4d27fb74ecbf77e2162c67726b8dc528d1f430`.
- **Status:** PRIMARY DATABASE RECORDS READ; LICENSE-RESTRICTED RAW FILES NOT
  REDISTRIBUTED; CALCULATION-GRADE SENSITIVITY ONLY

## License boundary

The downloaded files state that OPEN-ADAS data are for the user's personal use
and may not be redistributed with modeling code, inserted into a managed
database, or used commercially without written permission. Consequently petch
contains only a caller-path parser, record hashes, derived audit outputs, and
this provenance card. Loading requires an explicit personal-use
acknowledgement. The raw files are not committed.

## Claims table

| # | verified source/data claim | use and boundary |
|---|---|---|
| OA1 | The collision file contains 538 fine-structure levels, 137267 type-3 transitions, a 100--2,000,000 K temperature grid, A values, and effective collision strengths. | Deterministic direct-coronal Cl I spectral sensitivity. |
| OA2 | The raw calculated first excited term is at 114967.7 cm^-1, above both the separately supplied observed first excited term at 71958.36 cm^-1 and the measured 104591 cm^-1 ionization limit. | Raw calculated wavelengths/thresholds are rejected; NIST observed separations are mandatory. |
| OA3 | Configuration/multiplicity/L/J plus repeated-state energy rank joins 45 resolved observed levels to the independent collision file. | Every unmatched or unresolved NIST blend is excluded; no guessed state assignment. |
| OA4 | At the Mahorowala reactor's 2.1--2.5 eV temperature proxies, direct atomic-Cl emission inside 104.82--106.67 nm is about 10^-15 cm^3/s, while 106.67--112 and 112--120 nm are about 10^-11 cm^3/s. | Du's measured 105-nm Si/photon yield cannot be spread over the calculated shortwave spectrum. |
| OA5 | The dominant matched shortwave sensitivity is near 118.88 nm, followed by a group near 109--111 nm. | These wavelengths define new surface-yield and wafer-spectrum experiments; they do not grant depth closure. |

## Use decision

This provider is not a production atomic deck. Distorted-wave collision
strengths, theoretical branching, incomplete observed-state mapping, cascades,
quenching, and resonance escape must be independently bounded. Its immediate
value is falsification: it proves that a scalar “sub-120-nm” photon channel is
not equivalent to the band in which the large Du yield was measured.
