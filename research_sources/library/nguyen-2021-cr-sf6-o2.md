# nguyen-2021-cr-sf6-o2

**Nguyen et al., Cr and CrOx etching using SF6 and O2 plasma**

- **Citation:** V. T. H. Nguyen, F. Jensen, J. Hübner, E. Shkondin,
  R. Cork, K. Ma, P. J. Leussink, W. De Malsche, and H. Jansen,
  “Cr and CrOx etching using SF6 and O2 plasma,” *J. Vac. Sci. Technol. B*
  **39**, 032201 (2021).
- **DOI:** `10.1116/6.0000922`
- **Status:** PRIMARY OPEN FULL TEXT READ; VERIFIED EXCERPT ARCHIVED
- **Topic:** chromium-mask oxidation, fluorination, inhibitor formation, and
  ion-assisted directional removal

## Claims table

| # | verified source claim | use and boundary |
|---|---|---|
| Q1 | Mixed SF6/O2 Cr etching peaks near `0.75%` SF6/O2; at `300 W` the source reports rates up to approximately `400 nm/min`, a `420 V` self-bias, and selectivity below `5` toward Si. | Cr loss is strongly coupled to feed ratio and bias; no constant physical-sputter yield can transfer to Oxford. |
| Q2 | The switched process uses an oxygen-rich oxidation step followed by an F-rich step. The proposed sequence oxidizes Cr, forms volatile CrO2F2 under fluorine, leaves a CrFx-inhibited surface, and requires ion impact to rupture/remove that inhibitor. | Fixes the minimum Cr state topology: pristine/oxide/fluoride-inhibited states plus neutral conversion and energetic cleanup. |
| Q3 | At `30 W`, the switched process reports selectivity above `20`, good directionality, and about `1 nm/cycle` (`7 nm/min`). | Directionality and selectivity depend on temporal chemistry and ion energy, not only total dose. Different tool, temporal feed, and substrate prevent coefficient transfer. |
| Q4 | The paper states that directional continuation requires sufficient bias to remove CrFx while avoiding physical sputtering of underlying bulk Cr. | Rejects both a neutral-only law and a pure-physical-sputter-only law as complete mechanisms. |

## Use in petch

The moving-Cr Oxford sensitivity board uses Janissen's transferred net
TiO2:Cr selectivity only as an explicitly conditional rate-normalized law.
This source sets its model-form warning: a future microscopic Cr deck must
carry CrOx formation, F-driven volatile conversion, a CrFx inhibitor, and
ion-assisted inhibitor removal. None of Nguyen's numerical rates is an Oxford
coefficient.
