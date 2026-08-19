# choi-2013-tio2-cf4

**Choi et al., TiO2 bias and chemistry response in O2/CF4/Ar**

- **Citation:** K.-R. Choi, J.-C. Woo, Y.-H. Joo, Y.-S. Chun, and C.-I. Kim,
  “Dry etching properties of TiO2 thin films in O2/CF4/Ar plasma,” *Vacuum*
  **92**, 85–89 (2013).
- **DOI:** `10.1016/j.vacuum.2012.11.009`
- **Status:** PRIMARY FULL TEXT READ ONLINE; SOURCE-REPORTED FIGURES 2--6
  ENDPOINTS AUDITED; VERIFIED EXCERPT ARCHIVED; SOURCE PDF NOT LOCALLY ARCHIVED
- **Topic:** TiO2 fluorination, bias response, oxygen blocking, and
  ion-assisted product desorption

## Claims table

| # | verified source claim | use and boundary |
|---|---|---|
| Q1 | In `3/16/4 sccm O2/CF4/Ar`, `700 W` source power, and `2 Pa`, increasing DC-bias magnitude from `50` to `250 V` raises TiO2 rate from `130.9` to `197.2 nm/min` while TiO2:SiO2 selectivity falls from `0.65` to `0.56`. | Direct response-sign gate: a TiO2 law must respond to ion energy. The feed, ICP tool, and E-beam TiO2 film differ from Oxford and prohibit coefficient transfer. |
| Q2 | The source attributes the rate increase to enhanced bond breaking and sputter desorption of reaction products with increasing ion energy. | Requires an ion-assisted chemical removal channel; it does not identify its threshold, yield curve, or species partition. |
| Q3 | XPS reports Ti–F formation after fluorocarbon exposure and the oxygen sweep is nonmonotonic, with excess oxygen described as blocking fluorine adsorption. | Requires at least fluorinated and blocked/passivated surface states rather than an energy-independent removal scalar. |
| Q4 | The source reports no measurement uncertainty for the plotted response. | The two printed endpoints constrain sign and scale only; they cannot identify a unique square-root, thresholded, or coverage-dependent law. |
| Q5 | At fixed `16/4 sccm CF4/Ar`, `700 W`, `-150 V`, and `2 Pa`, the source reports TiO2 rates of `154.1`, `179.4`, and `137.5 nm/min` at `0`, `3`, and `9 sccm O2`. | The nonmonotonic response requires competing fluorine availability and oxygen blocking/cleanup in the model topology; it does not identify either probability. |
| Q6 | At fixed `3/16/4 sccm O2/CF4/Ar`, `-150 V`, and `2 Pa`, source power from `600` to `800 W` raises TiO2 rate from `136.0` to `208.3 nm/min` and TiO2:SiO2 selectivity from `0.66` to `0.83`. | This axis couples electron kinetics, ion density/energy, and radical production; it is a full reactor-to-surface response gate, not an isolated surface coefficient. |
| Q7 | At fixed feed, `700 W`, and `-150 V`, pressure from `1.2` to `2.8 Pa` lowers TiO2 rate from `187.7` to `138.7 nm/min` and selectivity from `0.60` to `0.50`. | A neutral-density-only model has the wrong qualitative closure; collisional sheath and ion-delivery response must remain coupled to neutral supply. |
| Q8 | AFM RMS roughness is reported as `36.5 Å` as deposited, `59.8 Å` after CF4/Ar, and `29.8 Å` after O2/CF4/Ar. | Requires a chemistry-dependent morphology score in addition to depth. Three endpoints do not uniquely identify roughening or smoothing physics. |

## Use in petch

`CHOI-2013-TIO2-BIAS-RESPONSE-R1` now records all source-reported endpoints
under the v2 multiaxis schema. A two-point square-root-of-bias decomposition is retained only as
an algebraic diagnostic and explicitly rejected as an identified surface law.
Together with the Ji same-gas morphology boards, this source fixes the minimum
future TiO2 topology: fluorination, neutral supply, ion-energy-dependent
desorption, competitive oxygen/passivation state, collisional delivery, and a
surface-morphology observable. It does not change the absolute Oxford prediction.
