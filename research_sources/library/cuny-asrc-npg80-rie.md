# cuny-asrc-npg80-rie

**ASRC CUNY Oxford PlasmaPro NPG80 RIE equipment record**

- **Owning-facility record:**
  `https://asrc.gc.cuny.edu/instruments/reactive-ion-etch-rie/`
- **Status:** PRIMARY FACILITY EQUIPMENT RECORD READ
- **Topic:** exact-tool capability and electrical/thermal configuration

## Claims table

| # | verified source claim | use and boundary |
|---|---|---|
| Q1 | The ASRC tool is an Oxford Instruments PlasmaPro NPG80 RIE for fluorine-based reactive-ion etching, accepting pieces and wafers up to 8 inches. | Fixes the tool identity and chemistry service class for the Zhu recipe. |
| Q2 | The listed source is 13.56 MHz RF, rated to 300 W, with a solid-state matching network. | Establishes the generator frequency and ceiling. The supplied 150 W recipe value remains forward-power demand, not absorbed plasma power. |
| Q3 | The active electrode has chiller/heater control from 0 to 80 C. | Supports the supplied 20 C table setpoint as a controlled boundary; it does not establish the actual wafer surface temperature. |
| Q4 | The facility describes this as a medium-density RIE system rather than an independently powered ICP-RIE. | The supplied recipe has one table-RF demand and no independent ICP power. An ICP source term must not be invented for this condition. |

## Use in petch

These are equipment facts, not plasma-state measurements. Electrode dimensions,
powered area, chamber volume, achieved self-bias, absorbed power, impedance,
species-resolved fluxes, and IEADs remain unknown. The exact-tool record therefore
constrains the reactor topology without closing its boundary conditions.
