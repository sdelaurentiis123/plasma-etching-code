# lam-direct-drive-patents

**Lam Research, DirectDrive and RF-reference patent family**

- **Primary patent records:**
  - `US10515781B1`, “Direct drive RF circuit for substrate processing
    systems”;
  - `US12165841B2` / `US20220199365A1`, “Dual-frequency, direct-drive
    inductively coupled plasma source”; and
  - `US12106947B2`, “RF reference measuring circuit for a direct drive
    system supplying power to generate plasma in a substrate processing
    system”; and
  - `WO2024015694A1`, “Plasma detection in semiconductor fabrication
    apparatuses.”
- **Primary host:** `https://patents.google.com`
- **Status:** PRIMARY PATENT FULL TEXT READ
- **Topic:** matchless RF source control, output-node sensing, and power
  calibration boundaries

## Claims table

| # | verified source claim | use and boundary |
|---|---|---|
| Q1 | The direct-drive circuit places current and voltage sensors at its output and calculates their phase offset to tune drive frequency. | Supports a declared RF-output measurement node and rapid state tracking. Phase control alone is not absorbed-power measurement. |
| Q2 | The patent contrasts direct frequency control with motor-driven variable capacitors that are too slow for rapid level-to-level pulsing. | Time-resolved reactor inputs must retain pulse states; a time-averaged setpoint can lose physically important transitions. |
| Q3 | The dual-frequency family describes direct-drive ICP without a traditional impedance match between source drive and coil/chamber, including separately driven source structures. | A direct-drive boundary is distinct from forward-minus-reflected matched RF, but both still terminate upstream of plasma absorption. |
| Q4 | The RF-reference patent switches the direct-drive output to an impedance-transforming reference circuit, RF power meter, and dummy load. | Provides a public calibration architecture for real source-output power. It does not measure coil/window/electrode heating during plasma operation. |
| Q5 | Lam’s official Akara page says DirectDrive is used in Akara and responds over 100 times faster than earlier sources. | Connects the patent family to a marketed product lineage only at the level Lam states publicly; patent embodiments are not asserted to be the exact production BOM. |
| Q6 | `WO2024015694A1` describes V/I or phase-magnitude sensing at a matchless-source output, at a conventional match output, or at upstream nodes combined with a lumped match model; it also describes per-coil current sensing. | Publicly establishes several Lam RF-observation nodes and multi-coil telemetry. Effective coil-load impedance detects plasma state but does not by itself separate plasma resistance from coil and downstream hardware resistance. |

## Use decision

The patents motivated `petch.reactor_global.power`:

- measured absorbed power is predictive directly;
- matched RF needs reflected-power and hardware-loss accounting; and
- DirectDrive may provide a cleaner calibrated output node, but still needs
  downstream coil/window/electrode loss subtraction.

The implementation intentionally rejects the shortcut
`DirectDrive output power == plasma absorbed power`.
