# nist-argon-persistent-lines

**Evaluated Ar I resonance wavelengths and transition probabilities**

- **Source:** NIST Physical Measurement Laboratory, *Persistent Lines of
  Neutral Argon (Ar I)*.
- **Official table:**
  `https://physics.nist.gov/PhysRefData/Handbook/Tables/argontable3.htm`
- **Status:** PRIMARY EVALUATED DATABASE TABLE READ

## Claims table

| # | verified source claim | use and boundary |
|---|---|---|
| NA1 | The ground-state Ar I resonance line at 1048.21987 A has `Aki = 5.32e8 s^-1`, lower/upper `J = 0/1`. | Combined with the NIST A--f relation, gives absorption oscillator strength 0.262903 for deterministic Voigt transport. |
| NA2 | The ground-state Ar I resonance line at 1066.65980 A has `Aki = 1.32e8 s^-1`, lower/upper `J = 0/1`. | Gives absorption oscillator strength 0.0675468; preserves the line-specific natural width instead of one broadband opacity. |
| NA3 | The table attributes the transition probabilities to the M03 compilation and line positions/classifications to VHU99/M73. | Database provenance only; these constants do not determine reactor temperature, absorber density, emitter location, quenching, or feature depth. |

## Landed use

`ResonanceLineData` converts the evaluated wavelengths, A-values, and
statistical weights through the NIST atomic-spectroscopy relation.  Tests pin
both resulting oscillator strengths.  These are atomic inputs to the
deterministic homogeneous and axisymmetric resonance-escape solvers; Tian's
equipment-model trapping factors remain held-out implementation targets, not
experimental validation.
