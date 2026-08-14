# turner-chabert-2014-rf-sheath

**Turner & Chabert, "A radio-frequency sheath model for complex waveforms"**

- **Citation:** M. M. Turner and P. Chabert, *Applied Physics Letters* 104,
  164102 (2014)
- **DOI:** 10.1063/1.4872172
- **Primary open manuscript:** https://arxiv.org/abs/1212.2612
- **Status:** PRIMARY FULL MANUSCRIPT READ; EQUATIONS 1--19 TRANSCRIBED AND
  EXECUTABLE
- **Topic:** reactor-sheath — arbitrary-current RF sheath closure

## Claims table

| # | source-grounded claim | consumed by |
|---|---|---|
| Q1 | Equations 1--8 close the time-averaged collisionless ion sheath from the Bohm current and Poisson equation, with the waveform-dependent mean/max voltage ratio `xi`. | `src/petch/reactor_global/current_driven_rf_sheath.py` |
| Q2 | Equations 10--14 give the instantaneous potential, electric field, displacement current, and moving electron-front position from the integrated sheath current. | `src/petch/reactor_global/current_driven_rf_sheath.py`; `moving_collisional_sheath_discrete_ordinates.py` |
| Q3 | Equation 15 defines `xi` as the cycle average of the normalized voltage waveform. | `TurnerChabertCurrentDrivenSheath.xi` |
| Q4 | For a single sinusoidal current, equation 18 gives `xi = 163/384`. | `tests/test_reactor_global_current_driven_rf_sheath.py` |
| Q5 | Equation 19 gives the sinusoidal-current maximum sheath-width scale.  The implementation uses the equivalent arbitrary-waveform equations 2--5 and 13, retaining the exact current-scale laws `s_max ~ J^3` and `V_max ~ J^4`. | implementation and exact JVP gate |

## Scope discipline

This source conditions on the current at the sheath boundary.  It does not
map generator forward power through a matching network and plasma impedance.
The implementation therefore accepts a de-embedded measured current or a
validated external-circuit result and rejects generator power as a substitute.
