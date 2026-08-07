# Wang–Olthoff 1999 Ar/Cl₂ ICP gate

**Status: frozen before a mixed Ar/Cl₂ solver or any model prediction exists**

Figure 9 provides 24 pixel-audited targets for a 13.56 MHz planar GEC ICP:
absolute total ion current and mass-resolved `Cl+` / `Cl2+` fluxes versus
pressure, for both pure chlorine and 20% chlorine in argon. Every point is
held out. No Wang–Olthoff observable may select a rate, wall coefficient,
power-transfer fraction, or species branch.

## Acceptance board

| observable | frozen acceptance |
|---|---|
| particle / charge residual | maximum normalized residual `≤ 1e-8` |
| electron-power residual | maximum normalized residual `≤ 1e-8` |
| total positive-ion flux | MAPE `≤ 30%`; every point `≤ 50%` |
| species fractions | mean absolute fraction error `≤ 0.10`; every species point `≤ 0.20` |
| pure-Cl₂ trends | total flux strictly increases; `Cl+` dominates at every pressure |
| mixed-feed trends | total flux strictly increases; `Cl2+` fraction remains below `0.05` |

The magnitude band is deliberately no tighter than inter-reactor
reproducibility discussed by the source. It is not an uncertainty claim:
Figure 9 publishes no chlorine error bars. The committed pixel allowance is
digitization uncertainty only.

## Power boundary and claim cap

The paper's 300 W is net power to the coil matching network. The authors state
that plasma-dissipated power is approximately 80% of the listed value. A 240 W
run is therefore a published-approximate reproduction boundary, not a direct
absorbed-power measurement.

Passing every gate would earn independent reactor/species validation
conditional on that approximate power boundary. It would not identify the
elementary `Cl2+`/`Cl+` electron-impact branch, reconstruct a Lam tool, transfer
to a dielectric CCP, or predict feature depth.
