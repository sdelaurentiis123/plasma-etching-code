# Malyshev 1998 Lam neutral-transport diagnostic

**Verdict: the missing Cl/Cl2 diffusion law is now reconstructed, but the
powered gas temperature and distributed wall state still prevent prediction.**

The source-parameterized Chapman--Enskog reconstruction, with Neufeld's
evaluated collision integral replacing the older Hirschfelder table cited by
Malyshev, gives
`D(298.15 K, 1 atm) = 0.151675 cm2/s`. This is
`1.12%` above the source's rounded
`0.15 cm2/s` room-temperature anchor after retaining its declared `1.25`
factor exactly. At the source-reported initial `333 K`, it gives
`D(1 atm) = 0.186215 cm2/s` and
`N D = 4.103975e+20 m-1 s-1`. The latter is
`33.91%` below
the old `6.21e20 m-1 s-1` Economou constant, which was not temperature-safe.

Each of the 23 measured-state Eq.-7 rows was then mapped through the exact
cylindrical Robin eigenmode at a declared `333 K` sensitivity. The powered
particle-pressure ledger retains the reported 95/5 Cl2/rare-gas inventory;
the rare gas is not incorrectly dissociated by Eq. 11. 22 rows admit a
physical effective wall probability, spanning
`0.01497--0.29455` with
median `0.05010`. One 11 cm / 10 mTorr /
200.704 W row is unattainable: its near-zero inferred dissociation requires a
wall-return frequency above the perfectly absorbing-wall limit. It is reported
as a physics failure, not assigned `gamma > 1`.

The source's fitted `gamma = 0.035`, replayed without refitting through the new
transport and measured electron-state board, has MAE
`9.878` percentage points and RMSE
`12.266` percentage
points in relative Cl2 density. This is not independent validation: Malyshev
fit that gamma to the same Figures 7--8 data using its older electron rate and
approximate transport closure.

The principal unresolved boundary is now explicit. The paper says the gas
starts at 333 K and heats with power, but publishes no powered gas-temperature
board. The incident-speed closure is also a thermalized sensitivity; Guha's
direct measurements warn that fresh Cl can remain nonthermal at low pressure.
Stafford's direct conditioned-wall data are at 300 K and cannot be silently
extrapolated to this distributed Lam wall. Consequently the effective
probabilities are sensitivity diagnostics only; no row supports a predictive
wall law, wafer flux, etch rate, or feature depth.
