# neufeld-1972-collision-integrals

**Lennard--Jones transport collision-integral correlations**

- **Citation:** P. D. Neufeld, A. R. Janzen, and R. A. Aziz, “Empirical
  Equations to Calculate 16 of the Transport Collision Integrals
  Omega(l,s)* for the Lennard-Jones (12-6) Potential,” *Journal of Chemical
  Physics* **57**, 1100--1102 (1972).
- **DOI:** `10.1063/1.1678363`
- **Official publisher record:**
  `https://pubs.aip.org/aip/jcp/article/57/3/1100/750261/`
- **Original-article scan inspected:**
  `https://www.scribd.com/document/832627226/Neufeld-Colision-Integral-Formula-ORIGINAL-PAPER`
- **Local implementation extract:**
  `research_sources/thesis_extracts/neufeld_1972_collision_integral.txt`
- **Status:** PRIMARY ARTICLE SCAN TEXT READ; EQ. 2 AND TABLE-I/II (1,1) ROW
  VERIFIED; LOCAL PUBLISHER-PDF PIXEL AUDIT STILL OPEN

## Claims table

| # | verified source claim | use and boundary |
|---|---|---|
| N1 | The authors calculate 16 reduced transport collision integrals for the Lennard--Jones (12-6) potential over `0.3 <= T* <= 100`. | The runtime refuses reduced-temperature extrapolation. |
| N2 | Eq. 2 is a four-term inverse-power/exponential correlation plus an optional sinusoidal term. | The `(1,1)` Table-I row has no sinusoidal coefficients, so binary diffusion uses only the four printed terms. |
| N3 | The `(1,1)` row is `A=1.06036`, `B=0.15610`, `C=0.19300`, `D=0.47635`, `E=1.03587`, `F=1.52996`, `G=1.76474`, `H=3.89411`. | These eight numbers are immutable source constants; no reactor or depth target selected them. |
| N4 | Table II gives 0.036% average deviation and 0.11% maximum deviation at `T*=0.75` for `(1,1)`. | The implementation records 0.11% as the correlation-fit ceiling; this does not cover Lennard--Jones parameter or experimental uncertainty. |
| N5 | The authors state their calculations are more accurate than the frequently used Hirschfelder--Curtiss--Bird values. | Neufeld is an explicit evaluated upgrade to the older table cited by Malyshev, not a claim that the source itself used Neufeld. |

## Executable decision

`neutral_transport.neufeld_1972_lennard_jones_omega_11` implements the exact
Table-I row and fails outside the source domain. The generic
`ChapmanEnskogBinaryDiffusivity` composes it with explicit masses,
Lennard--Jones parameters, mixing rules, evidence type, correction factor, and
uncertainty. A source correction remains provenance, not a hidden fit.

Malyshev 1998 cites the older Hirschfelder table for `Omega(1,1)*`. The petch
Cl/Cl2 provider uses Malyshev's printed equation, masses, Lennard--Jones
parameters, and immutable 1.25 measurement correction, but deliberately uses
Neufeld for the collision integral. It is therefore a source-parameterized
physics upgrade, not a source-identical table replay.

The original coefficient table has not yet been pixel-audited from a local
publisher PDF. The constant is therefore usable for transparent published-
model reconstruction, while predictive status additionally remains blocked by
the missing physical uncertainty of the Cl/Cl2 parameters and measurement
anchor.
