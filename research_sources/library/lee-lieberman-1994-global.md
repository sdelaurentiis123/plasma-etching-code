# lee-lieberman-1994-global

**Lee & Lieberman, conserved 0-D global balances for Ar/O2/Cl2 plasmas**

- **Citation:** C. Lee and M. A. Lieberman, “Global Model of Ar, O2, Cl2 and
  Ar/O2 High Density Plasma Discharges,” UCB/ERL M94/49, revised 28 November
  1994; journal version, *Journal of Vacuum Science & Technology A*.
- **DOI:** `10.1116/1.579366`
- **Primary record:**
  `https://digicoll.lib.berkeley.edu/record/134386`
- **PDF SHA-256:**
  `f2049e7041984d658d23688e8e8112a8d8e8a524172a8d2e335be8fde7fc2e23`
- **Project extract:**
  `research_sources/thesis_extracts/lee_lieberman_1994_global_model.txt`
- **Visual audit:** equations/assumptions at 240 dpi; Table 3 at 300 dpi;
  Figures 3 and 8 at native 180 dpi review. The source PDF is not
  redistributed.
- **Status:** PRIMARY FULL TEXT READ + EQUATIONS/TABLE 3 VISUALLY AUDITED
- **Topic:** chemistry-agnostic global-model structure, argon reaction deck,
  cylindrical wall losses, and independent reactor-validation route

## Claims table

| # | verified source claim | use and boundary |
|---|---|---|
| Q1 | The model is a volume-averaged cylindrical TCP reactor with `L = 7.5 cm`, `R = 15.25 cm`, and neutral temperature 600 K. | Fixed geometry/thermal inputs for the published-model reproduction. These are not portable Krüger dimensions. |
| Q2 | Total absorbed power is divided into electron volume loss, ion wall loss, and electron wall loss. Electron wall energy is `2 Te`; ion wall energy is stated as typically `5–8 Te`. | Requires an explicit power ledger. The `5–8 Te` statement is a model range, not a uniquely measured constant. |
| Q3 | Positive-ion continuity combines electron-impact production, Bohm wall loss, ion-ion recombination, and asymmetric charge exchange; sheath-edge quasineutrality closes electron loss. | Defines the conserved particle/current structure for molecular extensions. |
| Q4 | `A_eff = h_L 2πR² + h_R 2πRL`, with Eqs. 13–14 giving axial/radial edge factors; `1/Λ² = (π/L)² + (2.405/R)²`. | Landed in `petch.reactor_global.geometry`, with the high-pressure diffusion terms explicit rather than silently dropped. |
| Q5 | Table 3 gives the five argon volume rates: `1.23e-7 exp(-18.68/Te)`, `3.71e-8 exp(-15.06/Te)`, `2.05e-7 exp(-4.95/Te)`, `2.0e-7`, and `6.2e-10 cm3/s`. | Frozen, source-unit-preserving argon deck inputs. They are empirical/regressed literature coefficients, not first-principles cross sections. |
| Q6 | Table 3 gives ion Bohm wall loss and metastable diffusive wall loss as first-order rates. | The forms are importable only after the mean-free-path/diffusivity inputs carry primary provenance. |
| Q7 | Table 3 prints pooling as `Ar* + Ar* -> Ar + Ar+`. | Literal import violates charge. The closed ledger must make the physically emitted electron explicit and record this completion. |
| Q8 | Figure 3 reports `Te(p)` at 1000 W, 35 SCCM, zero recombination coefficient; argon `Te` decreases monotonically with pressure. | Published-model reproduction board, not independent validation. |
| Q9 | Figure 8 compares model curves with three argon experiments using each experiment’s geometry and operating conditions. | The symbols cannot be pooled anonymously; Ra/Mahoney/Oomori primary conditions must be recovered and preregistered separately. |

## Use decision

This source is the first executable verification authority for the independent
global-reactor stack. It supplies universal balance structure and a compact
argon deck, not a C4F6 mechanism and not Krüger’s missing wafer boundary.

The implementation therefore separates:

1. exact atom/charge/unit bookkeeping;
2. source-specific empirical rate decks;
3. geometry, wall, pumping, and power closures; and
4. independent experimental validation.

The report’s own Figure 3 can verify faithful reimplementation. A predictive
claim begins only on condition-specific experimental data that did not select
the deck inputs.
