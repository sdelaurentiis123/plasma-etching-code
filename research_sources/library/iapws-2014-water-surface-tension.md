# iapws-2014-water-surface-tension

**IAPWS R1-76(2014) — surface tension of ordinary water substance**

- **Citation:** International Association for the Properties of Water and Steam,
  *Revised Release on Surface Tension of Ordinary Water Substance*,
  IAPWS R1-76(2014), Moscow, Russia, June 2014 (6 pp.).
- **Primary record:** `https://iapws.org/documents/release/Surf-H2O`;
  PDF `https://iapws.org/public/documents/CH-L9/Surf-H2O-2014.pdf`
- **Underlying data:** the release states it is "a minor revision of the IAPWS
  release of 1994… The formulation itself is unchanged from the 1994 document",
  which in turn rests on the IAPS 1976 recommendation (Vargaftik lineage).
- **Status:** PRIMARY FULL TEXT ARCHIVED —
  `research_sources/thesis_extracts/iapws_r1_76_2014_water_surface_tension.txt`
- **Topic:** water–vapour surface tension vs temperature

## Claims table

| # | verified source claim | use and boundary |
|---|---|---|
| Q1 | Table 1, verbatim rows (t °C / experimental σ / uncertainty / calculated σ): 20 → 72.74 / 0.36 / 72.74 mN/m; 25 → 71.98 / 0.36 / 71.97 mN/m; 15 → 73.49; 30 → 71.19. | The γ used in every collapse-criterion evaluation. Pure water–vapour interface. |
| Q2 | Interpolating equation, verbatim (page 3, **visually verified at 150 dpi**): `σ = B τ^μ (1 + b τ)`, with `τ = 1 − T/Tc`, `Tc = 647.096 K`, `B = 235.8 mN/m`, `b = −0.625`, `μ = 1.256`; "This equation is valid between the triple point (0.01 °C) and reference temperature, Tc. It also provides reasonably accurate values when extrapolated into the supercooled region, to temperatures as low as −25 °C." | Use only if a temperature off the 5 °C table grid is needed. pdftotext drops the minus sign on `b` and the whole equation body — read the rendered page, not the text dump. |
| Q3 | Uncertainty at 20–25 °C: ±0.36 mN/m (≈0.5 %). | Sets the floor on how precisely any collapse threshold can be quoted from γ alone. |

## Use in petch

Sourced γ for `research_sources/RESEARCH_PATTERN_COLLAPSE_CRITERION_2026-08-20.md`.
Supersedes the "about 0.072 N/m" quoted by Mack 2006.
