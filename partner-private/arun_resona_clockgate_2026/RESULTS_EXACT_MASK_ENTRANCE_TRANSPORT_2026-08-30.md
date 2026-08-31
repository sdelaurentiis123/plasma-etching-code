# Exact-mask entrance transport and conditional etch explorer

Date: 2026-08-30
Branch: `codex/validation-first-multiphysics`

## Result

The supplied 30 um printed clock-gate mask is not one uniform aspect-ratio
feature. Its exact ten-polygon topology produces a spatially varying entrance
filter. A new deterministic operator integrates complete collision-free
characteristics through that topology without random sampling or a redrawn
proxy mask.

At the 0.4 um display grid, the exact-opening averages are:

| Incident population | Declared angular law | Mean direct transmission | 5--95% spatial interval |
|---|---|---:|---:|
| ions | 0.5 deg transverse-component sigma | 0.9603 | 0.7286--1.0000 |
| ions | 1.0 deg | 0.9249 | 0.7286--1.0000 |
| ions | 1.5 deg | 0.9165 | 0.7280--1.0000 |
| ions | 2.0 deg | 0.8814 | 0.6982--1.0000 |
| ions | 3.0 deg | 0.8082 | 0.5307--0.9991 |
| thermal neutral direct component | cosine incident-flux hemisphere | 0.06874 | 0.01860--0.11807 |

The conclusion is robust and useful: a narrow ion beam mostly clears the tall
entrance mask, while direct thermal-radical access is suppressed by roughly an
order of magnitude and varies strongly across the layout. The leading
unresolved material input is therefore the fate of F and O after encounters
with the printed-polymer walls. A scalar blanket rate cannot identify that
return/consumption law.

This is new information from Arun's actual mask rather than a transfer from a
rectangular-trench benchmark.

## Operator

The input is the unique periodic Boolean opening field from the exact GDS
polygons. For each deterministic angular ordinate, a floor node receives its
weight only if the complete straight characteristic from mask bottom to mask
top remains in gas. The path subdivision resolves at least two samples per
grid cell crossed in x, y, or z, so a one-cell track cannot be jumped.

The ion rule is tensor Gauss--Hermite quadrature over two independent signed
Gaussian transverse angles. The neutral rule is Gauss--Legendre quadrature in
polar cosine and paired midpoint azimuths for a cosine incident-flux
hemisphere. Both rules sum to one to machine tolerance.

The reusable implementation is `src/petch/extruded_mask_transport.py`. It is
not clock-gate-specific and accepts any Boolean footprint produced from GDS,
STL, analytic geometry, or another level-set source.

## Numerical checks

For the central 1.5 deg ion case:

- changing 5x5 to 7x7 Gauss--Hermite quadrature moves the mean transmission by
  0.0125 absolute;
- refining the footprint from 0.4 to 0.2 um with the same 5x5 rule moves the
  mean by 0.0146 absolute;
- the 0.2 um result gives 13 cells across the minimum 2.6 um opening.

The explorer therefore reports the display map as an exploratory transport
field and carries the numerical changes in its machine-readable receipt. It
does not promote sub-grid morphology or an atomic-accuracy claim.

## Conditional surface-law transfer

The explorer evaluates the common Belen SF6/O2 coupled F/O silicon equations
at every opening pixel. It uses the direct ion map and exposes separate F- and
O-wall recovery fractions. The open-feature rate is normalized to Miao et
al.'s independent 3.5--3.9 um two-minute result; Arun's target profile is not
used.

One central UI scenario uses 1.5 deg ions, 100 eV mean energy, 50% F wall
recovery, 10% O wall recovery, and a 3.7 um open-feature anchor. It produces a
mean conditional depth of 2.80 um with a 2.53--2.95 um 5--95% spatial interval.
Those recovery fractions are deliberately labeled interactive assumptions.
The result is a sensitivity example, not the frozen prediction.

This construction is still more physical than a hand-written recipe slider:

1. the exact mask geometry and direct transport are computed;
2. the local F/O/ion fluxes enter the actual coupled surface equations;
3. the independent process depth fixes only the open-rate normalization;
4. every unmeasured target-tool or polymer quantity remains exposed.

## What remains before an absolute 3-D profile claim

1. Add deterministic diffuse, reactive wall return for F and O. The current
   neutral map is the no-return direct component.
2. Evolve silicon and polymer sidewalls in the common material router, rather
   than treating the output as a depth heightfield.
3. Supply the target etcher's achieved ion-energy/angular distribution or DC
   self-bias, and either wafer fluxes or a same-run blanket Si depth.
4. Measure pre/post polymer height to constrain mask erosion and radical wall
   consumption.
5. Score the frozen result against a cross-section that was not used as input.

The first three measurements can be collected on one sacrificial run. They
would turn the current broad wall-recovery sensitivity into a narrow,
tool-specific profile interval.

## Reproduction

```text
PYTHONPATH=src python \
  partner-private/arun_resona_clockgate_2026/scripts/build_etch_explorer_data.py \
  --check

python partner-private/arun_resona_clockgate_2026/explorer/build.py

PYTHONPATH=src pytest -q tests/test_extruded_mask_transport.py
```

Machine-readable payload:
`results/etch_explorer_data.json`

Interactive fragment:
`explorer/arun_etch_explorer.html`
