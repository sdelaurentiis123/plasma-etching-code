# The Appendix-B angular classes: implemented, and they do NOT trim depth (2026-08-05)

`RESULTS_LIP_REMOVAL_AUDIT_2026-08-04.md` (`aca2aeb`) found three ion rows in
Krueger's Table B.0.1 carrying an angular-class marker that the module never
applied, and deliberately deferred them because "two of the three rows are the
oxide channels that set trench depth ... and must be graded by a confirmation
run rather than landed blind."

They are now implemented, gated, and **forecast before any run is spent**. The
forecast falsifies the reason they were expected to matter: with the etch front
as flat as it actually is, these factors **raise** the predicted depth rate by
~7 %, they do not trim it. The depth overshoot seen in ml19 is not this.

## The two classes, verbatim

Krueger's legend (thesis L5332) defines the marker but not the shape:

> The reaction probability `p0` is modified according to Eq. (2.40) if angular
> or energy dependence of the reaction is present. In that case, `E0,th` and `q`
> define the energy dependence and `∠` defines the nature of the angular
> dependence, with `∠=1` corresponding to the results obtained by [1] and `∠=2`
> corresponding to the results obtained by [2].

with `[1] = Kress et al., JVST A 17, 2819 (1999)`,
`[2] = Chang & Sawin, JVST A 15, 610 (1997)`.

The shapes are stated verbatim, and identically, in three theses of the same
MCFPM lineage:

**Huang** (thesis, Eq. 2.32 discussion, L2293-2296):

> For physical sputtering, `f(θ)` is an empirical function with a **maximum at
> 60°, reduced probability at normal incidence and zero probability at grazing
> incidence**. For chemically enhanced etching, `f(θ)` is **unity for normal
> incidence and angles up to 45°, with a monotonic roll-off to zero probability
> at grazing incidence**.

**Huard** (thesis, L2386-2391) — same two, plus the separability that licenses
multiplying our energy term by an angular factor:

> the MCFPM includes two angular dependent probability functions, `P(θ)`. One
> angular dependency typically has a maximum near a 60° angle of incidence, is
> less than unity at normal incidence, and drops to zero at grazing incidence,
> characteristics of physical sputtering. The second angular dependency function
> is unity at normal incidence, gradually dropping after 45° until reaching zero
> at grazing incidence, characteristic of chemical sputtering.

> the total yield of a sputtering reaction is given by `P(ε,θ) = P(ε)P(θ)`  (2.24)

**Qu** (thesis, L2698-2703) repeats the same two classes.

### Row assignment (Table B.0.1, column 5)

| row | class | petch kernel | before | now |
|---|---|---|---|---|
| `CF(s) + Ar+ -> EP` | 1 | `kernel_sputter` | applied | applied (unchanged) |
| `AC(s) + Ar+ -> C` | 1 | `kernel_ac` | **none** | class 1 |
| `SiO2(s) + Ar+ -> SiO2` | 1 | `kernel_bare` | **none** | class 1 |
| `SiO2CF(s) + Ar+ -> SiF + CO2` | 2 | `kernel_complex` | **none** | class 2 |

### Normalisation convention (load-bearing)

Both classes are normalised to `f(0) = 1`, so `p0` in the table is the
**normal-incidence** probability at the reference energy. This is the
convention the already-validated polymer row has used since the audit-corrected
campaign; mixing conventions across rows of one table would be incoherent. The
sources' "less than unity at normal incidence" for class 1 is a statement about
the peak/normal ratio (this form gives 4.17), not about `p0`.

Consequence, and the reason every prior result survives: **at normal incidence
nothing changes at all.** All 0-D validation (DEKNOB derived energy law, rung-0
Langmuir degenerate limit, selectivity, ledgers) is taken at
`cosine_incidence = 1.0` and is bitwise unchanged.

### The class-2 roll-off is `[VERIFY]`

The plateau edge (45°), both endpoints (unity at normal, zero at grazing) and
monotonicity between them are verbatim. The interpolation across the roll-off is
the minimal projected-flux choice `cos θ / cos 45°`, which introduces no
constant beyond the stated plateau edge. Chang & Sawin 1997 is paywalled
(fetch returned 403); the specific roll-off shape is therefore recorded
`[VERIFY]` while its stated properties are gated.

## Implemented shapes

```
class 1 (Kress)        f1(θ) = (1 + 9.3 sin²θ) cos θ          f1(0) = 1
class 2 (Chang&Sawin)  f2(θ) = min(1, cos θ / cos 45°)        f2(0) = 1
```

| θ | 0° | 30° | 45° | 60° | 80° | 89° |
|---|---|---|---|---|---|---|
| class 1 | 1.000 | 2.880 | 3.995 | 3.988 | 1.740 | 0.180 |
| class 2 | 1.000 | 1.000 | 1.000 | 0.707 | 0.246 | 0.025 |

## Forecast: this does not trim depth (the falsification)

`scripts/forecast_angular_classes.py` measures the real incidence-angle
distribution from the ml18 checkpoint's etch front, convolves each face with the
delivered beam (planar σ = 0.8334°, the digitised Krueger Fig-4 width), and
flux-weights the two classes.

| region | median tilt | class 1 (bare) | class 2 (complex) |
|---|---|---|---|
| floor (depth-setting) | 0.92° | **1.1445** | **1.0000** |
| wall (lateral) | 70.51° | 3.6484 | 0.7233 |
| all interface | 1.19° | 1.2067 | 0.9931 |

Both factors were 1.000 before this change, so these *are* the multiplicative
changes.

The SiO2 removal splits **49.7 % complex / 50.3 % bare** at the Krueger base on
a film-free floor (complex kernel 1.652e20 vs bare 1.671e20), giving

```
predicted depth-rate change = 0.497 x 1.0000 + 0.503 x 1.1445 = +7.3 %
```

**The expectation that motivated deferring these rows -- that they would trim
off-normal floor removal and pull an overshooting depth back toward 825 nm --
is falsified.** The etch front is nearly flat (median tilt 0.92°), both classes
are unity at normal incidence by construction, and the class-1 Kress form *rises*
away from normal (it peaks at 54.7°). So on a flat floor the change is either
neutral (class 2) or mildly enhancing (class 1), never trimming.

Where the classes do bite is the **sidewall**, at 70° median tilt: bare-oxide
sputter is enhanced 3.65x and complex removal suppressed to 0.72x. That is a
profile-shape effect (lateral etch, bowing, taper), not a depth effect.

## Status

- Implemented and gated; every previously validated normal-incidence result is
  bitwise unchanged.
- Forecast says depth rate **+7.3 %**, i.e. it moves an overshoot the wrong way.
  A confirmation run is therefore **not** warranted for the depth question.
- The remaining depth item is elsewhere. This pass closes the last unimplemented
  entry in the Appendix-B angular column; the mechanism table now has no
  unexamined angular markers.

## Gates

`tests/test_angular_classes.py`, 8 tests: `f(0) = 1` exactly for both classes;
the class-1 shape (peak 50-60°, peak/normal 4.17, zero at grazing, monotone
either side of the peak); the class-2 shape (unity through 45° verbatim,
monotone roll-off, zero at grazing); the two classes distinct off-normal;
normal-incidence steady state invariant; oblique incidence now moves the oxide
rows (the change is live); ledgers < 1e-9 across cos = 1.0/0.7/0.2/0.05; and
single-atom == scalar to 1e-12 at 65° so the per-event and mean paths stay
consistent now that three more kernels carry an angular factor.

Suite: **1152 passed, 1 skipped**.
