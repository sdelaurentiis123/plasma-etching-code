# Zhu/Oxford blind square-pillar board v1

This is the frozen pre-SEM feature board for the supplied Oxford PlasmaPro
NPG80 condition: `55/5/1 sccm CHF3/SF6/O2`, `30 mTorr`, `150 W` forward RF,
`20 C`, and `1200 s`, acting on `700 nm` ALD TiO2 under `45 nm` Cr.

The square geometry is a **sensitivity prior**, not a measurement of Freddie's
held-out layout.  It uses a same-group public `400 nm` pitch and a `80--320 nm`
width sweep.  The exact target GDS, SEM, depth, and post-etch result were not
used.

## What was solved

At six relief depths from `0` to `660 nm`, the common 3-D engine constructs one
periodic square-mask cell and deterministically gathers the declared ion
energy/angular measure onto every visible TiO2 and Cr triangle.  There is no
Monte Carlo position sampling.  The zero-relief flat limit is evaluated
analytically because a touching TiO2/Cr marching-cubes junction is not an
authoritative way to recover a flat projected area.

The incident ion measure closes to `1.98e-14`.  The `20 nm` production mesh
against `10 nm` sentinels at the two width extremes changes floor transmission
by at most `0.752%`.

At `660 nm` relief, the blanket-normalized ideal-floor dose factor is:

| inferred square width | solved floor-dose range |
|---:|---:|
| `80 nm` | `0.981--0.995` |
| `120 nm` | `0.969--0.991` |
| `160 nm` | `0.956--0.987` |
| `200 nm` | `0.938--0.982` |
| `240 nm` | `0.913--0.975` |
| `280 nm` | `0.872--0.964` |
| `320 nm` | `0.793--0.941` |

The range spans the preregistered `146.5--296.2 eV` impact-energy sensitivity
and `0--0.65` collisional-tail fraction.  It is transport uncertainty, not a
fit to a profile.

## Depth result

The deliberately broad reactor/density/yield/selectivity grid gives a
mask-pinned envelope of about `191--700 nm`.  That is not a probability
interval.  Its width is the direct consequence of leaving the target TiO2
surface yield unidentified.

A separate published-process analog slice asks the narrower conditional
question: what happens if Janissen's adjacent `34.125--43.467 nm/min` TiO2 rate
and `14--18.017` TiO2:Cr selectivity transfer?  It predicts:

| inferred square width | mask pinned | controlled while Cr survives |
|---:|---:|---:|
| `80 nm` | `675.6--700 nm` | `624.1--700 nm` |
| `120 nm` | `671.6--700 nm` | `620.8--700 nm` |
| `160 nm` | `666.8--700 nm` | `616.7--700 nm` |
| `200 nm` | `660.7--700 nm` | `611.5--700 nm` |
| `240 nm` | `652.1--700 nm` | `604.2--700 nm` |
| `280 nm` | `638.6--700 nm` | `592.7--700 nm` |
| `320 nm` | `612.7--700 nm` | `570.5--700 nm` |

Across that full-factorial analog grid, only one third of the rows retain Cr
for all 20 minutes; failing rows exhaust it at roughly `14.5--18.7 min`.
Fractions of a sensitivity grid are not probabilities.  The Janissen values
come from another machine, chemistry, and TiO2 material state and are never
installed as target coefficients.

## Verdict

The major Freddie uncertainty is now localized.  Direct feature shadowing is
modest over most of the inferred board and remains bounded even at the tightest
gap.  The exact depth/profile decision is dominated by the TiO2 removal law and
Cr survival, not by an arbitrary unknown floor-transmission factor.

This receipt does **not** claim a unique SEM, sidewall angle, bowing, final CD,
or collapse probability.  Those require a target-relevant sidewall surface law
and the exact GDS.  A blanket TiO2 loss plus residual Cr measurement from the
same run would collapse the dominant uncertainty without fitting the held-out
pillar SEM.

## Reproduce

```bash
python scripts/audit_zhu_npg80_square_pillar_blind_board.py --check
pytest -q tests/test_tio2_square_pillar.py \
  tests/test_zhu_npg80_device_geometry.py \
  tests/test_zhu_npg80_square_pillar_blind_board.py
```

The complete transport snapshots, all conditional trajectories, input hashes,
and certification flags are in `audit.json`.
