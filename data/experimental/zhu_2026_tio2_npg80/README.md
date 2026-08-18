# Zhu Oxford NPG80 TiO2/Cr condition — pre-SEM record

This directory freezes the exact Oxford PlasmaPro NPG80 RIE condition supplied
by the Nanfang Yu Group on 2026-08-18 before its SEM was available. It is a
prospective validation case, not a fitted profile. The application is a
visible-wavelength TiO2 metasurface-pillar process.

## Authoritative condition

- `20 min` etch at `150 W` table-RF demand and `20 C`
- `3.0e-2 Torr` pressure
- `55 sccm CHF3`, `5 sccm SF6`, `1 sccm O2`, no Ar during etch
- `700 nm` ALD TiO2 on fused silica
- `45 nm Cr` mask

The screenshot displays `4.5e-2 Torr` and `2 sccm SF6`. The operator corrected
those values after taking the screenshot; the corrected `3.0e-2 Torr` and
`5 sccm SF6` values are authoritative. The image bytes and transcription are
checksum-bound in `recipe_manifest.json`. The meeting transcript supplies the
ALD/fused-silica stack and identifies square, rectangular, and cross-section
pillar families, but it does not identify which exact layout the Monday SEM
will contain.

The same group has a 2026 Nature paper on TiO2/SRN metasurface optical-tweezer
arrays. Zezheng Zhu is a coauthor and is credited with metasurface design and
fabrication; Teng Qu is not an author on that paper. Its TiO2 example reports
`750 nm`-tall, `100–190 nm`-wide meta-atoms on a `290 nm` unit cell and shows
vertical, defect-free pillars. Those dimensions are recorded as adjacent
device evidence only: the supplied stack is `700 nm`, and no source establishes
that the Nature devices used this exact NPG80 recipe.

## What is already determined without the SEM

The film can clear in 20 minutes only if its effective patterned TiO2 removal
rate is at least `35 nm/min`. A 45 nm Cr mask can protect a full 700 nm relief
only if the effective TiO2:Cr selectivity is at least `15.56:1`, before adding
any engineering margin.

Janissen et al. measured about `40 nm/min` and `14:1` TiO2:Cr selectivity in an
adjacent CHF3/O2 TiO2 RIE process. That comparison makes clearance plausible,
but at `14:1` the mask supports only `630 nm` of protected TiO2 removal. It is
therefore a mask-survival warning, not a transferable prediction. The supplied
condition additionally contains SF6 and runs on a particular machine whose
achieved self-bias is unknown.

Run the target-free receipt with:

```bash
python scripts/audit_zhu_npg80_tio2_pre_sem.py --check
```

## What remains unidentifiable from the recipe alone

Forward power and feed flows do not uniquely determine absorbed power,
self-bias, ion energy distribution, species-resolved wafer fluxes, or the
TiO2/Cr surface response. The exact absolute depth and profile are therefore
not yet a defensible model output. The highest-value missing observations are:

1. achieved DC self-bias (and, if available, RF voltage/current waveforms),
2. blanket TiO2 rate and Cr loss or residual Cr for this condition,
3. TiO2 phase/density and the actual wafer-surface temperature,
4. mask geometry, pitch, local loading, chip/carrier dimensions and position,
5. the scale-bearing SEM for this exact condition.

Those measurements separate reactor-boundary error, mask exhaustion, and
feature-transport error. The SEM must not be used to retune the pre-registered
receipt.
