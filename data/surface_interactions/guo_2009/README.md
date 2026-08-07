# Guo 2009 C4F8/Ar--SiO2 translating-mixed-layer source board

Primary source: W. Guo, *Kinetics modeling and 3-dimensional simulation of
surface roughness during plasma etching*, MIT PhD thesis (2009), handle
[1721.1/46600](https://hdl.handle.net/1721.1/46600).

`table4_1_reaction_deck.csv` is a direct visual transcription of thesis Table
4.1 (PDF p. 100). It freezes the twenty printed reactions, rate expressions,
coefficients, thresholds, and atom-countable formulas. It is a source board,
not yet an executable surface mechanism.

The four adsorption/incorporation coefficients recovered from the page image
are:

- ion incorporation: `S_I = 1`;
- atomic-F adsorption: `S_F = 20`;
- generic-neutral adsorption on a C vacancy: `S_N_on_C = 3.5`;
- generic-neutral adsorption on an O vacancy: `S_N_on_O = 1.8`.

These are coefficients in the printed bond/site-rate expressions. In
particular, `20`, `3.5`, and `1.8` are **not probabilities** and must not be
clipped to `[0, 1]` or used as per-collision sticking probabilities.

The source says the coefficients were fitted to Yin's oxide etch yields over
C4F8/Ar conditions, while its table footnote says coefficients are assumed,
calculated, or experimentally fitted and identifies physical sputtering as
calculated. The table does not label every remaining row individually.
Accordingly this deck has an L1, beam/reactor-yield-regressed evidence ceiling;
atom accounting does not make it an atomistic potential or a first-principles
reaction network.

## Visual and PIL audit

The official MIT PDF SHA-256 is
`f5c78c0089fe4104019435c6fd34e9b8f284358dda1df0101ecec54c592750d2`.
PDF pages 100--103 were rendered with Poppler and inspected at original
resolution. Pillow 12.3.0 independently opened the rasters and verified their
RGB mode and dimensions. The render and source checksums are frozen in
`source_manifest.json`; source rasters and the 18 MB thesis PDF are not
redistributed.

The visual pass also preserved three contradictions in the source instead of
silently repairing them:

1. the species footnote omits atomic F from the neutral list although reaction
   2 explicitly consumes `F(g)`;
2. the nominal-neutral list prints `C3F3` with a superscript plus, even though
   it appears under `N`;
3. Table 4.2 prints both the `-421.98` and `+95.31` physical-sputtering terms
   with `cos²(theta)`. Read literally, that polynomial is negative over most
   off-normal angles.

The executable surface law declares a single inference for the third defect:
`+95.31 cos(theta)`, which completes the descending polynomial powers six
through zero and recovers the positive sputter-like curve described by the
source. The unmodified literal expression is also implemented as an audit
function and is never used as physics. This exponent repair has not yet been
traced to an independent original fit, so any feature-angle claim must retain
that sensitivity.

The translating-layer balance was independently checked against O. Kwon's 2004
MIT ScD thesis, *Surface kinetics modeling of silicon oxide etching in
fluorocarbon plasmas* (handle
[1721.1/28353](https://hdl.handle.net/1721.1/28353)). The Guo and Kwon equation
page renders, their dimensions, and their SHA-256 hashes are recorded in
`source_manifest.json`. Neither thesis PDF nor the audit renders are
redistributed.

An executable model may resolve the other inconsistencies only through a
declared mapping with its own evidence grade.

Run:

```bash
python scripts/audit_guo_2009_table4_1.py
```

The audit checks the data and manifest hashes, all twenty row identities,
thresholds, and exact elemental closure wherever the table prints a
composition-resolved reaction. It writes
`results/curated/guo_2009_table4_1/audit.json`.
