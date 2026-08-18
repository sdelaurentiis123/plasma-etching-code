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

The full supplementary tables are now visually audited rather than represented
only by that headline. The closest stack witness uses the same `45 nm` Cr mask,
`175 nm` diameter, `200 W`, `37.5 mTorr`, and a measured `-950 V` DC bias. A
within-batch power sweep gives `30/58/68 nm/min` at `100/165/200 W`, and two
feature batches on a second nominally identical RIE give `273 nm` after `8 min`
and `652 nm` after `15 min`. Those exact rows form a target-free process and
depth validation board. They do not transfer as an Oxford/ALD-TiO2 coefficient:
the source tool is a Fluor Z401S, the material is single-crystal rutile, and the
target adds SF6. Rebuild it with:

```bash
python scripts/audit_zhu_npg80_tio2_analog_board.py --check
```

Run the target-free receipt with:

```bash
python scripts/audit_zhu_npg80_tio2_pre_sem.py --check
```

## Blind reactor-dose clearance call

The shortest reactor-to-depth ledger is now frozen separately from the larger
chemistry build.  Conditional on the central conserved reactor state's
`2.29e19 m^-2 s^-1` global axial positive-ion flux, clearing `700 nm` in
`1200 s` requires `0.62--0.80` TiO2 formula units per incident positive ion
over the published `3.25--4.15 g cm^-3` ALD-density sensitivity.  A candidate
surface yield of `1` therefore needs a run-averaged feature-floor transmission
of `0.62--0.80`; a yield of `2` needs `0.31--0.40`.  These are exact atom/dose
requirements, not fitted yields.  The global axial flux remains a sensitivity,
not a validated local wafer diagnostic.

The independent Janissen feature board maps to `682.5--700 nm` after twenty
minutes when capped at the supplied film thickness; that pair is a
cross-machine comparison, not a target confidence interval.  Its mask rows
also straddle the target: the closest feature witness supports `630 nm` with
`45 nm` Cr, while its separate power-sweep selectivity supports `811 nm`.
Hegeman et al. independently measure the expected direction of the added SF6
term (`55 nm/min` TiO2 in SF6 versus `15 nm/min` in CHF3 at their common ICP
condition), but their rates and selectivity are not transferable to this CCP.

The preregistered binary call is therefore **full TiO2 film clearance expected
(`700 nm`)**, with Cr-mask exhaustion or feature-bottom transport as the
highest-risk route to a shorter or unusable pillar.  This is deliberately not
called an atomic-accuracy profile or a unique continuous uncertainty interval.
Rebuild the receipt with:

```bash
python scripts/audit_zhu_npg80_tio2_depth_gate.py --check
```

## Measured CHF3 ion-collision floor

The NIST Peko data now supply a checksum-replayed reactive-destruction kernel
for `CF3+ + CHF3`: a 21-point summed-DCT curve plus the reported `18 A2`
summed CID cross section. The implementation explicitly converts powered-ion
laboratory energy to the paper's relative collision energy. At the target
pressure, an ideal-gas/feed-fraction scale check gives order `0.2` reactive
optical depth per millimeter near a `200 eV` CF3+ laboratory energy. That is
large enough to reject the collisionless sheath as the final boundary.

This is not yet a full molecular IEAD. The comparison uses feed fraction as a
neutral-density proxy and a declared 1 mm slab; elastic/momentum-transfer and
angular scattering, other ion-neutral pairs, and the actual sheath composition
remain open. Rebuild the target-free scale receipt with:

```bash
python scripts/audit_zhu_npg80_cf3_collision_scale.py --check
```

Basurto and de Urquijo independently measured the mass-resolved `CHF2+` reduced
mobility in `CHF3`. The checksum-pinned Figure 1 series now provides a C1,
no-extrapolation swarm closure over `45.18--450.35 Td`. At `30 mTorr`, using
the `20 C` electrode setpoint only as a declared ideal-gas density proxy, the
derived drift-relaxation scale grows from roughly `0.1 mm` at `100 Td` toward
`1 mm` at `400 Td`. This is direct evidence that molecular momentum relaxation
belongs in the boundary model; it is not a unique elastic/angular cross
section or a target IEAD. Rebuild the receipt with:

```bash
python scripts/audit_zhu_npg80_chf2_mobility_scale.py --check
```

## Measured feed-electron kinetics

The exact `55/5/1 sccm` feed is now assembled from independently audited
CHF3, SF6, and O2 electron-collision data and solved with the deterministic
finite-frequency two-term Boltzmann operator at the measured `13.56 MHz` and
`30 mTorr`.  A field scan from `40` to `400 Td` brackets the frozen-feed
attachment-to-ionization transition between solved points at `200` and
`225 Td` (a linear diagnostic interpolation gives about `210 Td`).  This
corrects the earlier `40--100 Td` development scan, which sampled only a
net-attaching state and therefore could not participate in a sustaining
power/particle balance.

The bracket is not a fitted bulk field.  The feed fractions are not the
unknown steady plasma composition: CHF3 and SF6 dissociate, their daughter
species alter attachment and ionization, and wall/exhaust loss closes the
inventory.  The complete 38-reaction hydrogen-bearing CHF3 extension from
Sandia Table 9 is conservation-checked in code, with all copied and estimated
rates still labelled as such.  It is the next chemistry rung, not yet a
species-resolved wafer-flux prediction.  Rebuild the measured-feed receipt
with the licensed Song O2 workbook using:

```bash
python scripts/audit_zhu_npg80_feed_electron_kinetics.py \
  --source-workbook /path/to/o2_song_2026_supplement.xlsx
python scripts/audit_zhu_npg80_feed_electron_kinetics.py --check
```

The SF6 parent chemistry is now product-resolved rather than an aggregate
`SFx+`/`SFx-` source.  NIST Tables 24, 25, and 27 provide seven negative-ion
curves; Tables 3 and 16 constrain nine positive-ion branches.  Their executable
split preserves the evaluated total attachment and ionization curves exactly
and treats double ionization with the correct two-electron source.  The entire
three-gas deck is mapped through 46 atom- and charge-checked events in
`zhu_parent_collision_chemistry.py`; its shared species basis contains 45
entries including the daughter F2+ closure. Two source gaps
remain visibly graded: aggregate SF6 neutral dissociation currently uses the
literature-dominant `SF5 + F` branch, and aggregate O2 ionization currently
uses the `O2+` branch.  This closes the parent collision source term, not the
daughter-species, wall, exhaust, or absorbed-power balances. The supplemental
reactor network now contains 259 conserved daughter/heavy reactions over a
66-species basis, including the Sandia H/H2 loop, the Lim CHF3/O2 neutral
chain, and the Pateau SFx/O -> SOxFy titration chain that releases F,
while keeping target-machine wall physics and cross-family ion closures
explicit.

## First conserved open-reactor state

The fixed-pressure model now solves all 66 species, quasineutrality, throttle
loss, the represented-feed electron state, and absorbed-power balance in one
open system.  The central target-free sensitivity (`90 W` absorbed, `350 K`,
`170 mm` powered diameter, `30 mm` active height, published-compilation wall
probabilities) converges to a maximum normalized balance residual of
`3.90e-7`.  It produces `2.29e19 m^-2 s^-1` total axial positive-ion flux and
explicit thermal neutral-flux channels without using the held-out SEM or any
TiO2 depth.

That convergence exposes rather than hides the next gap.  Only `37.6%` of the
neutral inventory remains in the three gases represented by the exact
Boltzmann collision deck; HF and daughter molecules dominate the solved gas.
The reported `311.9 Td` is therefore `E/N` on the represented-parent basis.
Holding its dimensional field fixed gives `117.1 Td` on the total-neutral
basis, but neither number closes electron momentum and energy loss in the
unrepresented daughter gas.  The state is a conserved chemistry checkpoint,
not yet a certified wafer boundary.  Rebuild it with a previously converged
continuation state using:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
python scripts/run_zhu_open_reactor.py \
  --source-workbook /path/to/o2_song_2026_supplement.xlsx \
  --initial-state-json /path/to/prior-open-reactor-state.json \
  --output results/curated/zhu_npg80_open_reactor_v1/central.json \
  --maximum-evaluations 150
```

## Machine-family self-bias evidence

The two downloaded files are valid despite their generic/misleading names.
The NCSU thesis contains an exact Oxford PlasmaPro NGP80 CHF3 procedure whose
self-bias typically drifts from above 300 V at run start to below about 200 V
at run end. A separate Oxford PlasmaLab 80 thesis measures 276 V for the same
active-gas set (`SF6/O2/CHF3`) at exactly the target recipe's reported
`power/pressure = 5 W/mTorr`. Same-family measurements at adjacent chemistries
extend to 360--387 V at 30 mTorr and 400 V for CHF3/Ar.

Those values remove the fiction that the voltage is wholly unconstrained, but
they do not identify one exact target voltage. The committed transfer therefore
provides five deterministic sensitivity histories, preserves the exact-tool
start/end statements as censored inequalities, and forbids selecting a history
after viewing the held-out SEM. Rebuild it with:

```bash
python scripts/audit_zhu_npg80_self_bias.py --check
```

## What remains unidentifiable from the recipe alone

Forward power and feed flows do not uniquely determine absorbed power,
self-bias, ion energy distribution, species-resolved wafer fluxes, or the
TiO2/Cr surface response. The family evidence now supports a physics-constrained
voltage sensitivity ensemble, but exact absolute depth and profile are not yet
a defensible model output. The highest-value missing observations are:

1. achieved DC self-bias (and, if available, RF voltage/current waveforms),
2. blanket TiO2 rate and Cr loss or residual Cr for this condition,
3. TiO2 phase/density and the actual wafer-surface temperature,
4. mask geometry, pitch, local loading, chip/carrier dimensions and position,
5. the scale-bearing SEM for this exact condition.

Those measurements separate reactor-boundary error, mask exhaustion, and
feature-transport error. The SEM must not be used to retune the pre-registered
receipt.
