# The grazing lip-removal law: faithful to source, and not the missing removal (2026-08-04)

`RESULTS_O_CHANNEL_2026-08-04.md` (`63cfefa`) left one suspect standing for the
top-band closure excess (10.8x -> 9.4x after the O fix): the angular ion-removal
law at the 86-89 deg mask-lip faces, where `cos^2 theta (1 + B sin^2 theta)`
collapses ~190x and delivers 0.5 % of deposition against the ~0.71 x deposition
still unexplained.

This pass audited that law against Krueger's verbatim rows, his Eq. (2.40)
semantics, and the in-chemistry measured bounds. **The law is faithful and the
hypothesis is falsified: no angular function can balance the lip.** Nothing in
the removal law was changed.

## Verdict: none of (a) mis-transcription, (b) double-cosine, (c) measured shape

### (a) Transcription is correct -- verbatim rows

From `tmp/pdfs/krueger_thesis.txt`, Appendix B, Table B.0.1 (columns are
`p0  Eth  n  Er  angle-class`):

```
CF(s)        + Ar+  ->  EP    + Ar#                 0.9    20  0.5  500  1
AC(s)        + Ar+  ->  C     + Ar#                 0.001 200  0.4  250  1
SiO2(s)      + Ar+  ->  Ar#   + SiO2                0.0852 70  1    140  1
SiO2CF(s)    + Ar+  ->  SiF   + CO2  + Ar#          0.1471 35  1    140  2
```

and the legend (line 5332):

> The reaction probability `p0` is modified according to Eq. (2.40) if angular or
> energy dependence of the reaction is present. In that case, `E0, th` and `q`
> define the energy dependence and `∠` defines the nature of the angular
> dependence, with `∠=1` corresponding to the results obtained by [1] and `∠=2`
> corresponding to the results obtained by [2].

with `B.1 References [1] = Kress et al., JVST A 17, 2819 (1999)` and
`[2] = Chang & Sawin, JVST A 15, 610 (1997)`.

So our polymer-sputter row (`p0 = 0.9`, `Eth = 20`, `n = 0.5`, `Er = 500`,
Kress angular class) matches his table exactly, and `_threshold_power_yield`
matches his Eq. (2.40) energy form
`p(E,theta) = p0 f(theta) (E^n - Eth^n)/(Er^n - Eth^n)`.
**Kress-for-polymer is his citation, not our substitution** -- the "wrong
system" objection (Kress is Cu/Ar MD) is an objection to the source mechanism,
not to our transcription of it.

### (b) The double-cosine is not a bug, and would not close the gap anyway

Huang's thesis (same MCFPM lineage, Eq. 2.32 = Krueger's 2.40) defines the
intent verbatim:

> For physical sputtering, `f(theta)` is an empirical function with a maximum at
> 60 deg, reduced probability at normal incidence and **zero probability at
> grazing incidence**.

So `f(theta) -> 0` at grazing is *intended* in the source model, exactly as our
`(1 + B sin^2 theta) cos theta` does. In a per-particle MC the areal cosine is
supplied by the transport (fewer particles per unit area on a grazing face); in
petch it is supplied explicitly by `event_flux = flux * source_area * weight /
face_area` (`boundary_transport_3d.py`). Both codes therefore carry the same two
factors. Discarding ours would raise the lip ion share only to 0.134 x
deposition -- still far below the ~0.80 x needed.

### (c) The measured shape points the wrong way

Our `B = 9.3` gives peak/normal = 4.17 at 52.6 deg. The in-chemistry
measurements bound the peak near 1.3 (Cho 2000) and 1.33 (Schaepkens 1998).
Adopting a measured-bounded shape therefore *reduces* removal everywhere off
normal and closes the top **faster** -- the opposite of what is needed. (Noted
already in `RESEARCH_LIP_CERTAINTY_2026-08-04.md`; confirmed here.)

Shape footnote: `(1 + B sin^2) cos` peaks at `cos^2 = (1+B)/(3B)`, i.e. at
54.7 deg as `B -> inf`. It can *never* reach the 60 deg Huang describes, at any
`B`. Recorded as a bounded fidelity limit, not acted on.

## The falsification, quantitatively

Removal at a lip face = (areal ion flux) x (energy yield) x `f(theta)`.
With Table-I fluxes, `J_ion = 9.6e19`, 1500 eV (energy yield saturates at 1.0):

| theta | cos | `f` (Kress) | ion removal / deposition |
|---|---|---|---|
| 86.0 | 0.0698 | 0.7153 | 0.0287 |
| 87.0 | 0.0523 | 0.5377 | 0.0162 |
| 88.0 | 0.0349 | 0.3591 | 0.0072 |
| 88.7 | 0.0227 | 0.2336 | 0.0031 |
| 89.0 | 0.0175 | 0.1797 | 0.0018 |

To balance the lip (removal = deposition, with O supplying its exact 0.1953
share) the angular function would have to be

```
f_required(88.7 deg) = 61.6
```

against our own peak of 4.17 and a measured-bounded peak near 1.3. **The
requirement is ~15x beyond the maximum of any published per-ion yield and ~264x
our value at that angle.** No re-shaping, re-fitting, or re-citation of the
angular law can close the mask-top gap. Gated in
`tests/test_lip_removal_law_audit.py`.

## Where the residual must live instead

The lip budget is a near-cancellation of two large numbers: gross deposition at
Table-I fluxes is `1.668e20` cells/m^2/s = **~6 nm/s** of film, so an open
aperture requires removal to track deposition to within ~1 %. Our channels at a
near-vertical lip face supply 0.1953 (O) + ~0.003 (ions) = **0.198**. The
missing ~0.8 cannot come from the ion channel at any angle, so it is on the
**deposition side**: either the depositor flux delivered to near-vertical lip
faces, or the effective sticking there.

That redirects at the finding the certainty pass flagged and set aside: the
measured sidewall CFx sticking (Izawa 2007, ~0.004; KIOXIA/Hiwasa 2022 reading
it as the F-rich end of a 125x range against ~0.5 C-rich) is 20-25x below the
0.1 on-polymer class both petch and Krueger's mechanism apply. Not importable as
a constant -- it is model-inverted and species-resolved -- but it is the right
*magnitude and sign* to move a 0.198 removal-to-deposition ratio to near unity,
and it is the only anchor left that is.

## New finding: three angular markers are unimplemented

Appendix B marks four ion rows with an angular class. The module applies an
angular factor to exactly one:

| row | class | petch kernel | angular factor applied |
|---|---|---|---|
| `CF(s) + Ar+ -> EP` | 1 (Kress) | `kernel_sputter` | **yes** |
| `AC(s) + Ar+ -> C` | 1 (Kress) | `kernel_ac` | no |
| `SiO2(s) + Ar+ -> SiO2` | 1 (Kress) | `kernel_bare` | no |
| `SiO2CF(s) + Ar+ -> SiF + CO2` | 2 (Chang & Sawin) | `kernel_complex` | no |

`kernel_complex` does carry a cosine-dependent ZBL deposited-energy term, which
is a different treatment, not the `∠=2` function.

**Not implemented here, deliberately.** Two of the three rows are the oxide
channels that set trench depth, which is the one headline currently near its
target; adding angular factors there changes floor and sidewall removal and must
be graded by a confirmation run rather than landed blind. Expected direction:
`∠=1` on `kernel_bare` raises removal at 30-60 deg sidewall incidence and lowers
it near grazing; the `kernel_ac` effect is bounded at <= 0.6 nm of mask over
60 s (prior estimate) so it cannot move the mask match materially.
Recorded in `tests/test_lip_removal_law_audit.py` so the gap stays visible.

## Gates

`tests/test_lip_removal_law_audit.py`, 5 tests: the exact O share (0.1953), the
falsification bound (`f_required > 50`, and the no-areal-cosine variant still
leaving > 0.6 x deposition unexplained), the monotone grazing collapse, the
54.7 deg peak bound, and the angular-marker coverage record.

The `<= 1.3x` closure gate in the directive is **not met and was not chased**:
it is unreachable through this law, which is the result. The 200-270 nm band was
not touched -- no physics changed in this pass.

## Confirmation-run spec (do not spend it on this pass)

No run is warranted until a deposition-side change exists. When one does, the
graded run is unchanged from the O-channel doc:

```
python scripts/krueger_2024_trench_pilot.py \
  --dx-um 0.01 --radiosity-backend deterministic_extruded_2d \
  --transport-device cuda:0 --surface-state-remap-backend common_refinement \
  --topology-change-policy continue_gas_cavity --surface-model mixed_layer \
  --max-wall-s 86400 --duration-s 60 \
  --mixed-layer-volatilization-yield 1.0 \
  --grazing-ion-reflection literature_v1 --output <OUT>
```

graded on `top_cd_nm` / `neck_cd_nm` / `neck_depth_from_mask_top_nm` against
38.8 nm @ 271 nm (MCFPM) and 39.0 nm @ 200 nm (SEM), primary observable the
0-50 nm closure ratio (target <= 1.3x, currently 9.4x).

## Next step

Measure the per-face depositor flux and effective sticking on the top-band lip
faces against the source flux -- the same free, local decomposition that found
the O normalisation, now pointed at the deposition side. That is the only
remaining term large enough to carry the residual.
