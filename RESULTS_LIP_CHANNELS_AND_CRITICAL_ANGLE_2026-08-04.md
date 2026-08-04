# The lip budget by channel, and the critical wall angle (2026-08-04)

Track A's remaining defect after the O-channel fix (`63cfefa`), the removal-law
audit (`aca2aeb`) and the probe transport fix (`cbbd2d6`) is the mask-top band:
the evolution runs close it 10.8x faster than Krüger's profile
(`RESULTS_NECK_REGRADE_2026-08-04.md`), and the corrected probe reads
removal/deposition = 0.216 there against 0.88-0.99 in the bands that match.

This pass measures two things that had never been measured with correct
transport: **which channel supplies the removal** (by ablation, not by
re-deriving laws), and **how the balance depends on wall tilt** (the angle sweep
the previous pass ran on the truncated neutral field, invalidated by `cbbd2d6`).

## Prediction, recorded before the measurements returned

From the per-band audit table (`RESULTS_LIP_DEPOSITION_AUDIT_2026-08-04.md`) the
O share of deposition is geometry-free at 0.195 (isotropy 1.0000 face by face),
so the ion share is `removal/deposition - 0.195`:

| band | tilt (deg) | removal/dep | implied ion share |
|---|---|---|---|
| 0-50 | 0.47 | 0.216 | 0.021 |
| 50-100 | 2.24 | 0.513 | 0.318 |
| 100-150 | 6.78 | 0.993 | 0.798 |
| 200-270 | 4.79 | 0.884 | 0.689 |

Both petch ion-removal factors carry `cos(incidence)` -- the areal flux and the
angular yield -- so on a wall tilted `a` from vertical the ion removal scales as
`sin^2(a) (1 + 9.3 cos^2(a))`, a function that is 6.9e-4 at a = 0.47 deg and
0.142 at a = 6.78 deg: a 205x lever over the range the mask lip actually
occupies. Normalising on the top band (where the ion flux is largest) gives a
balance at

    removal/deposition = 1  ->  a_crit ~ 2.9 deg

which is a *sub-cell* tilt at dx = 10 nm (a 2.9 deg wall displaces 0.5 nm across
one cell). The prediction this sets up:

* the corrected angle sweep should show a **zero crossing near 3 deg**, not the
  "no crossing at any angle" the truncated-neutral sweep reported;
* Krüger's own digitised profile should **bracket** that angle -- his mask-top
  bands at 7.3 and 9.6 deg (removal dominated, stays open) and his neck band at
  1.86 deg (deposition dominated, closes) -- which is exactly the shape of his
  profile;
* our evolved profile brackets it too, but crosses **higher up** (17.3 deg at
  0-50 nm, 2.2 deg by 100-150 nm), which is exactly where our neck forms.

If that holds, the lip chemistry balance is not wrong: the difference between
the two profiles is *where the tilt crosses the critical angle*, which is a
profile-evolution and resolution question rather than a missing mechanism.

## Measurement 1: the critical wall angle, from the audit's own per-band data

Solving `ion_share(a_crit) = 1 - 0.1953` band by band, with each band's ion
share normalised on its *own* measured value (so the ion flux attenuation with
depth is carried by the data, not assumed):

| band (nm) | measured tilt | removal/dep | implied ion share | shape factor | a_crit |
|---|---|---|---|---|---|
| 0-50 | 0.47 | 0.216 | 0.021 | 7.0e-4 | **2.93** |
| 50-100 | 2.24 | 0.513 | 0.318 | 1.58e-2 | **3.58** |
| 100-150 | 6.78 | 0.993 | 0.798 | 1.42e-1 | **6.81** |
| 200-270 | 4.79 | 0.884 | 0.688 | 7.12e-2 | **5.18** |

Two things follow immediately.

**The balance is a critical-angle phenomenon, not a magnitude deficit.** The
shape factor spans 205x across the 0.5-7 deg range the mask lip occupies, so
removal/deposition passes through 1 within that range at every depth. The
100-150 nm band sits at tilt 6.78 deg against a_crit = 6.81 deg and reads
removal/deposition = 0.993 -- the probe geometry is, by coincidence, sitting
exactly on its own balance point there. Nothing about that is a tuned number:
it falls out of two published laws (`p_ox` on an isotropic O flux, and the
angular sputter yield on a directional ion flux) plus geometry.

**a_crit rises with depth** (2.9 -> 3.6 -> 6.8 deg over the first 150 nm)
because the ion flux attenuates faster than the neutral flux does. A wall that
is stable near the top is unstable further down unless it keeps tilting.

## Measurement 2: both profiles' neck locations are predicted by a_crit

Applying the criterion `tilt > a_crit -> opens` to the two measured profiles
(Krüger's digitised MCFPM profile and our evolved ml16a checkpoint, tilts from
`RESULTS_WALL_SLOPE_FALSIFICATION_2026-08-04.md`):

| band (nm) | a_crit | Krüger tilt | verdict | ml16a tilt | verdict |
|---|---|---|---|---|---|
| 0-50 | 2.93 | 7.29 | opens | 17.31 | opens |
| 50-100 | 3.58 | 9.59 | opens | 10.56 | opens |
| 100-150 | 6.81 | 5.96 | closes (marginally) | 2.20 | closes (hard) |
| 200-270 | 5.18 | 1.86 | closes | - | - |

The criterion reproduces the qualitative shape of *both* profiles: an open,
tapering top and a closing region beginning near 100 nm. It also reproduces the
difference between them, which is one of *margin*, not of sign:

* Krüger's 100-150 nm band sits 0.85 deg **below** its critical angle -- it
  closes slowly, over a broad depth range, and his neck ends up broad
  (~240 nm axially) with its minimum at 271 nm.
* Ours sits 4.6 deg below -- it closes hard and locally, and our neck is a sharp
  constriction at 120-130 nm.

So the lip chemistry is not missing a remover. What differs is **how far the
wall tilt has developed by 100-150 nm**: his profile still carries ~6 deg there,
ours has collapsed to 2.2 deg. Our taper is front-loaded (17.3 deg in the top
50 nm against his 7.3 deg) and exhausted by 100 nm.

That is a statement about the *depth profile of ion delivery to the wall in the
first 150 nm* -- a transport and profile-evolution observable -- and it is the
honest remaining item. It is not reachable by changing a chemistry constant:
every constant is depth-uniform, and a depth-uniform change moves the matched
100-270 nm bands out of agreement (they currently sit at 0.88-0.99).

## Measurement 3: the energetic populations per band (candidate (i) closed)

One corrected-transport gather, energetic flux and flux-weighted mean incidence
cosine split by population on the lateral mask faces:

| band (nm) | tilt | ion flux (m^-2 s^-1) | hot-neutral flux | hot/ion | cos ion | cos hot |
|---|---|---|---|---|---|---|
| 0-50 | 0.47 | 1.285e18 | 1.538e15 | **0.001** | 0.0226 | 0.0083 |
| 50-100 | 2.24 | 4.491e18 | 6.720e17 | 0.150 | 0.0516 | 0.0275 |
| 100-150 | 6.78 | 1.407e19 | 2.621e18 | 0.186 | 0.1268 | 0.0498 |
| 200-270 | 4.79 | 8.497e18 | 1.538e16 | 0.002 | 0.1262 | 0.0132 |

* **Hot neutrals cannot be the missing remover.** They are 0.1 % of the ion flux
  at the mask top and 15-19 % mid-mask -- and they are already counted in the
  removal the probe reports. Their distribution is physically sensible: the
  cascade creates them where ions graze the upper wall and they travel *down*,
  so they land at 50-150 nm rather than at the top. Candidate (i) is closed.
* **The mask top receives the least ion flux of any band, by 11x.** This is not
  a transport defect, it is the same cos(incidence) projection that drives the
  critical angle: a wall at 0.47 deg presents sin(0.47 deg) = 0.008 of its area
  to a vertical beam. Together with `cos ion` tracking tilt across the bands
  (0.023 -> 0.127), the transport is behaving exactly as the geometry demands.

The mask top is therefore *structurally* unable to clean itself: a vertical wall
receives no ion flux, and only the taper propagating down from the eroding
convex corner can lift it above the critical angle.

## What this makes the remaining item

Taper *propagation rate*, and there is a suggestive correlation in the archived
runs. Wall tilt at 100-150 nm, the band whose margin decides where the neck
lands:

| run | beam | tilt 0-50 | tilt 50-100 | tilt 100-150 | neck depth |
|---|---|---|---|---|---|
| ml13 | pre-P1a (narrow) | 14.28 | 11.99 | 4.64 | 228 nm |
| ml16a | post-P1a (sqrt-2 wider) | 17.31 | 10.56 | 2.20 | 130 nm |
| Krüger | his own IEAD | 7.29 | 9.59 | 5.96 | 271 nm |

The wider beam eroded the top *more* (17.3 vs 14.3 deg) and left the 100-150 nm
band *less* tilted (2.2 vs 4.6 deg), moving the neck up from 228 to 130 nm and
away from Krüger's 271 nm. Krüger's profile is the least front-loaded of the
three.

This is a correlation across two runs, not a demonstration, and it points at a
question P1a's own gate could not settle: the lift assumes the published IEAD
angle is a *planar marginal* of an axisymmetric distribution (in which case the
polar spread must exceed it by sqrt(2)). If the published angle is already the
polar angle with respect to the wafer normal, the lift is over-wide by that same
sqrt(2). P1a gated self-consistency -- the lifted planar marginal reproduces the
published width -- which holds under either reading, so it does not discriminate.
Settling it needs the HPEM convention, not another feature run -- and the
digitised source settles it **against** the over-wide reading: the published
table is `signed_angle_deg` spanning **-2.86 to +2.85 deg**, symmetric about
zero (`data/experimental/krueger_2024/digitized_figure4_iead.csv`, 878 rows,
438 negative). A polar angle with respect to the wafer normal is unsigned by
construction; a signed angle symmetric about zero is a planar projection. P1a's
reading is therefore supported by the data's own convention, the sqrt(2) lift
stands, and the front-loaded taper needs a different explanation than beam
width. The correlation above is recorded because it is real, not because it
survives this check.

## Initial geometry: settled

Krüger's initial feature is vertical-walled, so "initialise the shoulder from
the SEM" is **not** a faithful option and is formally withdrawn
(`RESEARCH_LIP_CERTAINTY_2026-08-04.md` rank 1). Verbatim, JVST A 42, 043008
(2024), Sec. IV: *"A SiO2 substrate is covered by a 850 nm thick AC film
patterned to ideally yield a straight walled opening with an initial width of
90 nm. The etch was performed for 60 s."* Both models therefore start from the
same vertical wall and must *develop* their taper; the taper's rate of
development, not its initial value, is what differs.
