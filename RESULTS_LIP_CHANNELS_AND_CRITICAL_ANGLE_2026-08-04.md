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

## Initial geometry: settled

Krüger's initial feature is vertical-walled, so "initialise the shoulder from
the SEM" is **not** a faithful option and is formally withdrawn
(`RESEARCH_LIP_CERTAINTY_2026-08-04.md` rank 1). Verbatim, JVST A 42, 043008
(2024), Sec. IV: *"A SiO2 substrate is covered by a 850 nm thick AC film
patterned to ideally yield a straight walled opening with an initial width of
90 nm. The etch was performed for 60 s."* Both models therefore start from the
same vertical wall and must *develop* their taper; the taper's rate of
development, not its initial value, is what differs.
