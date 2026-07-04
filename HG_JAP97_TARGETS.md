# HG JAP 82, 566 (1997) "Aspect-ratio-dependent charging" — digitized gate targets + W2 diagnosis

Read the full paper 2026-07-03 (Caltech record je8bd-j6v68, HWAjap97b.pdf) after the fine-AR
sweep showed a constant +0.10 floor-flux offset. The paper reframes the whole gate. This file is
the literature anchor for W2.

## What Fig. 4 actually plots (the thing petch gates against)

The "Normalized ion flux at trench bottom" curve (0.59 at AR1 -> 0.22 at AR4) is **NOT ion
transmission**. Verbatim: *"this curve is identical to the one of Fig. 3 that describes the
electron current to the trench bottom."* The bottom SiO2 is an insulator, so at steady state its
net current is zero and **ion flux to the bottom = electron flux to the bottom**. The curve is
therefore set by **electron shadowing** (isotropic electrons blocked by the long insulating PR
sidewalls + repelled by the negative entrance potential), not by ion optics.

Consequence for petch: floor flux is high (0.72->0.33 vs 0.59->0.22) because petch's **electron**
flux to the floor is high (measured 0.71->0.30 vs 0.59->0.22) — petch under-shadows electrons.
petch DOES enforce the balance (i_floor 0.723~=e_floor 0.710 at AR1; 0.354~=0.301 at AR4), so the
fix is to shadow/repel electrons correctly, which then pulls ion floor flux down with it.

## Digitized targets (all vs AR 1.0 -> 4.0)

Fig. 3 — normalized ELECTRON flux components (new sub-gates, sharper than floor-flux alone):
| component        | AR1  | AR4  | trend |
|------------------|------|------|-------|
| Total into trench| 0.92 | 0.82 | slow decline |
| Trench bottom    | 0.59 | 0.22 | steep decline (== Fig.4 ion) |
| PR sidewalls     | 0.20 | 0.45 | rises |
| poly-Si outer    | ~0.20| ~0.18| ~flat |
| poly-Si inner    | ~0.03| ~0.02| ~flat, tiny |

Fig. 4 — ion flux + energy at the inner poly-Si (foot) of the edge line:
| quantity                 | AR1  | AR4  | trend |
|--------------------------|------|------|-------|
| trench-bottom ion flux   | 0.59 | 0.22 | == electron |
| inner-poly (foot) ion flux| ~0.13| ~0.13| ~constant |
| avg foot ion energy (eV) | ~10  | ~28  | rises, no saturation |

Fig. 6 — poly-Si equipotentials (THE key W2 target — the split):
| line             | AR1 | AR2 | AR3 | AR4 | note |
|------------------|-----|-----|-----|-----|------|
| edge line        | 2   | 4   | 6   | 7   | stays LOW (outer wall fed by open-area electrons) |
| neighboring line | 6   | 20  | 31  | 39  | rises HARD, esp AR>2 (starved) |

Fig. 5 / text — bottom-center potential: 8 V (AR1) -> 33 V (AR4).

## The W2 defect, precisely

petch (post-W0, edge_open_model="none", AR4): Vedge 19.3, Vneigh 22.5 — nearly EQUAL. HG wants
Vedge 7, Vneigh 39 — a 32 V split. petch has no split because it does not deliver the open-area
electron supply to the **edge line's outer (open-facing) sidewall**. In petch electrons are
launched at the mouth plane moving DOWN, so they never reach a vertical outer sidewall that faces
the open area to the left; e_edgeOuter is only 0.037-0.053 (should keep the edge line pinned low).
HG mechanism (verbatim): *"The outer sidewall of the edge line is supplied by electrons from the
open space, which help prevent a significant increase in the corresponding potential. The same is
not possible for the poly-Si sidewalls of the neighbouring line, whose potential should increase
much more than that of the edge line."*

### W2 fix (this is the real task #45)
1. Model the open-area electron supply to the edge line's OUTER sidewall as a physical ballistic
   source from the open half-space (isotropic electrons arriving from the plasma above the open
   area, some with horizontal velocity toward the wall), NOT the mouth-plane down-going launch.
   The `line_of_sight` boundary was the right idea but it pinned Vedge to 0 (over-fed) and did not
   let the neighbor rise. Replace it with traced open-area electrons hitting the outer wall.
2. The neighbor line must stay starved so it rises to ~39 V. Verify it is not receiving stray
   open-area electrons through the periodic/buffer boundary.
3. Gate the SPLIT: Vedge 2->7 (+-30%) AND Vneigh 6->39 (+-30%) simultaneously (Fig. 6). Then the
   electron-flux components (Fig. 3) and floor flux (Fig. 4) should follow.

### W2 attempt log (so the next pass doesn't repeat dead ends)

- **FALSIFIED (2026-07-03): raise the electron launch plane to the sheath edge (z=1).** Hypothesis:
  open-area electrons launched from the top would illuminate the edge line's tall outer sidewall.
  Measured (AR1, W16, n4000/it400): outer-wall electron flux stayed 0.036 (unchanged, target 0.20);
  all potentials blew up to ~40 V and became equal (Vedge 39.2, Vneigh 40.7); residual 0.381 (broke
  convergence). Reason: the cos^0.6 MC launch is downward-biased, so few electrons ever move
  horizontally enough to strike a vertical wall regardless of launch height. Under-collection on
  vertical walls is intrinsic to the down-going MC source, not a launch-height artifact. Reverted.
- **Correct approach (next):** the open-facing outer wall should receive the *isotropic open
  half-space* electron flux directly, applied as a LOCAL per-cell boundary current on the outer-wall
  cells, self-consistently reduced by the local wall potential (electron flux ~ exp(V_wall/Te) for
  V_wall<0, saturated for V_wall>=0). The existing `edge_open_model="line_of_sight"` had the right
  physics (analytic open-area flux to the poly outer band) but applied it as a SCALAR lump on the
  whole edge conductor via `edge_boundary_net`, which over-fed and pinned Vedge=0. Reimplement it as
  a per-cell outer-wall flux term inside the conductor/insulator charge balance, so the edge line's
  outer wall is held low locally while its inner wall and the neighbor line stay starved and rise.

### New gates unlocked by reading the paper
- Fig. 3 electron-component gate: petch's diag["electron"] already returns floor/edge_outer/
  edge_inner/neighbor; add PR-sidewall electron flux and gate all five vs the table above.
- Fig. 6 split gate (above) — the single most diagnostic quantity.
- Fig. 4 foot energy 10->28 eV and foot flux ~0.13-constant (petch already tracks these).

Source: G. S. Hwang and K. P. Giapis, "Aspect-ratio-dependent charging in high-density plasmas,"
J. Appl. Phys. 82, 566 (1997). PDF: authors.library.caltech.edu/records/je8bd-j6v68
