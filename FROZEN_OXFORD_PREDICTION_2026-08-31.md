# FROZEN — Oxford NPG80 TiO2 blind prediction package, 2026-08-31

The pre-SEM prediction set for Freddie's run (Oxford PlasmaPro NPG80,
55/5/1 sccm CHF3/SF6/O2, 30 mTorr, 150 W forward RF, 20 C, 1200 s,
700 nm ALD TiO2 under 45 nm Cr) is complete, locally `--check`-certified,
and frozen as of this commit. The target SEM/GDS were not used anywhere;
from here they are the answer key only.

## The three frozen boards

1. `results/curated/zhu_npg80_square_pillar_blind_v1/` — reactor-dose
   board: blanket-normalized ideal-floor dose factors per width across the
   preregistered 146.5–296.2 eV / tail 0–0.65 envelope.
2. `results/curated/zhu_npg80_moving_cr_profiles_v1/` (v7, 112 endpoints,
   audit.json, check PASS) — legacy square width prior 80–320 nm:
   etch depth 684 nm (w80) falling to 627–654 nm (w320); Cr centre
   exhausted in every cell; sidewalls ~85–86.5 deg; w80 middle CD necks
   to 41–77 nm.
3. `results/curated/zhu_npg80_gds_square_profiles_v1/` (368 endpoints,
   audit.json, check PASS, envelopes rendered) — exact-GDS-inferred width
   ladder 105–250 nm: depth 684→656 nm with width; bottom CD tracks
   nominal GDS CD with a small positive bias; Cr exhausted at all widths.

Supporting frozen physics: the capillary pattern-collapse criterion
(`research_sources/RESEARCH_PATTERN_COLLAPSE_CRITERION_2026-08-20.md`)
gives intact full-height pillars an ~8–11x stiffness margin at drying, so
observed collapse implies etch-side weakening (necked waists, exhausted
caps, base damage), not capillary failure of healthy pillars.

## Claim boundary (unchanged from the board closeouts)

- The surface law is cross-machine rate-normalized (Janissen witnesses);
  self-bias/IED are interval inputs; mask polarity on the exact-GDS board
  remains operator-unconfirmed. This is a conditional blind interval
  prediction, not a certified absolute forecast of the target tool.
- `target_sem_used: false` and `coefficient_selected_from_target: null`
  in both audits.

## What happens next

1. Owner sends `ASK_FREDDIE_DATA_REQUEST_2026-08-20.md` — the self-bias
   reading and same-run blanket TiO2 loss narrow these intervals without
   touching the answer key.
2. When the SEMs unseal: score depth, CDs, sidewall angle, Cr survival,
   and collapse clustering against these boards as-is. No re-running with
   knowledge of the target.

Compute: closing box 49334941 (150 trajectories) destroyed after
retrieval; zero instances remain.
