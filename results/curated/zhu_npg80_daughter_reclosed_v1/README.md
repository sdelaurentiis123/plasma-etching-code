# Zhu NPG80 daughter-reclosed reactor/sheath board

This directory is the first nonlinear reclosure of the prospective Oxford
PlasmaPro NPG80 TiO2 condition after adding explicit electron collisions with
the dominant daughter gases.  It uses no SEM, TiO2 depth, blanket etch rate, or
wafer-flux target.

The solver closes 67 species, fixed pressure, exhaust, quasineutrality,
electron power, a two-term finite-frequency EEPF, multi-ion wall transport,
and a Maxwellian floating-potential sheath fixed point at four absorbed-power
sensitivities.  Every final state is replayed from the retained preceding power
node rather than from a temporary fixed-point iterate.

## Certified board

| absorbed power | represented collision-target fraction | total-neutral E/N | positive-ion flux | F thermal flux | global residual | sheath residual |
|---:|---:|---:|---:|---:|---:|---:|
| 60 W | 84.50% | 215.82 Td | `1.0221e19 m-2 s-1` | `1.0722e21 m-2 s-1` | `2.06e-9` | `1.59e-4 V` |
| 90 W | 81.28% | 199.84 Td | `1.0576e19 m-2 s-1` | `4.2205e21 m-2 s-1` | `2.09e-9` | `1.65e-4 V` |
| 105 W | 79.12% | 190.29 Td | `1.0263e19 m-2 s-1` | `6.6882e21 m-2 s-1` | `8.25e-9` | `1.19e-4 V` |
| 120 W | 76.91% | 182.15 Td | `1.0378e19 m-2 s-1` | `9.2210e21 m-2 s-1` | `3.10e-8` | `6.01e-5 V` |

The dimensional power ledgers also close.  The largest absolute difference
between absorbed power density and parent-collision, daughter-collision, and
charged-wall losses is `9.05e-4 W m-3` on an `8.84e4 W m-3` balance.

Adding HF and F2 collisions removes the old high-power coordinate pathology.
The collision deck now covers 76.9--84.5% of the neutral inventory instead of
shrinking to 13.0% in the old parent-only 120 W state.  The represented E/N is
therefore no longer inflated to 690 Td to describe an unrepresented daughter
gas.  The remaining missing fraction is explicit and bounded, not silently
folded into a parent-gas field.

The physical response is not a simple ion-current scaling.  From 60 to 120 W,
the total positive-ion flux changes by only 1.5%, while neutral F rises by a
factor of 8.60.  The centered ion mixture also evolves: CF3+ falls from 41.7%
to 20.3%, HF+ rises from 32.3% to a mid-board maximum near 47%, and H+ grows
from 2.2% to 17.9%.  A species-agnostic ion yield is therefore not an adequate
surface boundary.

## Authority boundary

The `276 V` self-bias magnitude is a same-family transfer anchor, not a
measurement on Freddie's target run.  The `60--120 W` values are absorbed-power
sensitivities corresponding to 40--80% of the reported 150 W forward demand;
the target coupling efficiency was not measured.  The HF collision set still
lacks vibration and dissociative attachment, and several heavy-particle rates
remain literature-family transfers.  Consequently this board supports a
conserved reactor sensitivity and a species-resolved wafer lift, not one unique
machine state or SEM profile.

Rebuild from licensed/local source payloads with:

```bash
python scripts/run_zhu_npg80_daughter_reclosed_board.py \
  --source-workbook /path/to/o2_song_2026_supplement.xlsx \
  --hcl-lxcat /path/to/lxcat_ne_xe_hcl.txt \
  --f2-lxcat /path/to/lxcat_siglo.txt
```

The raw LXCat and licensed O2 source bytes are deliberately not redistributed.
Their hashes are recorded in `audit.json` and every state.
