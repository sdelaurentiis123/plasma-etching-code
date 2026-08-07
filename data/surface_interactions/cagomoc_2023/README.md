# Cagomoc 2023 CF3 radical/ion SiO2 mechanism board

Primary source: C. M. D. Cagomoc, *Molecular Dynamics Simulation of SiO2
and SiN Etching for 3D NAND Memory Device Applications*, Osaka University
doctoral dissertation (2023), DOI
[10.18910/91922](https://doi.org/10.18910/91922).

`figure5_10_cf3_radical_ion_yield.csv` is a replayable digitization of
Figure 5.10. It records steady-state Si removal from flat SiO2 under
alternating normal-incidence 2000 eV CF3+ impacts and CF3 radical injections.
The yield is strongly non-monotone: an intermediate radical supply nearly
doubles removal, while the 200:1 and 300:1 cases form fluorocarbon film and
etch stop.

This is a **classical-MD mechanism constraint**, not an experimental
validation board and not a Krueger boundary measurement. It cannot select
Krueger's ion mixture, flux normalization, IEAD, or target depth. Its proper
use is stricter: an atom-balanced surface law for fluorocarbon oxide etching
must be capable of synergy, finite mixed-layer storage, product escape, and
polymer-driven etch stop without changing the underlying atom ledger.

## Visual/Pillow audit

The official 174-page PDF was checksum-pinned but is not redistributed. PDF
page 121 was rendered at 400 dpi and inspected at original resolution.
Pillow 12.3.0 independently verified the RGB image and its 3311 x 4682 pixel
dimensions. Dark-pixel ticks pin the yield axis; saturated-blue
distance-transform cores pin the unobscured marker centers. The two zero-yield
markers overlap the dashed zero guide, so the guide center is used and the
caption/body statement of etch stop is an independent check. All hashes,
coordinates, transforms, uncertainty, and scope limits are in
`source_manifest.json`.

PDF pages 117, 118, and 120 (Figures 5.7--5.9) were also inspected and hashed.
They constrain etch-stop versus dose, flat-versus-hole yield loss, preferential
oxygen removal, product identity, and redeposition. They have not been
digitized here because Figure 5.10 is the directly actionable surface-response
board.

## Source limitations retained

- The ion is neutralized at impact in the MD representation.
- Radicals use 0.5 eV normal incidence for computational speed, not a thermal
  isotropic feature distribution. Appendix D reports a similar deposited film
  for 0.026 eV radicals, but does not validate feature transport.
- Non-covalently bound products are removed at finite injection-cycle
  boundaries. The source explicitly warns that slow in-hole redeposition may
  therefore be underestimated.
- The curve is flat-wafer, 2000 eV CF3+/CF3, and cannot be transferred
  numerically to C4F6/Ar/O2 without an independently identified mapping.

Replay and optionally draw a QA overlay:

```bash
python scripts/digitize_cagomoc_2023_fig5_10.py --check
python scripts/digitize_cagomoc_2023_fig5_10.py \
  --check --overlay /private/tmp/cagomoc_fig5_10_overlay.png
```
