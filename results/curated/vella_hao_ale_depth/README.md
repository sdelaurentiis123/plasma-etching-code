# Vella–Hao Si/Cl₂/Ar⁺ absolute-depth board

This board tests a no-depth-fit transfer from the released Kounis-Melas DeepMD ALE trajectory to the
Vella–Hao experimental etch-per-cycle curve.

Regenerate the committed audit from the pinned DOE manuscript:

```bash
python scripts/audit_vella_hao_ale_depth.py \
  --paper-pdf /path/to/sha256-789bf503-vella-hao.pdf \
  --output results/curated/vella_hao_ale_depth/audit.json
```

The calculation uses exact atom-counted completed-cycle increments for the first 1000 simulated Ar
impacts, then the independently released bare-Si DeepMD sputter table for the remaining measured
positive-ion fluence. No measured depth sets a model coefficient or normalization.

The nominal 60/80/100 eV predictions are 0.4842/0.9453/1.4006 nm/cycle, versus
0.5558/1.0453/1.5427 nm/cycle from piecewise-linear interpolation of the checksum- and
vision-audited experimental curve. Maximum nominal relative error is 12.9%.

This is not a blind Tier-A result. The experimental mean energy was inferred, the positive-ion flux
was not species resolved, and no IEAD was measured. The committed audit lists these missing
uncertainties and the physical-tail assumption explicitly. It also records why the legacy transient
ROM cannot be used as an atom-conservative alternative.
