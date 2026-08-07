# Tinacba 2021 SF5+ atomistic depth-per-dose board

This board tests an atomistic surface provider against an independently
measured, fully identified beam boundary. Tinacba et al. calculated SF5+ impact
on Si and SiO2 with a DFT-informed modified-Stillinger–Weber MD potential, then
compared the result with mass-selected SF5+ bombardment. Ion energy was checked
with an energy-mass analyzer, dose with a Faraday cup at the sample position,
and depth with a contact profilometer. The MD potential was not fitted to beam
depth or yield.

At the four exact MD/experiment overlap points (Si and SiO2, 150 and 2000 eV),
mean absolute depth-per-dose error is 5.88%; the largest error is 15.04%, for Si
at 150 eV. No pass threshold was chosen after seeing the curve, so the board
reports the errors rather than manufacturing a PASS label.

The depth statement follows directly from the source equation `Y=dN/D`.
At common dose and material number density, predicted/measured depth equals
MD/measured yield. The audit also reports illustrative nanometers per
`1e16 cm^-2`, using model-film densities derived only from the source's printed
Figure-10 MD slopes. That conversion cannot change any error ratio.

This is stronger surface evidence than a fitted reactor profile, but narrower
in scope. It validates a normal-incidence SF5+ provider, not an SF6 plasma,
neutral-F synergy, angular response, or feature transport. The source
intentionally suppresses S-S, S-Si, and S-O chemistry; this must remain visible
whenever the board is cited.

Replay:

```bash
python scripts/audit_tinacba_2021_sf5_depth.py --check
```
