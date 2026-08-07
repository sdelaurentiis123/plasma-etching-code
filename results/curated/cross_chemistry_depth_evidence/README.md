# Cross-chemistry absolute-depth evidence

The common scientific claim is narrower—and stronger—than “validated across
chemistries.” Two non-Krüger chemistry families have absolute surface
depth-per-dose evidence with no fit to the compared target:

- mass-selected SF5+ on Si and SiO2: four points, 5.88% MAPE, 15.04% maximum
  error, retrospective and normal-incidence only;
- Cl2/Ar+ silicon ALE: three absolute depths, 12.88% maximum relative error,
  retrospective cross-source transfer with measured fluence but inferred mean
  ion energy.

Neither is a feature-profile test. Krüger is an evaluated evolving-feature
test and currently misses absolute depth by 58% under the published aggregate
boundary.

Mahorowala 1998 adds a different and valuable class: 11 absolute-rate points
at one fixed 75 s exposure, spanning 112.5–459.375 nm dimensionally, with all
13 run-matched SEM panels in pure Cl2/oxide-mask poly-Si. The table and montage
have passed checksum-bound, original-resolution PIL audits. This fixes the old
chlorine-corpus problem of missing exposure time, but it is not yet a model
pass: the source does not measure the species-resolved wafer flux or per-run
IEAD/IAD, and we have completed zero predictions against it.

Yoshie provides 49 genuinely value-blind held-out feature-depth targets with
full visual/PIL reconciliation, but zero mechanism predictions have been
revealed or graded because the phase-resolved wafer boundary and persistent
C/F/S surface-memory closure are incomplete.

Therefore the count of completed formal held-out feature/profile predictions
is **zero**, not two. `audit.json` binds these classifications to the source
receipts and prevents a future summary from silently promoting surface
depth-per-dose or a prepared dataset into a held-out feature pass.

Replay:

```bash
python scripts/audit_cross_chemistry_depth_evidence.py
```
