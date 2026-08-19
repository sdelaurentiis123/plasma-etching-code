# Zhu NPG80 daughter-EEDF sensitivity v1

This board freezes each conserved Oxford reactor state and repeats only the
electron Boltzmann solve after adding the two highest-leverage missing neutral
targets.  It preserves the dimensional electric field and 13.56 MHz drive as
the represented collision density changes.

The result is decisive about model form:

- the literature-authorized partial HF replay lowers mean electron energy by
  27.7--45.0% across the 60--120 W board;
- the frozen electron-growth coefficient is positive in every parent-only
  replay and negative after HF is added;
- F2 changes mean energy by less than 0.7% beyond HF, but its attachment makes
  the electron-balance change more negative.

Therefore the parent-only heavy-particle state cannot be certified as the
self-consistent solution once daughter collisions are admitted.  The next
physics step is a nonlinear reactor reclose with the enlarged collision basis.
This board is not that reclose and supports no wafer flux, TiO2 depth, or SEM
claim.

HF is a deliberately partial literature replay.  Huang et al. transfer HCl
momentum to HF and threshold-shift HCl dissociation/ionization.  Their separate
HF vibrational and attachment sources are not present here, so completeness is
false.  The SIGLO F2 EFFECTIVE row is deconvolved as elastic plus the explicit
inelastic set over 0--120 eV; raw use would double-count momentum loss.

Raw LXCat bytes are not redistributed.  The receipt was generated from local
hash-locked exports found through a secondary mirror and retains the official
Hayashi/SIGLO database references.  Certification requires replacing those
files with fresh official exports and comparing selected-process hashes.

Full local replay:

```bash
python scripts/audit_zhu_npg80_daughter_eedf_sensitivity.py \
  --source-workbook /path/to/Song-2026-O2-supplement.xlsx \
  --hcl-lxcat /path/to/user-supplied-HCl-LXCat-export.txt \
  --f2-lxcat /path/to/user-supplied-F2-LXCat-export.txt
```

The committed scalar receipt can be structurally checked without the
rights-restricted source bytes:

```bash
python scripts/audit_zhu_npg80_daughter_eedf_sensitivity.py --check
```

Primary source trail:

- Huang et al., *J. Vac. Sci. Technol. A* **38**, 023007 (2020), DOI
  `10.1116/1.5125568`.
- Pitchford et al., *Plasma Processes and Polymers* **14**, 1600098 (2017),
  DOI `10.1002/ppap.201600098`.
- Hayashi and SIGLO databases on LXCat; source bytes remain user-supplied.
