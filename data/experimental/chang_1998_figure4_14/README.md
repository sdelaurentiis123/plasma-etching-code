# Chang 1998 Figure 4.14 SiO2 angular yield

This dataset replays the seven filled-square oxide markers in Chang's MIT
thesis Figure 4.14 at the original 600 dpi source resolution.  The condition is
100 eV Ar+ with atomic-Cl/Ar+ flux ratio 90.  The SiO2 response is physical:
it rises to approximately twice its normal-incidence yield near 60 degrees and
falls to zero at grazing incidence.

The source PDF is version-controlled at `research_sources/chang_thesis.pdf`;
the rendered thesis pixels are checked locally and are not duplicated here.
Reproduce the checksum, pixel localization, CSV, and a visual overlay with:

```bash
python scripts/digitize_chang_1998_figure4_14.py --check \
  --overlay /tmp/chang_figure4_14_overlay.png
```

This curve supports the Ar+/Cl/SiO2 beam mechanism only.  Reusing it for Cl+
or Cl2+ is a projectile-transfer sensitivity, not a measured oxide law.
