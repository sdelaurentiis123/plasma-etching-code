# Hong 2023 TiO2 feature-response board

This board records the independent nanoscale TiO2 measurements in Hong et al.,
*Materials Science in Semiconductor Processing* **164**, 107617 (2023), DOI
`10.1016/j.mssp.2023.107617`.

The experiment used a 300 mm ICP etcher at 1000 W source power, 150 W bias,
2 mTorr, and 70 C. It reports blank and 80 nm-pattern TiO2/ACL rates,
TiO2-to-ACL selectivity, and the P1/P3 depth ratio for three continuous or
asynchronously pulsed modes. Figure 3 separately prints the remaining ACL and
P1 TiO2 depth on the cross-section SEMs. The source PDF and pixels are not
redistributed; their checksums and reproducible render/crop recipes are in the
manifest.

The strongest transfer-safe result is the feature response: under continuous
C4F8/SF6/Ar, the 80 nm P1 pattern etches about 1.48 times as deep as the 220 nm
P3 pattern. Asynchronous pulsing moves that ratio to about 1.09. Pattern scale
can therefore produce spatially clustered outcomes under one nominal reactor
condition. That does **not** establish a radial flux gradient.

This is an external mechanism/transport target for the Zhu work, not an
absolute calibration. Hong used an ICP, C4F8/SF6/Ar, 2 mTorr, ACL, and circular
patterns. Zhu uses an Oxford capacitively coupled RIE, CHF3/SF6/O2, 30 mTorr,
chromium, and square/rectangular metasurface pillars. No Hong absolute rate or
surface coefficient may be transplanted into the Zhu deck.

Regenerate the committed tables and manifest with:

```bash
python scripts/digitize_hong_2023_tio2.py --write
python scripts/digitize_hong_2023_tio2.py
```
