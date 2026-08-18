# Zhu NPG80 wall-resolved sheath/global coupling v1

The original source-geometry state assigned the same 250 eV-per-charge ion
wall energy to the powered lower electrode, grounded upper electrode, and
grounded sidewall. This rung resolves those surfaces. The powered-electrode
drop is the target-free 276 V Oxford-80 same-chemistry family anchor plus a
Maxwellian floating plasma potential; grounded surfaces see only that plasma
potential. The reactor and sheath are iterated to a voltage residual below
0.01 V without using an SEM, etch depth, or surface yield.

At the fixed point, charged-wall power falls from 37.43 to 26.29 kW/m3,
electron density rises by 1.879x, positive-ion flux rises by 1.161x to
1.452e19 m^-2 s^-1, and F thermal flux rises by 2.648x. The exact conditional
700 nm blanket-clearance requirement becomes 0.985--1.257 TiO2 formula units
per positive ion before feature attenuation.

This is a stronger physical boundary, not a unique equipment prediction. The
276 V bias is a machine-family transfer rather than a target-tool diagnostic;
90 W remains an absorbed-power sensitivity; molecular ion-neutral angular
cross sections and the TiO2/Cr surface response are not yet measured on the
target condition. The receipt therefore keeps absolute profile-depth support
false.

The O2 workbook is the CC BY-NC Song et al. dataset
`10.60893/figshare.jpr.30850013.v1`, SHA-256
`6f98ac82e169d25d0a4328b1a3703f733668539adb8141d736209d199013c860`.
The licensed source file is intentionally not committed.

Rebuild the final fixed-point state from the preceding iterate (or iterate the
plasma potential from the source-geometry state), then check the receipt:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
python scripts/run_zhu_open_reactor.py \
  --source-workbook /path/to/song_2026_o2.xlsx \
  --initial-state-json /path/to/preceding_iterate.json \
  --powered-electrode-sheath-drop-V 296.47484368556405 \
  --grounded-surface-sheath-drop-V 20.474843685564036 \
  --output results/curated/zhu_npg80_sheath_coupled_v1/central_276V.json

python scripts/audit_zhu_npg80_sheath_global_coupling.py --check
```
