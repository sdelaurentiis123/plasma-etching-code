# Oxford conditional-profile numerical convergence sentinel

This audit isolates timestep and grid error for one target-free conditional
profile: a 200 nm square pillar, 360 s etch, 43.4667 nm/min cross-machine TiO2
rate endpoint, and the high-energy / 0.65 angular-tail transport sensitivity.
The gates were committed in `939af19` before the refined result was revealed.

| case | mesh (nm) | maximum step (s) | depth (nm) | top CD (nm) | middle CD (nm) | bottom CD (nm) |
|---|---:|---:|---:|---:|---:|---:|
| coarse timestep | 20 | 8 | 259.95 | 196.67 | 200.29 | 235.39 |
| refined timestep | 20 | 4 | 258.58 | 196.83 | 200.78 | 239.21 |
| refined grid | 10 | 4 | 260.79 | 197.24 | 198.18 | 215.03 |

## Verdict

- Timestep depth change: 0.529%, passing the frozen 2% gate.
- Maximum timestep CD change: 3.81 nm, passing the frozen 5 nm gate.
- Grid depth change: 0.844%, passing the frozen 5% gate.
- Grid top/middle CD changes: 0.42 / 2.60 nm.
- Grid bottom-CD change: **24.17 nm, failing the frozen 10 nm gate**.
- Fine-grid maximum x/y difference: 0.00092 nm, passing the frozen 5 nm
  symmetry gate.
- Particle balance closes exactly at the reported precision; the largest
  conservative remap residual is `3.61e-15`.

The grid refinement therefore supports a stable conditional depth and stable
top/middle CDs for this sentinel. It does not certify the bottom CD or the
coarse-grid bottom flare. The fine solution still predicts physical lower-wall
flare (bottom CD 215.03 nm versus top CD 197.24 nm), but its magnitude needs a
5 nm mesh rung before it can be quoted quantitatively.

This is numerical verification of one rate-normalized conditional law. It does
not validate that law for the Oxford target, identify the missing TiO2/Cr
surface coefficients, support a unique target SEM, or support atomic accuracy.

## Reproduction

```bash
python scripts/audit_zhu_npg80_profile_convergence.py --check
pytest -q tests/test_zhu_npg80_profile_convergence.py
```

`cases/` retains every full profile and its exact execution specification;
`audit.json` binds them by SHA-256 to the frozen gates.
