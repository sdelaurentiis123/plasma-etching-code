# Oxford 5 nm profile-refinement verdict

This is the target-free refinement of the preregistered `200 nm`-wide,
`360 s`, high-energy/broad-tail conditional Oxford sentinel. It compares an
existing `10 nm / 4 s` CPU case, an exact CUDA replicate, `10 nm / 2 s`, and
the completed `5 nm / 2 s` CUDA case.

## Frozen-gate verdict

Eight of nine numerical checks pass. CPU/CUDA replay is exact; timestep depth
and CD pass; 10-to-5 nm depth, x/y symmetry, particle balance, and conservative
state remap pass. The one failure is the 10-to-5 nm CD gate: the maximum CD
change is `8.612 nm`, above the frozen `5 nm` tolerance. Therefore bottom CD
and the quantitative footing magnitude are **not numerically certified**.

The `5 nm` endpoint has depth `260.797 nm`, top CD `197.311 nm`, middle CD
`197.593 nm`, and the existing bottom metric `208.814 nm`. Its positive
bottom-minus-top value (`11.503 nm`) cannot be promoted to a quantitative
target prediction.

## Why this is not nonphysical material growth

The mechanism permits only finite, nonnegative TiO2 removal. Cr and fused
silica are pinned; growth and redeposition are disabled; negative removal
velocity is rejected. The remaining shape can therefore widen downward only
because the shadowed lower wall recedes more slowly than the upper wall.

The failed change is also localized. Every sampled width from the top through
`75%` of relief changes by at most about `2.51 nm` between the 10 and 5 nm
meshes. The frozen `5 nm` CD threshold is first exceeded at `85%`, and the
largest change occurs at `90%`. This sentinel etches only `260.8 nm` into a
`700 nm` film, leaving about `439.2 nm` of continuous TiO2 below the floor.
There is no independent pillar bottom at this time; the deepest samples lie
in the physical sidewall-to-unetched-film junction.

That post-result localization does not alter the preregistered failure. It
means the body profile is stable, the relief junction is physical, and the
exact flare magnitude remains grid-unresolved. It does not certify a final
Oxford pillar base or validate the omitted TiO2/Cr surface law.

```bash
python scripts/audit_zhu_npg80_profile_convergence_refinement.py --check
pytest -q tests/test_zhu_npg80_profile_convergence_refinement.py
```
