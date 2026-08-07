# Krüger eight-ray QMC receiver-area singularity

Status: **diagnosed by an exact frozen-checkpoint paired control; QMC excluded
from moving-profile authority.**

At accepted step 85 (`t=10.249480050516949 s`) of the 10 nm Guo/Krüger
trajectory, the general-3D diffuse-neutral operator produced a maximum face
speed of `3.5893e-6 m/s`. The exact same geometry, boundary, energetic
transport, surface state, chemistry, and seed was reevaluated at zero
duration with diffuse neutral exchange disabled. Its maximum speed was
`2.3985e-8 m/s`.

The anomalous receiver is mask face 652:

- area `1.05866e-21 m2`, only `2.11277e-5` of the active-face median;
- no incident energetic-ion flux;
- an eight-ray categorical source tally assigns finite diffuse-neutral weight
  to it and then divides by its vanishing receiver area;
- the resulting C3F4 flux is `6.16805e24 m-2 s-1`;
- five narrow-band nodes inherit the face value, making it profile-active.

With diffuse exchange off, every target-face neutral flux and both etch and
growth velocities are exactly zero. Material ledgers remain exact in both
arms. The observed spike is therefore an estimator/mesh-receiver
singularity, not evidence for a new yield, flux scale, chemistry branch, or
surface-smoothing rule.

For Krüger's translationally invariant trench, the candidate authority is the
deterministic extruded 2-D exchange operator. General 3-D may re-enter only
after replicated physical-patch-local exchange/velocity convergence or an
error-controlled hierarchical deterministic estimator. Integrated flux
convergence alone is insufficient.

Reproduce from the pinned local checkpoint:

```bash
python scripts/audit_krueger_qmc_receiver_singularity.py \
  --input-directory /private/tmp/krueger_guo_transient_dt125_dx10 \
  --output /private/tmp/krueger_qmc_receiver_singularity_audit.json
```

The compact committed receipt records the checkpoint/audit/config hashes and
the complete paired-control result.
