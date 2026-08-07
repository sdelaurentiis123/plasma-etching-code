# DFT-trained NNP/ZBL transfer to Karahashi beams

This board tests the released An et al. atomistic SiO2 outputs on all nine
exact-overlap Karahashi CF+, CF2+, and CF3+ beam points. Target, ion identity,
energy, and normal incidence must match exactly; the audit performs no
interpolation, extrapolation, or parameter fit.

Run:

```bash
python scripts/audit_an_2026_nnp_karahashi_transfer.py --check
```

The result is a real advance, not yet a universal surface closure. On identical
nine-point support the DFT-trained NNP/ZBL calculation cuts the error sharply
relative to the frozen Guo/Kwon chemistry and gives the correct etch sign for
every point. At 1000 eV all three species are within 10%, and the post-hoc
`>=750 eV` subset has about 10.6% MAPE.

The full overlap is not called validated: two low-energy points miss by
67–72%, and the authors' comparison was not blind. Their own transient test
locates the defect in product escape and film evolution: deleting stable
products from any depth delays the etch-to-deposition transition and permits
more than 10 nm of removal where the cited experiment gives about 2 nm.

The code-level audit is pinned too. At author commit `4bcd0350`,
`SM_codes/PlasmaEtchSimulator/calc/byproduct.py` deletes every disconnected
cluster whose composition matches the configured byproduct list without a
formation-depth or residence-time test. The optional depth-gated path fixes
the region to 20 Å and labels that value arbitrary in source. The committed
audit stores the source and released-transient checksums; no unlicensed author
implementation or weights were copied into petch.

The physics consequence is architectural. Use atomistic calculations to
derive the prompt, species/energy-resolved event kernel; couple that kernel to
a finite atom-balanced mixed layer, fluorocarbon film, and explicit
depth-/time-dependent product escape. Do not turn these outputs into a scalar
Krüger depth knob. The released data stop at 1000 eV, include no radical
coflux, and do not identify Krüger's unpublished ion mixture.
