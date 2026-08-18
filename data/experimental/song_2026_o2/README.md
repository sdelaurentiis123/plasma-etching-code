# Song et al. 2026 O2 electron collisions

This directory records the provenance and ingestion contract for the 2026
AIP evaluation of electron collisions with molecular oxygen. The supplemental
workbook is licensed CC BY-NC 4.0 and is intentionally **not** copied into the
repository. `source_manifest.json` locks the exact workbook by SHA-256 and
records its official Figshare download.

`petch.reactor_global.o2_electron_collisions` reads a caller-supplied matching
workbook as OOXML, without changing it, and constructs a deterministic
Boltzmann deck. The derived deck includes momentum transfer, vibration, five
electronic states, neutral dissociation, total positive ionization, and
dissociative attachment. Every extrapolation beyond the tabulated support is
named and hashed.

For the Zhu/NPG80 target feed, the committed audit is
`results/curated/zhu_npg80_feed_electron_kinetics_v1/audit.json`. That audit
uses no SEM or etch-depth outcome. It establishes only the local
species-resolved electron-kinetic rung over a declared bulk-field scan; the
global particle/power balance remains responsible for selecting a reactor
state.
