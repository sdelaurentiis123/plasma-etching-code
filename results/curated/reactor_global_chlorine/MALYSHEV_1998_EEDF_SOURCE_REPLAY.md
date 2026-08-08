# Malyshev 1998 molecular-only EEPF failure receipt

## Verdict

The molecular-only collision deck is **retired under the corrected Malyshev
Eq. 11 pressure closure**. It creates a dissociated chlorine population but
contains no atomic-Cl electron collision target, so the coupled particle and
charge system does not close. This failure is recorded rather than preserving
the stale source-replay numbers. The atomic-Cl replay is the minimum operative
deck.

| absorbed fraction | source W | failure | interpretation |
|---:|---:|---|---|
| 0.50 | 300 | RuntimeError | molecular-only electron collision deck cannot close a dissociated chlorine state; atomic-Cl electron collisions are mandatory |
| 0.50 | 500 | RuntimeError | molecular-only electron collision deck cannot close a dissociated chlorine state; atomic-Cl electron collisions are mandatory |
| 0.30 | 300 | RuntimeError | molecular-only electron collision deck cannot close a dissociated chlorine state; atomic-Cl electron collisions are mandatory |
| 0.30 | 500 | RuntimeError | molecular-only electron collision deck cannot close a dissociated chlorine state; atomic-Cl electron collisions are mandatory |
| 0.70 | 300 | RuntimeError | molecular-only electron collision deck cannot close a dissociated chlorine state; atomic-Cl electron collisions are mandatory |
| 0.70 | 500 | RuntimeError | molecular-only electron collision deck cannot close a dissociated chlorine state; atomic-Cl electron collisions are mandatory |

No reactor state, wafer flux, or feature-depth claim is supported by this
negative board. No feature observable selected any coefficient.
