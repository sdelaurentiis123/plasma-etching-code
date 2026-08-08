# Native chlorine elastic-energy gate - 2026-08-07

## Result

The electron-energy infrastructure now separates a collision-exact
stationary-target momentum-transfer moment from the Kemaneci Equation-18
source approximation. No chlorine elastic cross section has been promoted to
evaluated physical authority.

For electron mass `me`, target mass `M`, incident energy `E`, and
center-of-mass deflection angle `chi`, two-body elastic kinematics gives the
electron-to-target energy transfer

`deltaE = [2 me M/(me+M)^2] E (1-cos(chi))`.

After angular integration, `sigma_m = integral (1-cos(chi)) dsigma`, so the
stationary-target Maxwellian loss coefficient is

`[2 me M/(me+M)^2] <sigma_m v E>`.

It must use the energy moment from the same momentum-transfer table. Kemaneci
Equation 18 instead uses

`3 Te (me/M) <sigma_m v>`

under an explicit `Te >> Th` assumption. Both are executable, named modes. For
a constant cross section, collision weighting gives
`<sigma v E>/<sigma v> = 2 Te`; in the heavy-target limit the exact form is
therefore `4 Te me/M <sigma v>`, while the Kemaneci form is `3 Te me/M`, about
25 percent lower.

The physical mode is exact only for a stationary heavy target. A future
two-temperature Boltzmann operator must account for target thermal motion when
`Th` is not negligible; the current name and provenance forbid that extension
from happening silently.

## Recovered official implementation assets

The COMSOL 6.2, 6.3, and 6.4 chlorine model downloads contain byte-identical
momentum-transfer tables. Their identities are pinned in
`research_sources/digitized/comsol-6.4-chlorine-global-assets-manifest.md`; raw
bytes are not redistributed.

| target | numeric records | energy support | SHA-256 prefix |
|---|---:|---:|---|
| Cl2 | 49 | `0.01983--29.17 eV` | `6bf31f960206...` |
| Cl | 28 | `0--25 eV` | `6c1391eda01d...` |

The Cl2 file's last record lacks a newline, so `wc -l` reports 48 despite 49
numeric records. Parsers and audits must count records, not newline bytes.

## Maxwellian-support audit

The table ceiling matters more for the exact energy moment than for the event
rate. The fractions below are the constant-cross-section Maxwellian kernel
above the last tabulated energy; they diagnose missing support without
inventing a cross-section extrapolation.

| target | Te (eV) | missing rate kernel | missing energy kernel |
|---|---:|---:|---:|
| Cl2 | 3 | 0.0642% | 0.3472% |
| Cl2 | 5 | 1.9999% | 6.9799% |
| Cl2 | 10 | 21.1893% | 44.2040% |
| Cl | 3 | 0.2243% | 1.0590% |
| Cl | 5 | 4.0428% | 12.4652% |
| Cl | 10 | 28.7297% | 54.3813% |

The raw tables therefore do not close the full Kemaneci `0.5--10 eV` fit
domain without a sourced high-energy extension. The runtime's strict tail
gates correctly refuse unsupported temperatures.

## Evidence verdict

- The NIST chlorine review supplies evaluated **total elastic** Cl2 cross
  sections but identifies molecular momentum-transfer measurements as missing.
  Total elastic cannot replace momentum transfer without an angular-scattering
  assumption.
- Gregorio--Pitchford's 2012 recommended Cl2 momentum-transfer row is
  Rescigno's close-coupling calculation with a low-energy extrapolation. The
  paper obtains poor agreement with the old drift and characteristic-energy
  data, declines to retune momentum transfer because those measurements are
  unreliable, and explicitly calls for new measurements.
- Gonzalez-Magana--de Urquijo supplied those newer direct pulsed-Townsend
  measurements in 2018 over `3--450 Td`. Kawaguchi et al. published an updated
  cross-section set in 2020 and validated transport against measured data. The
  full papers/curves are not locally archived, so neither source can yet land a
  number or a pass.
- Griffin et al. and Wang et al. provide calculated atomic-Cl elastic and
  momentum-transfer data. Griffin validates method quality indirectly against
  Ar; Wang reports strong low-energy model dependence near the Ramsauer
  minimum and no resolving experiment.
- The COMSOL tables can support an exact implementation replay after their
  node mapping is verified. The Cl2 table is not measured evidence, and neither
  table can yet support an evaluated physical electron-power claim.

## Next gates

1. Archive and pixel-audit the full Kawaguchi 2020 and
   Gonzalez-Magana--de Urquijo 2018 articles; grade the Gregorio and Kawaguchi
   sets against the direct swarm markers without tuning.
2. Audit the full Griffin 1995 paper against the COMSOL atomic-Cl table nodes
   and provenance.
3. Recover author-supplied Wang fine-structure arrays rather than digitizing a
   lossy plot if possible.
4. Extend the high-energy support from primary/evaluated data or keep each
   temperature outside support fail-closed.
5. Only then add the two elastic rows to the 44-row nonelastic COMSOL replay
   and compare the source approximation with the exact stationary-target
   energy moment side by side.
