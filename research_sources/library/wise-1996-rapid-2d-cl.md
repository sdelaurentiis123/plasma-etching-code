# wise-1996-rapid-2d-cl

**Rapid quasineutral two-dimensional chlorine ICP reference**

- **Citation:** R. S. Wise, D. P. Lymberopoulos, and D. J. Economou,
  “Rapid two-dimensional self-consistent simulation of inductively coupled
  plasma and comparison with experimental data,” *Applied Physics Letters*
  **68**, 2499--2501 (1996).
- **DOI:** `10.1063/1.116171`
- **Official author/university PDF:**
  `https://www.chee.uh.edu/sites/chbe/files/faculty/economou/apl_mpres_1996.pdf`
- **PDF SHA-256:**
  `25bfedd2ebf45bcb62b237ea0a5e27d1e495850673d73f3ded4c8b8321999e65`
- **Local verified extraction:**
  `research_sources/thesis_extracts/wise_1996_rapid_2d_cl_verified_excerpt.txt`
- **Status:** PRIMARY FULL TEXT READ; EQUATIONS/BOUNDARIES TEXT-VERIFIED;
  FIGURE 3 DIRECT MARKERS PIL-DIGITIZED FROM A SHA-PINNED 300-DPI RENDER

## Claims table

| # | verified source claim | use and boundary |
|---|---|---|
| W1 | The solver is modular: electromagnetic deposition feeds electron energy/rates, which feed separate charged- and neutral-species transport modules; the modules iterate to self-consistency. | Primary architecture for the deterministic reactor tiers and block iteration. It does not validate the present chemistry deck or a Lam geometry. |
| W2 | The bulk is quasineutral; positive- and negative-ion particle balances are solved and electron density is recovered from electroneutrality. | Direct basis for the charged spatial tier. A scalar electropositive reaction-diffusion equation is not an admissible chlorine-ion substitute. |
| W3 | Electrons are taken in Boltzmann equilibrium and the bulk electric field is recovered from their pressure gradient. | Basis for the current quasineutral potential closure. Spatial electron temperature and nonlocal EEDF coupling remain a higher tier. |
| W4 | At the wall the positive-ion loss is local density times local Bohm velocity, negative-ion density is set to zero, and the sheath is treated outside the quasineutral bulk. | Basis for the species-resolved wall boundary. This does not supply biased-wafer IEADs; those remain in the separate sheath provider. |
| W5 | The source board is a GEC ICP at 20 mTorr, 180 W, 20 sccm pure chlorine, 13.56 MHz, with no wafer/bias/etching; the deposited power is toroidal and has an approximately 1 cm axial skin depth. | Supplies source-shape sensitivities and a future independent spatial board, not a Mahorowala or Lam boundary calibration. |
| W6 | Predicted line-integrated electron density, central negative-ion density, and radial electron-density/temperature/potential profiles were compared with GEC-cell measurements without adjusted rates; negative-ion agreement was within the stated factor-two measurement uncertainty, while predicted electron temperatures were about 50% below probe values. | A mixed quantitative validation target. The electron-temperature miss must remain visible; the paper cannot be cited as universal validation of reactor depth. |
| W7 | Figure 3 contains 21 directly measured Langmuir-probe markers: seven radii each for electron density, electron temperature, and plasma potential at the 180 W / 20 mTorr / 20 sccm base condition. | The checksum-pinned board is `data/experimental/wise_1996_gec_icp/figure3_radial_measurements.csv`; it can grade spatial closure without selecting any feature-depth parameter. The source prints no marker error bars, so the board does not support an uncertainty-weighted pass. |

## Use decision

This is the primary architecture source for the new spatial charged tier.  The
implemented rung retains quasineutrality, Boltzmann electrons, separate
positive/negative ion transport, exponential-fit drift fluxes, a Bohm
positive-ion wall, and a zero-density negative-ion wall.  It conditions
distributed sources on an already solved 0-D inventory, so it is presently a
wafer-partition prediction rather than a full 2-D reactor-state prediction.

The GEC hardware geometry and measured radial profiles are now landed from
the primary Wise and Miller records.  Promoting the present conditional lift
to a reactor-state prediction still requires a measured or independently
solved coil field, spatial electron energy and electron-impact chemistry, and
local nonlinear mutual neutralization.  No feature depth may select any of
those quantities.
