# kemaneci-2014-chlorine-global

**Detailed continuous/pulsed chlorine global-model reproduction source**

- **Citation:** E. Kemaneci, E. Carbone, J.-P. Booth, W. Graef, J. van Dijk,
  and G. Kroesen, “Global (volume-averaged) model of inductively coupled
  chlorine plasma: Influence of Cl wall recombination and external heating on
  continuous and pulse-modulated plasmas,” *Plasma Sources Science and
  Technology* **23**, 045002 (2014).
- **DOI:** `10.1088/0963-0252/23/4/045002`
- **Official university full text:**
  `https://pure.tue.nl/ws/files/3931833/844662334177360.pdf`
- **PDF SHA-256:**
  `fda9f3da209e31993b2c1405cad86a5838185d0e0dc297a1131e744243466c75`
- **Local extraction:**
  `research_sources/thesis_extracts/kemaneci_2014_chlorine_global.txt`
- **Status:** PRIMARY FULL TEXT + EQUATIONS 13--18, TABLE 4, AND FIGURE 10
  VISUALLY AUDITED AT 500 DPI
- **Visual audit:** publisher-PDF pages 6, 8, and 15 were rendered at 500 dpi.
  The energy-balance page, chemistry table, and Figure-10 render SHA-256 values
  are, respectively,
  `30ddb3a4ba32625b9a2447d7ff5b68e04347370744f26631a9fafa25cf43f69b`,
  `65bb47a5f7da3b2642165fa47d48542f8660719dd188ef8f59b28380e371f78d`,
  and `3416ab6b360bbe16a97be93f55152077ecd50924c42fe7c243db80a446ff9a13`.

## Claims table

| # | verified source claim | use and boundary |
|---|---|---|
| K1 | The model is volume averaged and Maxwellian, evolves particle balances plus electron energy density `p_e = 3/2 n_e T_e`, and supplies heavy-particle temperature externally. | A detailed source-reproduction target, not evidence that spatial transport, EEDF shape, gas heating, or RF absorption is closed. |
| K2 | The particularly studied chamber is a hard-anodized-Al cylinder of radius 0.275 m and length 0.10 m, with a four-turn 13.56 MHz planar coil and 50 sccm pure Cl2 feed. | Defines one equipment board. It is neither the Malyshev Lam Alliance geometry nor a generic Lam reactor. |
| K3 | Table 4 contains a ten-species chlorine mechanism, states a fit domain of `T_e = 0.5--10 eV`, and says most electron-rate fits were adapted from Thorsteinsson while selected rates were recomputed from cited cross sections. | Every landed row must remain a published-compilation/source-reproduction rate with its temperature domain; the table is not a uniform measurement set. |
| K4 | Equations 13--18 split absorbed power, chemical reaction exchange, charged-particle wall loss, and elastic electron loss. Equation 17 computes a reaction-energy scalar from the printed species internal energies. | Supplies the ledger topology and elastic-loss form. It does not justify replacing a collision-energy moment with an Arrhenius fit exponent. |
| K5 | Figure 10 prints Cl2 vibrational levels at 0.07, 0.14, and 0.21 eV; Cl fine/excited levels at 1.25, 1.35, and 10.17 eV; Cl2+ at 11.50 eV; Cl+ at 14.25 eV; and Cl- at -2.36 eV. | These values can reproduce the source's level convention. NIST independently supplies the predictive molecular thresholds and dissociation energy. |
| K6 | Figure 10 places both ground-state Cl2 and ground-state Cl at zero, while Table 4's ground-state dissociation fit has an `8.84 eV` exponential parameter. | **Inference from the audited source:** the Figure-10 level convention cannot encode the Cl2 bond-dissociation energy, and the fit exponent is not a thermochemical event energy. A fundamental ledger must source those quantities separately. |
| K7 | The paper imposes gas temperature over 300--1500 K rather than solving heavy-particle heating, and states that self-consistent gas heating requires a coupled heavy-particle energy balance. | No gas-temperature or depth prediction may be claimed from this reproduction rung. |
| K8 | The paper compares several published chlorine measurements and reports generally good agreement, while also noting that updated vibrational cross sections materially change `Cl2(v=1)`, electron temperature, and electron density. | Supplies preregistration candidates and a mechanism-sensitivity warning; it is not one clean held-out board because source conditions and wall closures vary across comparisons. |

## Use decision

Kemaneci is the strongest available full-text blueprint for a detailed native
chlorine source-reproduction tier. It provides the expanded species/reaction
topology, published fit domain, wall and pulse structure, and the separation of
chemical, wall, and elastic energy terms.

It does **not** close the predictive electron-energy ledger. In particular,
its printed level convention sets both molecular and atomic ground states to
zero, so Equation 17 cannot recover the measured `D0(Cl2)=2.4793 eV` stored in
`christophorou-olthoff-1999-cl2`. Dissociative attachment also removes an
electron at a collision-selected kinetic energy; that loss requires the
Boltzmann energy moment `<sigma v E>`, not the Table-4 rate-fit exponent or a
single Figure-10 level difference.

The native implementation must therefore keep two tiers side by side:

1. a clearly labeled Kemaneci source reproduction for code verification and
   comparison with the paper's boards; and
2. an evaluated physical tier whose reaction thresholds and energy moments
   come from measured/evaluated collision data and which fails closed wherever
   those data are incomplete.

Neither tier may tune a reactor constant to feature depth.
