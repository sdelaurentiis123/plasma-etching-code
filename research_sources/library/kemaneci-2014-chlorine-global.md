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
- **Official implementation cross-check:**
  `https://doc.comsol.com/6.4/doc/com.comsol.help.models.plasma.chlorine_global_model/chlorine_global_model.html`
- **PDF SHA-256:**
  `fda9f3da209e31993b2c1405cad86a5838185d0e0dc297a1131e744243466c75`
- **Local extraction:**
  `research_sources/thesis_extracts/kemaneci_2014_chlorine_global.txt`
- **Status:** PRIMARY FULL TEXT + EQUATIONS 13--18 VISUALLY AUDITED AT 500 DPI
  + TABLE 4 AND FIGURE 10 VISUALLY AUDITED AT 600 DPI
- **Visual audit:** publisher-PDF page 6 was rendered at 500 dpi; pages 8 and
  15 were rerendered at 600 dpi for the coefficient-sign and level audits.
  The energy-balance page, chemistry table, and Figure-10 render SHA-256 values
  are, respectively,
  `30ddb3a4ba32625b9a2447d7ff5b68e04347370744f26631a9fafa25cf43f69b`,
  `fa796f2f85ffb8b98a50f87d06d18af7f0c703826720420f64d80b4b61e3363a`,
  and `b2492597fd938b5980aedfefd2a7baa63710a1320998b83259cdf0d181d0a898`.

## Claims table

| # | verified source claim | use and boundary |
|---|---|---|
| K1 | The model is volume averaged and Maxwellian, evolves particle balances plus electron energy density `p_e = 3/2 n_e T_e`, and supplies heavy-particle temperature externally. | A detailed source-reproduction target, not evidence that spatial transport, EEDF shape, gas heating, or RF absorption is closed. |
| K2 | The particularly studied chamber is a hard-anodized-Al cylinder of radius 0.275 m and length 0.10 m, with a four-turn 13.56 MHz planar coil and 50 sccm pure Cl2 feed. | Defines one equipment board. It is neither the Malyshev Lam Alliance geometry nor a generic Lam reactor. |
| K3 | Table 2 resolves to ten heavy species plus electrons when unqualified `Cl2` is treated as `Cl2(v=0)`. Table 4 states a fit domain of `T_e = 0.5--10 eV` and says most electron-rate fits were adapted from Thorsteinsson while selected rates were recomputed from cited cross sections. | Every landed row must remain a published-compilation/source-reproduction rate with its temperature domain; the table is not a uniform measurement set. |
| K4 | Equations 13--18 split absorbed power, chemical reaction exchange, charged-particle wall loss, and elastic electron loss. Equation 17 computes a reaction-energy scalar from the printed species internal energies. | Supplies the ledger topology and elastic-loss form. It does not justify replacing a collision-energy moment with an Arrhenius fit exponent. |
| K5 | Figure 10 prints Cl2 vibrational levels at 0.07, 0.14, and 0.21 eV; Cl fine/excited levels at 1.25, 1.35, and 10.17 eV; Cl2+ at 11.50 eV; Cl+ at 14.25 eV; and Cl- at -2.36 eV. | These values can reproduce the source's level convention. NIST independently supplies the predictive molecular thresholds and dissociation energy. |
| K6 | Figure 10 places ground-state Cl2 at zero and ground-state Cl at `1.25 eV` per atom, so its two-atom asymptote is `2.50 eV`. The primary evaluated `D0(Cl2)` board is `2.4793 eV`, while the official COMSOL reproduction separately assigns `ediss=4 eV`. Table 4's dissociation-rate fit contains an `8.84 eV` exponential parameter. | The Figure-10 convention does encode an approximate bond asymptote; the COMSOL event-energy input and the fit exponent are distinct, conflicting quantities. A fundamental ledger must use the evaluated dissociation threshold rather than silently selecting any of the three. |
| K7 | The paper imposes gas temperature over 300--1500 K rather than solving heavy-particle heating, and states that self-consistent gas heating requires a coupled heavy-particle energy balance. | No gas-temperature or depth prediction may be claimed from this reproduction rung. |
| K8 | The paper compares several published chlorine measurements and reports generally good agreement, while also noting that updated vibrational cross sections materially change `Cl2(v=1)`, electron temperature, and electron density. | Supplies preregistration candidates and a mechanism-sensitivity warning; it is not one clean held-out board because source conditions and wall closures vary across comparisons. |
| K9 | Table 4 labels charge exchange as reactions `(28)--(32)` while printing `Cl2(v=0--3)`, which expands to four channels, not five. The official COMSOL reproduction implements four channels. Table 4 stars eight excitation rows (10--15 and 17--18), and the raw COMSOL 6.4 model contains 38 forward rows (including two elastic channels) plus all eight reverses = 46 volume features. | The printed range is quarantined as an off-by-one source defect. The native forward replay expands exactly four charge-exchange channels and records the missing source label rather than inventing a fifth species or event. A 44-count is valid only for the non-elastic forward+reverse subset. |
| K10 | A 600-dpi visual audit confirms reaction 17 is `4.55e-14 Te^-0.46 exp(-2.01/Te - 0.001/Te^2)`. The first native replay and its formula-shaped test both carried `Te^+0.46`; the official COMSOL model independently exposes the negative exponent. | Corrected before any Kemaneci solve. The regression now uses an independently precomputed 50-digit value at `Te=2.3 eV`, rather than repeating the implementation formula. |
| K11 | The official COMSOL 6.4 model uses the Table-4 reaction-20 fit with `exp(-13.29/Te)`, whereas the primary paper prints `exp(-13.19/Te)`. | The native primary-paper replay retains `-13.19`; the COMSOL divergence is quarantined and must be selected explicitly in any implementation-reproduction mode. |

## Use decision

Kemaneci is the strongest available full-text blueprint for a detailed native
chlorine source-reproduction tier. It provides the expanded species/reaction
topology, published fit domain, wall and pulse structure, and the separation of
chemical, wall, and elastic energy terms.

It does **not** close the predictive electron-energy ledger. In particular,
its printed level convention gives a `2 x 1.25 = 2.50 eV` atomic asymptote,
while the COMSOL reproduction separately uses `ediss=4 eV`; neither replaces
the measured `D0(Cl2)=2.4793 eV` stored in
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

The first source-reproduction rung is now executable as
`build_kemaneci_2014_forward_chlorine_network()`: 36 printed non-elastic
forward reactions, ten heavy states plus electrons, exact atom/charge closure,
and strict enforcement of the `0.5--10 eV` fit domain. The two elastic channels
and eight detailed-balance reverses remain explicit next gates, so this
forward network cannot yet be called the paper's complete 44-reaction
non-elastic model (46 features when the two elastic rows are included).
