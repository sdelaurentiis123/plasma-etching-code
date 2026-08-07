# mahoney-1994-planar-icp

**Mahoney et al., spatial electron state in a planar argon ICP**

- **Citation:** L. J. Mahoney, A. E. Wendt, E. Barrios, C. J. Richards, and
  J. L. Shohet, “Electron-density and energy distributions in a planar
  inductively coupled discharge,” *Journal of Applied Physics* **76**,
  2041–2047 (1994).
- **DOI:** `10.1063/1.357672`
- **Primary repository:** University of Wisconsin MINDS, record `1793/9564`.
- **PDF SHA-256:**
  `acd59d5def6373a81f1ec73248608b7c67a8769b010f2287fa05b04fd8cc61b7`
- **Project extract:**
  `research_sources/thesis_extracts/mahoney_1994_planar_icp.txt`
  (SHA-256
  `1636904c776fd88dbc98d0b15cb30ecd36788cf0fa925b4059e2b33d4e216d7e`).
- **Visual audit:** Table I on PDF page 4 and Figure 11 on PDF page 6 at
  200 dpi; checksums are in the digitization manifest.
- **Status:** PRIMARY FULL TEXT READ + TABLE I/FIGURE 11 VISUALLY AUDITED
- **Topic:** independent pure-argon ICP plasma-state validation and RF-power
  boundary failure

## Claims table

| # | verified source claim | use and boundary |
|---|---|---|
| Q1 | The source uses a four-turn planar induction coil at 13.6 MHz, a quartz window, a Faraday shield, and a grounded metal liner with 22.8 cm inner diameter and 13.7 cm length. | Defines a condition-specific planar ICP geometry. It is not a generic Lam chamber and not Krüger’s CCP. |
| Q2 | Electron density and EEDFs were measured with RF-filtered Langmuir probes; the retained peak observable is on axis at `z = 5.0 cm`. | The global-model center density is compared to a spatial peak, not a volume average. |
| Q3 | Table I prints five 100 W argon rows spanning 10–100 mTorr, including two independent 20 mTorr pumping configurations. | All rows are retained in the frozen independent board; the duplicate is not averaged away. |
| Q4 | Reported net RF power is incident minus reflected power, with reflected power below 5% after matching. The paper explicitly warns of additional induction-coil and matching-network dissipation. | `100 W` is an upper bound on absorbed plasma power, not a measured absorbed-power input. |
| Q5 | The paper reports that its electron-density result can be two to five times below ion-density determinations in similar electropositive planar ICPs. | The independent board froze model/measured center-density ratio `[1,5]`; it was not widened after the run. |
| Q6 | The source discusses saturation of density with increasing power and decreasing coupling efficiency as plasma resistance changes. | A constant generator-to-plasma efficiency is not source-authorized over arbitrary power or state. |
| Q7 | The unchanged five-reaction argon closure passed conservation, pressure trends, and normalized density-shape gates, but failed the frozen absolute-density interval when `100 W` was treated as absorbed. Across the four temperature/wall-energy corners it predicted 4.05–19.09 times the measured electron density. | Independent board verdict is **FAIL**. It isolates the missing RF-to-absorbed-power and gas-temperature boundaries; it is not repaired by tuning reaction rates. |
| Q8 | A preregistered target inversion found one constant transfer-fraction interval satisfying all five density rows inside each sensitivity corner: 21.9–26.2%, 26.3–30.4%, 41.4–50.4%, and 49.7–59.0%. None overlaps Hopwood's separate 70–90% apparatus context. | Diagnostic only. The inferred interval changes more than twofold with unresolved gas temperature/wall energy and is prohibited as a production input. |

## Use decision

Mahoney is the first independent condition-specific plasma-state board for the
new reactor kernel. It is deliberately preserved as a failed gate. The result
does not falsify conservation or the measured pressure response; it proves
that forward-minus-reflected generator power cannot be silently passed to the
global power balance as plasma absorption.

The source supplies no condition-specific calorimetry or hardware-loss map.
Any later inference of an absorbed fraction from its measured density is a
diagnostic inversion, not a predictive validation.
