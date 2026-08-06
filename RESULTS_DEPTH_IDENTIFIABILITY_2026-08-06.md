# Absolute depth: the ceiling proof is retracted; the boundary is underidentified

## Bottom line

The aggregate-boundary 60 s simulation remains **346.833 nm**, versus the
measured **825 nm**. That numerical MISS is real. The former proof that 825 nm
was physically unreachable from Krüger's published inputs is not.

The proof failed because it compared petch's species-agnostic `_complex_yield`
formula to one CF3+ beam marker and treated the marker's rounded
`1.5 SiO2/ion` value as a universal ceiling. Direct end-to-end replay showed
that the default mechanism actually produced the same `0.380584 SiO2/ion` for
F+, CF+, CF2+, and CF3+ at 1000 eV: energetic ion identity was discarded.
Karahashi's full Figure 4 instead contains a strong species ladder and reaches
`1.8736` for CF3+ at 1500 eV.

An independent molecule/ion beam experiment then removes the premise
completely. Takada, Toyoda, and Sugai measured stable C5F8 molecule + Ar+
co-incidence yields of `~1.2` at 400 eV and `2.4–2.5` at 900 eV, both at
molecule/ion ratio 1. The range reflects a real source discrepancy: the open
TMRSJ paper says `2.5`, while the related *Journal of Applied Physics*
publisher abstract says `2.4`. That is 95.2–99.2% of the `2.521` run-average
wafer-ion normalization demanded by 825 nm. It is not a C4F6 law, but it proves
that a pure-CF3+ value is not a ceiling on fluorocarbon-plasma surface removal.

The corrected verdict is therefore:

- the existing aggregate-boundary run misses absolute depth by 58%;
- measured surface mechanisms omitted from that run are large enough to make
  the old impossibility claim false;
- Krüger publishes neither the positive-ion species composition nor the stable
  C4F6 flux needed to evaluate those mechanisms;
- exact depth is **underidentified**, not proven impossible and not yet
  predicted.

The deterministic arithmetic is replayed by
`scripts/audit_krueger_depth_identifiability.py` and frozen in
`results/curated/depth_identifiability/audit.json`.

## 1. What the Karahashi pixels actually say

The source PDF was rendered at 600 dpi. PIL/NumPy axis localization,
marker-center transcription, error-cap retention, a checksum manifest, and a
full-resolution visual overlay are replayed by
`scripts/digitize_karahashi_2007_fig4.py`.

At 1000 eV and normal incidence:

| ion | digitized SiO2/ion |
|---|---:|
| F+ | 0.3232 |
| CF+ | 0.6751 |
| CF2+ | 1.1957 |
| CF3+ | 1.4703 |

The CF3+ series rises to `1.8736` at 1500 eV and remains `1.7549` at
2000 eV. Karahashi's statement that yield “gradually saturated” above 1000 eV
does not turn the 1000 eV rounded value into a hard maximum, and the experiment
contains no stable-molecule co-incidence.

The code now provides an opt-in
`Karahashi2007ReactiveIonYieldTable`. It linearly interpolates only within each
species' digitized positive-yield support, only at normal incidence, and
refuses unknown species, angles, and energy extrapolation. The tabulated
points are fitted/reproduction evidence, not an independent validation of the
closure. The aggregate Krüger run cannot enable it because that paper labels
the entire energetic population only as `Ions`.

## 2. The stable-parent channel exists in direct measurement

Takada's radical-free beam experiment used mass-selected Ar+ at normal
incidence and a C5F8 molecular beam at 45 degrees. Figure 3's SiO2 series,
digitized and visually audited at 600 dpi, is:

| C5F8/Ar+ ratio | 400 eV yield, SiO2/Ar+ |
|---:|---:|
| 0.25 | 0.6697 |
| 0.50 | 1.0162 |
| 1.00 | 1.1969 |
| 2.50 | 1.0591 |
| 10.0 | 0.7945 |

The source text independently reports `0.67` at ratio 0.25 and `~1.2` at ratio
1. At 900 eV and ratio 1, the open TMRSJ text reports `2.5` while the related
JAP publisher abstract reports `2.4`; both are retained. The non-monotone
400 eV curve matters as much as its magnitude: more fluorocarbon is not a free
speed multiplier. Excess supply builds a fluorocarbon layer and pushes the
system back toward deposition/etch stop.

This source is in the internal literature library as
`takada-2005-tmrsj`, with the archived PDF, full text, verbatim claims, reverse
index, digitized table, raw pixels, and claim boundary.

It is deliberately **not** imported as a C4F6 mechanism. Takada identifies the
C5F8 ring/double-bond adsorption route as potentially unique. Copying that
curve to C4F6 would replace one false certainty with another.

## 3. Correct arithmetic: run average is not final-floor yield

With SiO2 formula density `2.2e28 m^-3`, 60 s, and Krüger's published aggregate
wafer ion flux `1.2e20 m^-2 s^-1`:

| quantity | value | interpretation |
|---|---:|---|
| 825 nm target | 2.5208 SiO2 / wafer ion | depth-integrated lower-bound normalization |
| 346.833 nm simulation | 1.0598 SiO2 / wafer ion | same normalization |
| difference | 1.4611 SiO2 / wafer ion | unresolved effective removal |
| target / simulation | 2.3787 | numerical depth-rate gap |
| target / 0.70 final ion delivery | 3.6012 SiO2 / delivered floor ion | counterfactual final-geometry diagnostic only |

The old documents called `2.52` a floor yield and then also divided by 0.70 to
derive a `2.40x` wafer-flux bound. Those are different normalizations. A real
feature spends most of its time shallower than its final geometry; applying
the final delivery to the full 60 s is not a history integral.

Takada's reported `2.4–2.5` range is 95.21–99.17% of the target's
wafer-ion-normalized lower bound. That comparison establishes plausibility,
not closure: it is C5F8 at 900 eV, while Krüger used C4F6 with a broad
aggregate IEAD and an evolving feature.

## 4. How much omitted parent flux could exist?

Krüger reports 10 mTorr and feed flows C4F6/Ar/O2 =
140/100/105 sccm, but Table I omits stable C4F6 at the wafer. As a scale
reference, assigning the inlet flow fraction to the total pressure and
assuming an undissociated 300 K gas gives:

- C4F6 partial-pressure reference: `0.5410 Pa`;
- one-sided molecular-impingement reference: `6.465e21 m^-2 s^-1`;
- parent/ion ratio reference at the wafer: `53.88`.

This is not a reactor prediction or a strict upper bound: species-dependent
residence time and pumping can move the chamber fraction away from the inlet
flow fraction, even before plasma chemistry and surface coupling. Under the
existing final-geometry diagnostic deliveries (neutral 0.10, ion 0.70),
`3.25%` of this declared reference would produce floor parent/ion ratio 0.25,
and `13.0%` would produce ratio 1.

Those percentages must not be called measured survival fractions. They are
only a sensitivity scale. They show why omitting the parent boundary prevents
an impossibility proof: a few-percent availability relative to a plausible
feed scale is already an order-one surface variable in a measured analogous
system.

The reactor literature reinforces the missing-boundary diagnosis with direct
C4F6 measurements. Kim et al. put a mass/energy analyzer at the powered
electrode of a C4F6/Ar CCP. Their neutral/radical spectrum contains a distinct
C4F6 parent signal, and their positive-ion ordering is CF+, C3F3+, CF3+, CF2+,
then C3F5+. The energy distributions also depend on ion mass. This is not
Krüger's boundary—the experiment used C4F6/Ar at 20 mTorr and 300 W, without
O2, and reports detector count rates rather than calibrated absolute fluxes.
It does directly confirm that both missing boundary classes exist in a C4F6
plasma and that one aggregate IEAD is not species resolved.

Huang's similar fluorocarbon CCP simulation reports only about 24% feedstock
dissociation in its base case and finds large C2F4+ and C3F5+ ion populations.
Huard's lineage explicitly treats those large energetic ions as fragmenting
on impact and coupling energy delivery to fluorocarbon incorporation. The
current aggregate Krüger path implements neither a measured large-ion mixture
nor stable-parent co-incidence.

## 5. What changed in petch

1. Karahashi Figure 4 now has a checksum-bound, species-resolved data table,
   digitization manifest, PIL replay, and visual QA overlay.
2. An opt-in species-resolved reactive-ion beam closure reproduces the full
   F+/CF+/CF2+/CF3+ ladder end to end and refuses unsupported use.
3. The default Ar/aggregate path is bitwise unchanged.
4. The false 4.7% “independent validation” gates were removed. The remaining
   supply-bound gates apply only to the incumbent neutral-assisted complex
   channel and explicitly do not bound missing channels.
5. Takada Figure 3 now has the same archived-source, full-text, library,
   digitization, checksum, loader, and visual-audit chain.
6. The `ion_flux_normalization` option remains available only as an explicit
   fit to the 825 nm feature target. It is not a blanket calibration: the
   source publishes no blanket etch datum.

## 6. Atomic-accuracy path to an actual depth prediction

Atomic-level accuracy is a convergence target, not a label earned by adding
more adjustable rates. For this reactor the minimum clean closure is:

1. species-resolved positive-ion flux and IEAD at the wafer, including the
   large CxFy+ population;
2. stable C4F6 molecular flux at the wafer;
3. direct C4F6 + representative-ion beam data across energy, angle, and
   molecule/ion ratio, including the deposition transition;
4. species-aware impact fragmentation and C/F incorporation with conservative
   product routing;
5. a reactor model that predicts items 1–2 from knobs and is calibrated on
   diagnostics other than the feature depth;
6. full evolving-feature prediction held out from that calibration.

Sean's 0-D reactor model is structurally the right way to generate a measured
boundary from machine knobs. It resolves the dependency, but it does not close
this C4F6 case until its chemistry is extended and graded against
species-resolved diagnostics.

No ion-flux scale, C5F8 analog law, or Karahashi species mixture has been fitted
to 825 nm in this correction. The science now stops exactly where the evidence
stops.
