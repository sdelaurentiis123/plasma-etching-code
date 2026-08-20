# Capillary-force pattern collapse: published criterion, sourced constants, and the square-pillar adaptation

Date: 2026-08-20
Question: given a predicted TiO2 pillar geometry (square pillars, 80–320 nm
wide, up to 700 nm tall, 400 nm pitch, on fused silica), what does the published
physics say about whether the pillars survive the rinse/dry step, and where is
the collapse threshold as a function of width?

**Doctrine note.** Every formula and number below is quoted from a source that
was fetched and archived this session, or is explicitly marked
`[unquoted — verify on next use]`. Anything we derived ourselves is in a section
headed **OUR ADAPTATION** and is never presented as published physics.
Retrieval routes are given so the next agent does not re-hunt.

---

## 0. Bottom line (read this first)

1. The classic criterion is **Tanaka, Morigami & Atoda, JJAP 32, 6059 (1993)**,
   derived for **lines**. We could not read the original (paywalled); we carry
   its equations verbatim from **Mack (2006)**, and independently cross-check
   the transcription against **Chandra's UPenn thesis Eq. 5.11**, whose bracket
   is algebraically identical.
2. A **published pillar variant exists**: Chandra thesis Eq. 5.11 (= Chandra &
   Yang, *Langmuir* 25, 10430 (2009)), for **circular** pillars in a square
   lattice collapsing as 4-pillar diagonal clusters. It is the only closed-form
   capillary pillar criterion we located.
3. For **square** cross-sections there is a clean geometric identity
   (Section 5.1): a square pillar of side `a` loaded by a full-face bridge
   deflects **exactly** like a Tanaka line of width `a`. So the line criterion
   transfers to square pillars with **no change of prefactor** — only the choice
   of gap and the one-sided-load assumption are ours.
4. **Numerical result for the device geometry** (h = 700 nm, p = 400 nm, water
   at 20 °C, worst-case θ = 0°): required modulus for stability peaks at
   **E_crit ≈ 13.8 GPa at a = 80 nm** and falls to ≈ 3–3.5 GPa at 160–320 nm.
   Measured amorphous ALD TiO2 is **107–156 GPa**. Margin **≈ 8–11×** at the
   tightest width. **Prediction: no capillary collapse anywhere in 80–320 nm.**
5. The threshold is crossed at **a ≈ 33–37 nm** (h = 700 nm) or at
   **h ≈ 1.2–1.3 µm** (a = 80 nm, AR ≈ 15–16). Those are the numbers to quote
   as "the collapse edge".
6. Two large caveats: (a) the criterion is **elastic-instability only** — no
   sourced treatment of fracture or of base adhesion for TiO2-on-fused-silica
   was found; (b) Chandra's own experiment shows the isolated-bridge criterion
   **over-predicts the required modulus by ≈14×** for 2-D pillar arrays dried
   from a *continuous* liquid body, i.e. it is conservative in that regime.

---

## 1. FORMULA BOX — published expressions, verbatim

### 1.1 Laplace pressure of a bridge between two parallel walls

> **F1** `R = ws / (2 cosθ)`  — Mack (2006), Eq. (1), verbatim
> **F2** `ΔP = 2 γ cosθ / ws` — Mack (2006) Eq. (2) denominator; identical to
> Chandra thesis Eq. 1.8 and Eq. 2.6 (`ΔP = P1 − P2 = 2γ_lv cosθ / x`)

Chandra thesis Sec. 2.3 states the two principal radii of the bridge between
parallel plates are `x/(2cosθ)` and `∞`, which is the derivation of F2.

Full interfacial-energy minimisation gives an **additional line term** that
Tanaka's criterion retains but the simple Laplace argument drops
(Chandra thesis Eq. 2.11, verbatim):

> **F3** `F = −dU/dx = −(2 γ_lv cosθ / x)·L·H − (γ_lv sinθ)·L`
> with `H = h + (x / 2cosθ)(1 − sinθ)`

The first term is Laplace pressure × wetted area `L·H`; the second is the
surface-tension pull on the contact line of length `L`. **This is why `sinθ`
appears in the collapse criterion alongside `cosθ`.**

### 1.2 Capillary load and beam bending for a LINE (Tanaka, via Mack 2006)

Verbatim from Mack (2006), *The Lithography Expert*, November 2006:

> **F4** (Eq. 2) `F = σH/R = 2 σ H cosθ / ws`
> — "For a resist thickness H, the resulting capillary force per unit length of line"
>
> **F5** (Eq. 3) `δ = (3/2) (F/E) (H/wl)³`
>
> **F6** (Eq. 4) — **THE CLASSIC CRITERION**
> ```
> E/σ  ≤  (4√3 · Al³ / ws) · [ 3 As cosθ + sinθ + √( 9 As² cos²θ + 6 As cosθ sinθ ) ]
> ```
> "where As is the aspect ratio of the space and Al is the aspect ratio of the
> line (H/wl)". Collapse occurs when the inequality holds (E/σ **below** the RHS).
>
> **F7** (Eq. 5) small-angle form, valid "If the contact angle is less than
> about 80º and the aspect ratio of the space is high":
> ```
> E/σ  ≤  (8√3 · Al³ / ws) · [ 3 As cosθ + sinθ ]
> ```

**Equivalent "critical modulus" form** (pure algebra on F6, multiplying the
bracket through by `ws`; this rearrangement is ours, the content is Mack's):

```
E_crit(line) = 4√3 · γ · h³ · B / ( wl³ · ws² )

B ≡ 3 h cosθ + ws sinθ + √( 9 h² cos²θ + 3 h ws sin 2θ )        [units: length]
```
(using `6 As cosθ sinθ · ws² = 3 h ws sin2θ` after the multiply-through, since
`sin2θ = 2 sinθ cosθ`). **Collapse if E < E_crit.**

### 1.3 Capillary criterion for CIRCULAR PILLARS (Chandra thesis Eq. 5.11)

Verbatim, thesis p. 77 (visually verified from the PDF at 170 dpi — `pdftotext`
drops the radicals and the √2):

> **F8** "the critical elastic modulus, E'crit, in case of Laplace pressure
> difference due to isolated capillary bridge between the four pillars is given
> by^{8,15}"
> ```
>            128 γ h³ ( 3h cosθ + w sinθ + √( 9h² cos²θ + 3hw sin(2θ) ) )
> E'crit =  ───────────────────────────────────────────────────────────
>                                3 π d³ w²
> ```
> "Here, `w = √2 p − d` is the spacing between the diagonally opposite
> micropillars."

Refs 8 and 15 of that chapter are Chandra/Taylor/Yang, *Soft Matter* 4, 979
(2008) and **Tanaka/Morigami/Atoda, JJAP 32, 6059 (1993)** — i.e. this *is* the
published pillar-ised Tanaka criterion. Its bracket `B` is **identical** to the
line bracket above; only the prefactor and the gap definition differ:

| geometry | prefactor | stiffness term | gap that enters |
|---|---|---|---|
| line (F6) | `4√3 = 6.928` | `wl³` | `ws` = space width (face-to-face) |
| circular pillar (F8) | `128/(3π) = 13.581` | `d³` | `w = √2p − d` (diagonal) |

**Transcription verified numerically:** Chandra states that for h = 9 µm,
d = 0.75 µm, p = 1.5 µm, θ = 60° "E'crit is calculated as 27 GPa". Our
transcription with γ = 72.7 mN/m gives **26.6 GPa** (−1.5 %). ✔

### 1.4 The competing branch — lateral meniscus interaction (Chandra Eq. 5.9)

Verbatim, thesis pp. 74–75:

> **F9** `Ecrit = (32 √2 γ cos²θ h³ / (3 d⁴)) · f(r)`, `r = p/d`
> `f(r) = 1/(r − k) · ( √(2/(k²−1)) + √(1/(2k²−1)) )`   (Eq. 5.10a)
> `r = (1/k)·[ (√2(k²−1)^(−1/2) + (2k²−1)^(−1/2)) / (√2(k²−1)^(−3/2) + 2(2k²−1)^(−3/2)) ] + k` (Eq. 5.10b)
> Fig. 5.4b: **f(2) ≈ 3.2**.

This applies when the pillars are still surrounded by a **continuous** liquid
body (no isolated bridges). Chandra states the resulting Ecrit for the same
geometry is **2 GPa**, i.e. ~13× smaller than the bridge criterion; our
transcription reproduces 2 GPa with f(2) ≈ 3.1–3.2. ✔

> ⚠ **Our observation, flagged:** as printed, thesis Eq. 5.3 has units of N/m
> while Eq. 5.7 has units of N; they cannot both be right. Eq. 5.9 is
> dimensionally sound (`γ h³/d⁴` → Pa) and reproduces the paper's own number, so
> we quote 5.9 and do **not** use 5.3/5.7 quantitatively.
> `[verify Eqs. 5.3/5.7 against Chandra & Yang, Langmuir 25, 10430 (2009)]`

### 1.5 Adhesive (dry-state) collapse of pillars — the SQUARE-pillar formula

Verbatim, Chandra thesis Eqs. 4.1–4.3 (p. 56, visually verified), attributed to
Glassmaker/Jagota/Hui/Kim, *J. R. Soc. Interface* **1**, 23 (2004) (ref. 8) and
to ref. 26 for ground collapse:

> **F10** ground collapse: `E*g = 2^(11/2) · 3^(3/4) · (1−ν²)^(1/4) · h^(3/2) · W / (π d)^(5/2)`
> — "where, W is work of adhesion and ν is the Poisson ratio"
>
> **F11** lateral collapse, **square pillars**: `E*L = 8 h⁴ γs / (3 w² d³)`
>
> **F12** lateral collapse, circular pillars:
> `E*L = 2^(1/4) · 32 h³ γs (1−ν²)^(1/4) / ( 3^(3/4) π d^(5/2) w^(3/2) )`
>
> "where γs is the surface energy of the pillar material… `w` is the spacing
> between two adjacent pillars" (thesis Ch. 4 convention, p. 53).

These are **adhesion**, not capillarity: they govern whether pillars that have
been brought into contact **stay** stuck after the liquid is gone (stiction).
They are the correct check for the post-dry state, and F11 is the only
square-cross-section pillar formula we found published.

### 1.6 Symbols and units

| symbol | meaning | SI unit | note |
|---|---|---|---|
| `γ`, `σ` | liquid–vapour surface tension | N m⁻¹ | Mack writes σ, Chandra γ |
| `θ` | contact angle of the rinse liquid on the structure sidewall | rad (quoted in °) | **not** on the substrate |
| `h`, `H` | structure height | m | |
| `wl`, `d`, `a` | line width / pillar diameter / square-pillar side | m | |
| `ws`, `g` | face-to-face gap (space width) | m | line criterion |
| `w` | in Eq. 5.11: **diagonal** spacing `√2p − d`; in Ch. 4: adjacent spacing | m | **two different conventions — do not conflate** |
| `p` | pitch (centre-to-centre) | m | |
| `As = h/ws`, `Al = h/wl` | space and line aspect ratios | – | |
| `E` | Young's modulus of the structure material | Pa | |
| `ν` | Poisson ratio | – | |
| `γs` | **solid** surface energy of the structure material | J m⁻² | adhesion branch only |
| `W` | work of adhesion structure↔substrate | J m⁻² | ground-collapse branch |
| `δ` | lateral deflection of the top of the structure | m | |
| `ΔP` | Laplace pressure deficit in the bridge | Pa | |
| `B` | `3h cosθ + gap·sinθ + √(9h²cos²θ + 3h·gap·sin2θ)` | m | the shared bracket |

---

## 2. SOURCED CONSTANTS TABLE

| quantity | value | conditions | source | status |
|---|---|---|---|---|
| Water surface tension, 20 °C | **72.74 ± 0.36 mN/m** (calc. 72.74) | pure water–vapour, ITS-90 | IAPWS R1-76(2014), Table 1 | **VERIFIED, archived** |
| Water surface tension, 25 °C | **71.98 ± 0.36 mN/m** (calc. 71.97) | " | IAPWS R1-76(2014), Table 1 | **VERIFIED, archived** |
| Water surface tension, 15 / 30 °C | 73.49 / 71.19 mN/m | " | IAPWS R1-76(2014), Table 1 | **VERIFIED, archived** |
| σ(T) correlation | `σ = B τ^μ (1+bτ)`, τ=1−T/Tc, Tc=647.096 K, B=235.8 mN/m, b=−0.625, μ=1.256 | 0.01 °C → Tc | IAPWS R1-76(2014), p. 3 | **VERIFIED (page rendered)** |
| Water CA on TiO2, dark, single crystal | **74°** | rutile single crystal, before UV | Janczarek/Hupka/Kisch, PPMP 40, 287 (2006), Table 1 — *relay of Wang et al. 1999* | **RELAY** (verbatim from a fetched PDF, but the measurement is Wang's) |
| Water CA on TiO2, after UV | **0°** | 40 mW cm⁻² for 30 min | same, Table 1 + text | **RELAY** |
| Water CA on TiO2 thin film, dark-stored | **53–58°** | TiO2 on Pilkington Activ™ glass, own measurement | Janczarek et al. 2006, Tables 3–5 text | **VERIFIED, archived** |
| UV-amphiphilic mechanism | qualitative only; **no number in the paper** | rutile TiO2(110) + TiO2-coated glass | Wang, Hashimoto, Fujishima et al., *Nature* **388**, 431–432 (1997), DOI 10.1038/41233 | abstract verbatim; **no CA value obtainable** |
| Water CA on **ALD** TiO2 | — | — | — | **NOT FOUND this session** |
| E, amorphous ALD TiO2, 80 °C, 97 nm | **156 ± 19 GPa** (nanoind.) / **107.4 ± 0.4 GPa** (LSAW) | TiCl4/H2O ALD on Si, amorphous | Ylivaara et al., *Thin Solid Films* **732**, 138758 (2021), Table 3 | **VERIFIED, CC-BY full text archived** |
| E, amorphous ALD TiO2, 110 °C, 98 nm | **152 ± 5 GPa** / **126.7 ± 0.5 GPa** | " | same | **VERIFIED** |
| E, amorphous ALD TiO2, 150 °C, 98 nm | **149 ± 4 GPa** / **129.0 ± 0.5 GPa** | " | same | **VERIFIED** |
| E, anatase ALD TiO2, 300 °C, 102 nm | **165 ± 16 GPa** / **166.5 ± 0.6 GPa** | anatase | same | **VERIFIED** |
| Hardness, amorphous ALD TiO2 | **6.7–7.3 GPa** | 80–150 °C, ~100 nm | same, Table 3 | **VERIFIED** |
| Poisson ratio, ALD TiO2 | **ν = 0.27** | used by Ylivaara for both methods | Ylivaara Sec. 2, citing Borgese et al., *Surf. Coat. Technol.* **206**, 2459 (2012) | **VERIFIED as the value used**; Borgese itself not read |
| E, ALD TiO2 (Borgese direct) | ~151 GPa | ALD TiO2, XRD-based | Borgese et al. 2012, DOI 10.1016/j.surfcoat.2011.10.050 | **[unquoted — verify on next use]** |
| Residual stress, ALD TiO2 | 230–967 MPa tensile | 80–350 °C, 19–370 nm | Ylivaara Table 3 | **VERIFIED** |
| ALD TiO2 density (context) | 3.25–3.68 g cm⁻³ (38–125 °C) | TiCl4/H2O | Piercy/Leng/Losego, JVST A 35, 03E107 (2017) — repo entry `piercy-2017-ald-tio2-density` | pre-existing repo entry |
| γs (solid surface energy) of TiO2 | — | — | — | **NOT SOURCED** — blocks F11/F12 |
| W (work of adhesion) TiO2↔fused silica | — | — | — | **NOT SOURCED** — blocks F10 |
| Fracture strength / toughness, amorphous TiO2 film | — | — | — | **NOT SOURCED** |
| IPA surface tension (~21–23 mN/m at 20–25 °C) | — | — | — | **[unquoted — verify on next use]**; no primary source fetched |
| Resist E for context | 2–6 GPa (ArF low, Novolak high); resist–water CA 50–70° | photoresists | Mack (2006) | VERIFIED (context only, not transferable to TiO2) |

---

## 3. UNIT ARITHMETIC (worked, at the tightest geometry)

Take the recommended square-pillar form (Section 5.1), a = 80 nm, h = 700 nm,
p = 400 nm ⇒ face-to-face gap `g = p − a = 320 nm`, θ = 0°, γ = 72.74 mN/m.

**Step 1 — the bracket B** (units of length):
```
B = 3h cos0 + g sin0 + √(9h² cos²0 + 3 h g sin 0)
  = 3(700 nm) + 0 + √(9 · (700 nm)²)
  = 2100 nm + 2100 nm  =  4200 nm  =  4.200 × 10⁻⁶ m
```

**Step 2 — numerator** `4√3 · γ · h³ · B`:
```
4√3 = 6.9282
γ h³ = 0.07274 N m⁻¹ × (700 × 10⁻⁹ m)³ = 0.07274 × 3.43 × 10⁻¹⁹ m³
     = 2.4950 × 10⁻²⁰ N m²
× B  = 2.4950 × 10⁻²⁰ N m² × 4.200 × 10⁻⁶ m = 1.0479 × 10⁻²⁵ N m³
× 6.9282                                    = 7.2600 × 10⁻²⁵ N m³
```

**Step 3 — denominator** `a³ g²`:
```
a³ = (80 × 10⁻⁹ m)³ = 5.120 × 10⁻²² m³
g² = (320 × 10⁻⁹ m)² = 1.024 × 10⁻¹³ m²
a³ g² = 5.2429 × 10⁻³⁵ m⁵
```

**Step 4 — divide; check the units:**
```
[N m³] / [m⁵] = N m⁻² = Pa                                        ✔
E_crit = 7.2600 × 10⁻²⁵ / 5.2429 × 10⁻³⁵ = 1.3847 × 10¹⁰ Pa
       = 13.85 GPa
```

**Step 5 — compare:** amorphous ALD TiO2 `E = 107–156 GPa`
⇒ margin **7.8× (LSAW low end) to 11.0× (nanoindentation)**. Stable.

---

## 4. What is PUBLISHED vs what is OURS

| element | status |
|---|---|
| ΔP = 2γcosθ/gap for a bridge between parallel walls | **published** (Mack Eq. 1–2; Chandra Eq. 1.8, 2.6) |
| the extra `γ sinθ` contact-line term | **published** (Chandra Eq. 2.11) |
| the line stability criterion F6/F7 and its bracket | **published** (Tanaka 1993 via Mack 2006 verbatim) |
| the circular-pillar 4-cluster criterion F8 with `w = √2p − d` | **published** (Chandra thesis Eq. 5.11 = Chandra & Yang, Langmuir 2009) |
| the meniscus-interaction alternative F9 | **published** (Chandra Eq. 5.9) |
| square-pillar **adhesive** lateral collapse F11 | **published** (Glassmaker 2004, relayed verbatim by Chandra Eq. 4.2) |
| square-pillar **capillary** criterion | **NOT PUBLISHED — ours** (Section 5.1) |
| choice of face-to-face vs diagonal gap for a square lattice | **ours** (Section 5.2) |
| the one-sided (asymmetric-drying) load assumption for a pillar | **ours**, following Mack's explicit two-line worst-case argument |
| all numbers in Sections 3, 6 | **ours**, computed from the published formulas + sourced constants |

---

## 5. OUR ADAPTATION — how to apply this to square pillars at 400 nm pitch

### 5.1 The geometric identity that makes the line criterion exact for square pillars

The line criterion depends on the material and geometry only through the
combination (bending stiffness) ÷ (width over which the pressure acts). Compare:

**Line**, per unit length of line (take L = 1 m):
- distributed lateral load per unit height `q = ΔP × 1 m` [N m⁻¹]
- second moment `I = (1 m)(wl³)/12` [m⁴]
- tip deflection of a cantilever under a uniform load: `δ = q h⁴ / (8 E I)`
  ⇒ `δ = 1.5 · ΔP · h⁴ / (E wl³)`

**Square pillar** of side `a`, bridge covering one full face:
- `q = ΔP × a` [N m⁻¹]
- `I = a⁴/12` [m⁴]
- `δ = ΔP a h⁴ / (8 E a⁴/12) = 1.5 · ΔP · h⁴ / (E a³)`

**Identical.** (Note also that a square section has `I = a⁴/12` about *every*
centroidal axis, so the result does not depend on the bending direction.)
Consequently the Tanaka/Mack criterion F6 applies to a square pillar with
`wl → a` and **no change of prefactor**. This is a derivation, not a citation;
but the only new physical assumption is *which* gap and *how many* faces are
loaded — see 5.2/5.3.

Consistency check that the deflection formula matches the source: Mack's Eq. 3
`δ = (3/2)(F/E)(H/wl)³` with `F = q·H` is exactly `δ = q h⁴/(8EI)`, `I = wl³/12`. ✔

### 5.2 Which gap enters — two options, both reported

For a square lattice of pitch `p = 400 nm` and pillar side `a`:

- **face-to-face (nearest-neighbour) gap** `g = p − a`: 320 nm at a = 80 nm down
  to **80 nm at a = 320 nm**. This is the smallest gap present, hence the
  largest Laplace suction. **Recommended, and numerically the more conservative
  of the two at every width in our range.**
- **diagonal gap** `w = √2 p − a` (Chandra's convention for the 4-pillar
  cluster): 485.7 nm at a = 80 nm down to 245.7 nm at a = 320 nm.

Two criteria, both evaluated below:

```
(A)  E_crit,NN   = 4√3 · γ h³ B(h, g) / ( a³ g² ),      g = p − a
                   [line criterion, exact for square section per 5.1]

(B)  E_crit,diag = 8 · γ h³ B(h, w) / ( a³ w² ),        w = √2 p − a
                   [Chandra Eq. 5.11 with the circular stiffness term
                    πd³/64 replaced by the square-section a³/12 — OURS]
```
The prefactor 8 in (B) comes from rewriting Eq. 5.11 as
`E'crit = (2/3)·γh³B/((I/b)·w²)` with `I/b = πd³/64` for a circle, then
substituting `I/b = a³/12` for a square: `(2/3)·12 = 8`. **This constant
transfer is ours and is not validated by any source.** Note it differs from the
line prefactor 4√3 = 6.93 by a factor 2/√3 = 1.155 — i.e. Chandra's 4-pillar
diagonal derivation is not exactly the line derivation, and we do not know
which of the two is right for a square lattice. Use (A) as the working answer
and (B) as the sensitivity.

### 5.3 The assumptions we are making (state these whenever the number is used)

1. **One-sided load.** A pillar with symmetric bridges on all four faces feels
   zero net force; collapse requires symmetry breaking. Both criteria assume the
   asymmetric worst case. Mack states this explicitly for lines: "The two-
   dimensional model used here is really a worst case."
2. **Bridge spans the full pillar height.** Chandra notes bridges that form late
   in the dry-out sit "near the base of the micropillars and thus will exert
   much less torque" — so this is an upper bound on the load.
3. **Parallel-plate ΔP.** A bridge between two finite square faces has a second
   principal radius that is *not* infinite; F2 ignores it. Chandra's published
   pillar criterion makes the same approximation. Unquantified.
4. **Rigid base, no fracture, no plasticity.** Elastic cantilever only.
5. **θ is the contact angle on the pillar sidewall**, i.e. on plasma-processed
   TiO2, not on an as-deposited film.

### 5.4 Contact angle must be carried as a band, not a value

`cosθ` is the dominant sensitivity. Sourced band for TiO2:

| condition | θ | cosθ |
|---|---|---|
| UV-activated / freshly plasma-cleaned | **0°** | 1.00 |
| dark-stored TiO2 thin film | **53–58°** | 0.60–0.53 |
| dark/aged rutile single crystal | **74°** | 0.28 |

A post-etch, post-ash TiO2 sidewall is plausibly at or near the **θ ≈ 0** end —
which is also (per Mack's Q7: the worst case is at `tan⁻¹(1/3As) < 10°`) the
worst case. **Report the threshold as a band from θ = 0° to θ = 74°, and use
θ = 0° for any go/no-go call.** The θ = 0 → 74° swing moves `E_crit` by ≈ 2.4×.

---

## 6. NUMBERS FOR THE DEVICE GEOMETRY

h = 700 nm, p = 400 nm (square lattice), γ = 72.74 mN/m (water, 20 °C).
Values are the **modulus required for stability**; the pillars are stable when
the film modulus exceeds them.

### 6.1 E_crit vs pillar width (GPa)

| a [nm] | g = p−a [nm] | w = √2p−a [nm] | (A) NN, θ=0° | (A) NN, θ=74° | (B) diag, θ=0° | (B) diag, θ=74° |
|---|---|---|---|---|---|---|
| 80 | 320 | 485.7 | **13.85** | 5.66 | 6.94 | 3.28 |
| 100 | 300 | 465.7 | 8.07 | 3.24 | 3.87 | 1.80 |
| 120 | 280 | 445.7 | 5.36 | 2.11 | 2.44 | 1.12 |
| 160 | 240 | 405.7 | 3.08 | 1.16 | 1.24 | 0.55 |
| 200 | 200 | 365.7 | 2.27 | 0.82 | 0.78 | 0.33 |
| 240 | 160 | 325.7 | 2.05 | 0.71 | 0.57 | 0.24 |
| 280 | 120 | 285.7 | 2.30 | 0.75 | 0.47 | 0.19 |
| 320 | 80 | 245.7 | 3.46 | 1.08 | 0.42 | 0.16 |

Compare with **E(a-TiO2) = 107–156 GPa**. Worst case (a = 80 nm, θ = 0°,
criterion A) leaves a margin of **7.8×** against the lowest measured modulus and
**11.0×** against the highest. **No collapse predicted at any width in 80–320 nm.**

Note the non-monotonicity: `E_crit` has a minimum near a ≈ 240 nm. Narrow
pillars are compliant (a³ in the denominator) and wide pillars have a small gap
(g² in the denominator); the safest width in this layout is around 240–280 nm.

### 6.2 Where the threshold actually is

**Critical width at h = 700 nm** (collapse if a < a_c):

| criterion | θ | E = 152 GPa | E = 107.4 GPa |
|---|---|---|---|
| (A) NN | 0° | a_c = **32.8 nm** (AR 21.3) | a_c = **37.2 nm** (AR 18.8) |
| (A) NN | 74° | a_c = 24.4 nm (AR 28.7) | a_c = 27.6 nm (AR 25.4) |
| (B) diag | 0° | a_c = 26.7 nm (AR 26.2) | a_c = 30.1 nm (AR 23.3) |
| (B) diag | 74° | a_c = 20.9 nm (AR 33.4) | a_c = 23.6 nm (AR 29.7) |

**Critical height h_c (nm) [aspect ratio] at p = 400 nm** — collapse if h > h_c:

| a [nm] | (A) θ=0°, E=152 | (A) θ=0°, E=107.4 | (A) θ=74°, E=152 | (B) θ=0°, E=152 | (B) θ=0°, E=107.4 |
|---|---|---|---|---|---|
| 80 | 1274 [15.9] | 1168 [14.6] | 1676 [21.0] | 1514 [18.9] | 1388 [17.4] |
| 120 | 1615 [13.5] | 1481 [12.3] | 2155 [18.0] | 1966 [16.4] | 1803 [15.0] |
| 160 | 1856 [11.6] | 1701 [10.6] | 2496 [15.6] | 2328 [14.5] | 2134 [13.3] |
| 200 | 2003 [10.0] | 1836 [ 9.2] | 2709 [13.5] | 2612 [13.1] | 2395 [12.0] |
| 240 | 2054 [ 8.6] | 1883 [ 7.8] | 2790 [11.6] | 2827 [11.8] | 2592 [10.8] |
| 280 | 1997 [ 7.1] | 1831 [ 6.5] | 2722 [ 9.7] | 2972 [10.6] | 2725 [ 9.7] |
| 320 | 1802 [ 5.6] | 1652 [ 5.2] | 2464 [ 7.7] | 3046 [ 9.5] | 2793 [ 8.7] |

The 700 nm stack sits at roughly **half** the critical height everywhere. Because
`B ∝ 3h` at θ = 0, `E_crit ∝ h⁴`, so a modulus margin of M corresponds to a height
margin of only M^(1/4). Concretely, under the most conservative column
(criterion A, θ = 0°, E = 107.4 GPa) the height margin `h_c/700 nm` runs from
**1.67×** (a = 80 nm) to **2.69×** (a = 240 nm), while the modulus margin over the
same widths is 7.8×–52×.

### 6.3 Temperature and rinse-liquid sensitivity

- 20 °C → 25 °C: γ 72.74 → 71.98 mN/m; `E_crit` at a = 80 nm goes 13.85 → 13.70
  GPa (**−1.0 %**). Negligible.
- Replacing water with a ~21.7 mN/m alcohol rinse would give `E_crit` = 4.13 GPa
  at a = 80 nm — a 3.35× relief, exactly the ratio of surface tensions (`E_crit`
  is linear in γ). **The IPA surface tension used here is
  `[unquoted — verify on next use]`.** Tanaka's abstract names precisely this
  lever: "the use of a low-surface-tension rinse liquid… is effective."

---

## 7. PUBLISHED EXPERIMENTAL DATA — what exists, and what does not

### 7.1 TiO2 or oxide nanopillars vs aspect ratio

**Not found.** We searched OpenAlex (OA-filtered), Semantic Scholar and the open
web for TiO2 / SiO2 / oxide nanopillar collapse-vs-aspect-ratio datasets and
retrieved **none** in full text. Documented negative result. Nearest analogues:

| dataset | what it says | status |
|---|---|---|
| Chandra thesis Ch. 6 (epoxy/PMMA-co-MMA micropillars, two geometries, dried from liquid) | "the micropillar arrays in both the geometries A and B are stable at elastic modulus of 1.2 GPa… and unstable for elastic modulus of 745 MPa"; Eq. 5.11 predicted 16.7 GPa and 4.1 GPa for the same geometries | **VERIFIED** (archived full text) — **the only quantitative calibration of the pillar bridge criterion we have** |
| Zhang, Lo, Taylor & Yang, *Langmuir* **22**, 8595–8601 (2006) — replica-moulded polymer nanopillars; PDMS collapses above AR ≈ 6, stiffer polyurethane/epoxy above AR ≳ 12 | closest AR-vs-material dataset for pillars | **[unquoted — verify on next use]**; DOI `10.1021/la061372+`, closed on OpenAlex, no OA copy found |
| Si nanopillar wet-clean test structures, AR 20, h = 600 nm, d = 30 nm, pitch 90 nm, HF/DIW/IPA sequence | closest *inorganic* nanopillar collapse test vehicle | **[unquoted]** — abstract-level only (Scientific.Net SSP 314, 167) |
| Ghosh et al., "Preventing the Capillary-Induced Collapse of Vertical Nanostructures", *ACS Appl. Mater. Interfaces* (2022), DOI `10.1021/acsami.1c17781` | Si nanostructure collapse prevention | **NOT FETCHED** — the U. Twente OA mirror returns HTTP 403 to curl (tried twice, with cookies and browser headers) |
| Duan & Berggren, *Nano Lett.* **10**, 3710 (2010) and Duan et al., *Small* (2011), DOI `10.1002/smll.201100892` | HSQ (SiOx-like) 10-nm-scale controlled collapse | closed; **NOT FETCHED** |
| "Dimension and process effects on the mechanical stability of ultra-small HSQ nanopillars", *J. Nanopart. Res.* (2021), DOI `10.1007/s11051-021-05364-5` | HSQ **nanopillar** stability vs dimension — the best available oxide-like pillar dataset | closed on OpenAlex; **NOT FETCHED** — **highest-value next fetch** |

### 7.2 Base adhesion / fracture for TiO2 on fused silica

**No published treatment found.** One timeboxed search returned only unrelated
thin-film adhesion work (W/Pt/Ti on amorphous SiO2; SiO2 on polymers). We have
neither `W` (work of adhesion) nor `γs` (TiO2 surface energy) nor a fracture
strength for amorphous TiO2 films, so **F10, F11 and F12 cannot be evaluated**
and the elastic-vs-fracture regime question (Section 8.1) stays open.
Ylivaara's measured **tensile residual stress of 230–967 MPa** is the only
relevant mechanical datum we have for the film–substrate system, and it is three
orders of magnitude below E, so it does not change the elastic criterion.

---

## 8. APPLICABILITY CAVEATS

### 8.1 Which failure regime?

Three distinct mechanisms appear in this literature and only the first is
covered by F6/F8:

1. **Elastic capillary instability** (Tanaka/Chandra): the meniscus force grows
   faster with deflection than the elastic restoring force ⇒ runaway. This is
   what Sections 3 and 6 compute.
2. **Adhesive stiction** (Glassmaker/Roca-Cusachs, F10–F12): once two surfaces
   touch, they stay stuck if the adhesion energy beats the stored bending
   energy. **Cannot be evaluated — γs and W unsourced.** Because 1 is not
   triggered at our geometry, 2 is moot unless something else brings pillars
   into contact.
3. **Fracture at the base**: for a brittle ceramic at E ~ 150 GPa the pillar may
   snap before it bends far enough to touch. Tanaka's abstract names "a rigid
   and highly adhesive resist material" as the fix, implying base adhesion
   matters, but **no fracture criterion for TiO2 was sourced.** Since our
   deflections are far below the instability point, this branch is also not
   currently binding — but it becomes the controlling branch first if the height
   goes above ~1.2 µm.

### 8.2 Isolated bridge vs continuous liquid — the biggest model risk

Chandra thesis, verbatim (p. 71): *"for 1D array of line patterns, Laplace
pressure argument is applicable because in that geometry isolated liquid between
the lines could exist, resulting in different Laplace pressures."* For **2-D
pillar arrays** he argues (and shows experimentally) that while the array is
still immersed in a continuous liquid body, the driving force is the much weaker
**lateral meniscus interaction**, not the isolated-bridge Laplace pressure —
Eq. 5.6 puts the bridge torque at "at least 12 times greater" for θ = 60°,
AR 10, and Ch. 6 shows Eq. 5.11 over-predicting the required modulus by ≈14×
(geometry A) and ≈3.4× (geometry B) relative to experiment.

**Consequence for us:** the E_crit numbers in Section 6 are a **conservative
upper bound** for a drying sequence in which the array remains wetted by a
continuous film. They become the *right* number only in the late dry-out stage
when isolated bridges actually form — and Chandra notes those late bridges sit
near the base and exert less torque still (Q7). Both corrections point the same
way: **the real threshold is at least as favourable as Section 6 says.**

### 8.3 Which gap enters the meniscus force

- For the **line** criterion it is unambiguously the space width `ws` (Mack's
  Eq. 1–2 and Fig. 2).
- For **pillars**, Chandra uses the **diagonal** spacing `w = √2p − d` because
  4 pillars close on a diagonal cluster. We report both (Section 5.2).
- **Pitch never enters directly** — only through the gap. A 400 nm pitch with
  320 nm pillars (80 nm gap) is a very different capillary problem from a 400 nm
  pitch with 80 nm pillars (320 nm gap).
- Edge pillars at the boundary of a metasurface see an **unbalanced** meniscus
  and are the natural first failures; the line literature makes the same point
  (Mack: the isolated-pair geometry "is probably the worst case"). We have no
  sourced quantitative treatment of the edge case.

### 8.4 Other boundary conditions

- Mack notes short lines "connected to another resist feature at one or both
  ends" are stabilised. Free-standing pillars have **no such support** — the
  cantilever model is the right one for us, with no relief.
- Ylivaara's moduli are measured on **planar ~100 nm films on Si**, not on
  700 nm pillars on fused silica. A columnar/porous 700 nm growth could be
  softer. The nanoindentation-vs-LSAW spread (≈45 % at 80 °C) is already a
  warning that the number is method-dependent.
- The criterion is **linear in γ** and therefore linear in every surfactant or
  solvent choice; supercritical CO2 drying removes the liquid–vapour interface
  altogether (the Namatsu-lineage fix) and voids the criterion.

---

## 9. GAPS AND THE NEXT FETCHES (ranked)

1. **Tanaka, JJAP 32, 6059 (1993), pp. 6059–6064** — read the derivation and
   confirm F6 from the original. Currently `[unquoted]`.
2. **HSQ nanopillar stability, *J. Nanopart. Res.* 23 (2021), DOI
   10.1007/s11051-021-05364-5** — the closest oxide-like *pillar* dataset for a
   sanity check of the criterion at the nanoscale.
3. **Chandra & Yang, *Langmuir* 25, 10430 (2009)** — resolve the Eq. 5.3/5.7
   dimensional inconsistency and confirm Eq. 5.11's provenance.
4. **γs of TiO2 and W for TiO2/fused silica** — unblocks F10–F12 and the
   stiction branch.
5. **Water contact angle on a plasma-processed *ALD* TiO2 sidewall** — the
   single largest sensitivity in the criterion (2.4× across the sourced band).
6. **Namatsu, APL 66, 2655 (1995)** — the k·AR² experimental relation and the
   value of k.
7. **IPA / surfactant rinse surface tension from a primary source.**
8. **Zhang et al., Langmuir 22, 8595 (2006)** — the polymer nanopillar AR data.

---

## 10. PROVENANCE — what was fetched, and how

| source | archived extract | retrieval route |
|---|---|---|
| Mack (2006), *The Lithography Expert*, Pattern Collapse | `thesis_extracts/mack_2006_pattern_collapse_tutor55.txt` | direct curl of `lithoguru.com/scientist/litho_tutor/Tutor55 (Nov 06).pdf` |
| Chandra PhD thesis (UPenn) | `thesis_extracts/chandra_2010_upenn_thesis_capillary_pillars.txt` | **DSpace-7 API**: `repository.upenn.edu/server/api/core/bitstreams/16dc42c7-5336-4f5a-8d0f-369a8a0cd218/content` (the `/bitstreams/<uuid>/download` path returns an HTML shim). Eqs. 4.1–4.3, 5.7–5.11 and Fig. 5.4b verified by rendering PDF pages 67/70/88/89/91 at 170 dpi. |
| Ylivaara et al., Thin Solid Films 732, 138758 (2021), CC-BY | `thesis_extracts/ylivaara_2021_ald_tio2_mechanical.txt` | **Aaltodoc API**: `aaltodoc.aalto.fi/server/api/core/bitstreams/bbf794b3-ad46-4621-9bbe-1588075509e3/content` (ScienceDirect 403; JYX metadata-only). OA location found via `api.openalex.org/works/doi:…` |
| IAPWS R1-76(2014) | `thesis_extracts/iapws_r1_76_2014_water_surface_tension.txt` | direct curl of `iapws.org/public/documents/CH-L9/Surf-H2O-2014.pdf`; equation page rendered at 150 dpi |
| Janczarek/Hupka/Kisch, PPMP 40, 287 (2006) | `thesis_extracts/janczarek_2006_tio2_hydrophilicity.txt` | direct curl of `journalssystem.com/ppmp/pdf-79371-15436` |
| Tanaka 1993 (record + abstract only) | `thesis_extracts/tanaka_1993_jjap_iopscience_record.txt` | curl of the IOPscience article page; body paywalled |
| Que & Gianchandani, JVST B 18, 3450 (2000) | `thesis_extracts/que_gianchandani_2000_jvstb_resist_mechanics.txt` | direct curl from the U. Michigan group page (context: resist E = 7.0 ± 0.5 GPa UV6, 21.5 ± 5.0 GPa APEX-E) |

**Blocked this session:** AIP/pubs.aip.org (Cloudflare 403), ScienceDirect (403),
Springer/link.springer.com (empty body), research.utwente.nl OA mirror (403 even
with cookies + browser headers), nature.com bronze-OA PDF (redirects to paywall),
RSC, ACS.

Library entries created alongside this document (not committed):
`library/tanaka-1993-jjap.md`, `library/mack-2006-pattern-collapse.md`,
`library/chandra-2010-thesis.md`, `library/glassmaker-2004-jrsi.md`,
`library/namatsu-1995-apl.md`, `library/iapws-2014-water-surface-tension.md`,
`library/ylivaara-2021-ald-tio2-mechanical.md`,
`library/janczarek-2006-tio2-hydrophilicity.md`,
`library/wang-1997-nature-amphiphilic.md`.
