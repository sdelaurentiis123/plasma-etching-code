# petch — product accessibility + science roadmap

Current evidence supports a fast feature-evolution codebase, selected ViennaPS transport/parity gates,
several reduced chemistry/process modules, an ALE gradient proof, and an experimental charging/adjoint
stack under active convergence work. It does **not** yet support the blanket claim that the unified 3-D
SiO2 engine is experimentally validated, end-to-end differentiable, or always 14x faster than ViennaPS.
Those are gates to earn, not premises for product copy.

## Product truth contract

The product separates four things that must never be blurred:

1. **Governing solver physics:** conservation, electrostatics, trajectories/collisions given cross
   sections, geometric transport, surface balances, and interface kinematics.
2. **Physical mechanism inputs:** boundary distributions, cross sections, energy-angle-state yield
   tables, reaction networks, dielectric/conductive properties, uncertainty, units, and provenance.
3. **Unknown closures:** constrained calibration, atomistic/ML surrogates with validity diagnostics, or
   explicit refusal. A closure may not silently replace conservation or invent a reaction.
4. **Validation and transfer:** calibrate uncertain inputs on declared cases, predict held-out cases,
   and report numerical, measurement, parameter, and structural-model errors separately.

The first product slice is a SiO2-focused **forward** engine: config and physical closure data in; profile,
surface/charge state, uncertainty, and a validity decision out. Gradients enable calibration and a later
inverse target-profile-to-recipe loop, but do not mean that inverse design is already shipped. Until
surface feedback to a reactor model is closed, the accurate description is “feature-scale with physical
reactor boundary conditions,” not a fully coupled multiscale digital twin.

Material generality is earned only when a second chemistry runs through the unchanged core. Existing
Si-Cl2-Ar+ ALE, SF6/O2, cryogenic, and Bosch modules are valuable starting mechanisms, but several remain
separate reduced models rather than proof of one universal stateful chemistry engine.

---

## Part 1 — make it the most accessible tool (ranked by leverage)

**1. A Colab notebook (the single highest-leverage move).** A "try petch in your browser, free GPU"
notebook: install → run a hole etch → `result.plot()` → run `vps_sweep` → see 14×. For a GPU tool this
removes the #1 barrier (no NVIDIA box needed to evaluate). One link in the README/tweet = anyone can
reproduce the headline in 60 s. *Build: a `notebooks/petch_intro.ipynb`.*

**2. Visual-first output.** Scientists adopt what they can see. We now have `result.plot()` (cross-section)
and `result.save(.vtk)` (ParaView). Next: a `result.animate()` of the etch evolving (record φ each step →
GIF) — that GIF is the tweet. *Build: depth-history → frames → GIF.*

**3. `pip install petch` on PyPI.** Right now it's `pip install -e .`. Real PyPI = one-line install,
discoverability, version pins. *Build: `python -m build` + `twine upload` (after a clean `pip install .`
test in a fresh env).*

**4. The ViennaPS migration guide.** A side-by-side "ViennaPS → petch in 5 lines, 14× faster" page. The
API already mirrors `Domain`/`Process`, so this is mostly a table + a worked example. This is the *adoption
funnel*: their users are the ones who already need this. *Build: `docs/migrate_from_viennaps.md`.*

**5. Examples gallery.** `examples/`: hole, trench, ARDE sweep, the differentiable inverse-design demo —
each a runnable script + a figure. The gallery is the documentation people actually read.

**6. Lead with differentiability only after the gradient gates close.** The defensible near-term wedge is
physical calibration and transfer with explicit uncertainty. “The only differentiable etch simulator”
requires finite-difference agreement through transport, the converged charging fixed point, chemistry,
and profile motion—not only the existing ALE proof.

**7. Docs site.** The `docs/` HTML explainers → a real site (GitHub Pages) with API reference +
the physics writeups + the benchmark. Trust + SEO.

Quick wins to bank first: Colab notebook (#1), animate() (#2), migration guide (#4).

---

## Part 2 — the science we can uniquely do

The combination **fast + differentiable** opens experiments that are impossible in ViennaPS (not
differentiable) and impractical by finite-difference (too many slow runs). Two are genuinely novel:

### A. Differentiable inverse design of process recipes  ★ the headline science
The intended second act is to gradient-descend process knobs to hit a **target profile** —
vertical sidewalls, a depth, zero microtrenching — by backpropagating through the simulator (`wp.Tape`).
*"Given the structure you want, compute the recipe."* The existing proof is limited to the reduced ALE
model. Scaling it to a multi-knob 3-D profile target is open work and must be described separately from
the forward product.

### B. Differentiable calibration to experimental SEM data  ★ the useful science
Fit the model's physical parameters (sticking coefficients, yield prefactors, rate) to a published SEM
cross-section (Belen 2005, de Boer, Gomez) **by gradient descent on the profile mismatch** — instead of
hand-tuning. *"Auto-calibrate the simulator to a wafer dataset."* This is how you'd actually deploy it on
a real process, and it's a real research contribution (differentiable physics + experimental data → fitted
constants). It also turns our honest "not validated to wafers yet" into a *method* for getting there.

### C. ARDE / RIE-lag landscapes at scale
14× speed → sweep the full design space (aspect ratio × flux × energy × feature width × O₂ fraction) and
emit the design-rule response surfaces a process engineer wants. The kind of map that's too slow to make
with ViennaPS. Pairs with the gradients for **global sensitivity**: which knob most controls sidewall
angle? bow? microtrench depth? — read straight off ∂profile/∂param.

### D. Optimal *time-varying* recipes (differentiable optimal control)
Most etches use constant conditions. With gradients through time, optimize a **recipe schedule**
(flux/energy vs time) for a profile constant recipes can't reach — a principled generalization of
Bosch-style alternation. Genuinely new capability.

### E. Uncertainty quantification / process-window mapping
Fast forward model → Monte-Carlo over input tolerances → which process variations blow the CD budget?
Map the robust process window. Cheap once a run is ~1 s.

### F. Studies of the physics ViennaPS omits (as it gets calibrated)
Systematic profile studies of charging (notching/twisting), the bimodal-IEDF tail effect, and
redeposition/taper — each an effect ViennaPS can't model. (Charging needs the proper field solve first.)

**The two to lead with: (A) inverse design and (B) differentiable calibration.** Both are unique to a
fast+differentiable simulator, both make great papers/threads, and (B) is the credible path from
"ViennaPS-accurate" to "real-wafer-accurate."
