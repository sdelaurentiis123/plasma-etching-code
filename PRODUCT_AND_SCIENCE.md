# petch — product accessibility + science roadmap

What we have: a fast (≈14× ViennaPS-GPU), ViennaPS-accurate, **differentiable**, modular 3D plasma-etch
simulator with a clean `Domain`/`Process`/`Result` API. Two questions: how to make it a *great, accessible
product*, and what *science* it uniquely enables.

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

**6. Lead with differentiability in the pitch.** "Fast" gets attention; "the only *differentiable* etch
simulator" is the moat — it unlocks the science below. The README/tweet should headline both:
*14× faster AND gradient-based (inverse design, calibration) — neither possible in ViennaPS.*

**7. Docs site.** The `docs/` HTML explainers → a real site (GitHub Pages) with API reference +
the physics writeups + the benchmark. Trust + SEO.

Quick wins to bank first: Colab notebook (#1), animate() (#2), migration guide (#4).

---

## Part 2 — the science we can uniquely do

The combination **fast + differentiable** opens experiments that are impossible in ViennaPS (not
differentiable) and impractical by finite-difference (too many slow runs). Two are genuinely novel:

### A. Differentiable inverse design of process recipes  ★ the headline science
Gradient-descend the process knobs (ion flux, energy, O₂ fraction, sticking) to hit a **target profile** —
vertical sidewalls, a depth, zero microtrenching — by backpropagating through the simulator (`wp.Tape`).
*"Given the structure you want, compute the recipe."* ViennaPS can only guess-and-check. This is the demo
that defines the tool. *We already have a 1-parameter inverse-design proof; scale it to the multi-knob,
profile-target case.*

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
