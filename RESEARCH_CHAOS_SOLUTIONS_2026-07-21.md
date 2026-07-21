# Solution-space literature synthesis: calibrating a chaotic feature-scale simulator

(Deep-research agent output, 2026-07-21; full citations inline. Companion to the R4 chaos
verdict in the validation protocol and to NEXT_STEPS.md.)

Ranked architecture the combined literatures support:

1. **Analytic / semi-analytic form factors (do FIRST, highest ROI).** Baum-Rushmeier-Winget
   (SIGGRAPH '89) evaluate unoccluded form factors analytically (Stokes' theorem) and use
   sampling ONLY for the occluded fraction; Schroder-Hanrahan (SIGGRAPH '93) give the exact
   polygon-to-polygon closed form. Recast our estimator as analytic-unoccluded x ray-estimated
   visibility fraction (a ratio estimator; Ramamoorthi et al. give the variance theory):
   sampling noise then scales with the bounded [0,1] visibility fraction, not the full kernel
   -- expect 1-2 orders of magnitude variance reduction, collapsing most of the 30-46 nm
   endpoint scatter AT THE SOURCE and diagnostically separating numeric noise from real
   front-feedback physics.
2. **Ensemble-median/quantile observables + stochastic-kriging/hetGP calibration.** The
   stochastic-simulator literature (Ankenman-Nelson-Staum 2010; Binois-Gramacy hetGP +
   replication-vs-exploration, arXiv:1710.03206; Fadikar quantile-emulation, SIAM/ASA JUQ)
   says: the output is a random variable; calibrate to its median/quantiles with planned
   replication. CRN helps gradients but can decorrelate near the bifurcation.
3. **Wide cells restore self-averaging.** Aharony-Harris (PRL 1996) non-self-averaging near
   instabilities is the precise theory of our 20-nm-cell finding; LER correlation lengths are
   8-24 nm, so cells of 3-5x xi (~60-100 nm) restore self-averaging at linear cost -- cheaper
   than large ensembles at fixed narrow cell.
4. **Physically-scaled Mullins/curvature regularization** (Sethian-Adalsteinsson III) with
   regularization length below xi -- damps grid noise, never tuned to hide the instability.
5. **Shadowing (LSS/NILSS/NILSAS) + EKI + randomized-smoothing AD** -- the machinery for
   gradients/calibration of chaotic engines; reserve for the differentiable roadmap. Lesson:
   backprop through ensemble statistics, never a single seed's endpoint.

Field practice notes: Kushner MCFPM fights pseudo-particle roughness with local surface
smoothing + statistical weighting (Huard thesis); nobody characterizes seed-to-seed endpoint
distributions or their effect on model-experiment comparison -- our finding plus the
cell-length-vs-xi convergence study are two publishable gaps.

(Agent's full per-thread writeup with all URLs preserved in the session log; key primaries:
Baum 1989 doi 10.1145/74333.74367; Schroder-Hanrahan graphics.stanford.edu/papers/formfactor;
Hanrahan hierarchical radiosity SIGGRAPH '91; Ankenman Oper.Res. 58(2); arXiv:1710.03206;
SIAM/ASA doi 10.1137/17m1161233; Aharony-Harris PRL 77:3700; Guo-Sawin J.Phys.D 42:194014;
arXiv:1204.0159, 1611.00880, 1801.08674; arXiv:2104.03384.)
