# Orchestration protocol (Opus-led, Fable-advised)

Operating model for the petch build loop.

## Roles
- **Opus (orchestrator + executor)** — owns git/gates/decisions and does the primary coding.
  Dispatches Opus sub-agents for heavy or parallel builds. Reviews all results critically before
  committing. This is the default actor.
- **Fable (advisor)** — invoked via sub-agent (`model: fable`) for: scoping a hard problem,
  critiquing a plan *before* execution, evaluating a research idea / physical hypothesis,
  decomposing an ambiguous question, sanity-checking a surprising result. Rules for Fable:
  high-level analysis only, **concise** (no long write-ups), **no web search**, no strenuous code.
  Fable's output is a recommendation to Opus, not a deliverable.
- **Opus executors** — sub-agents (`model: opus`) for delegated builds. Brief them with: exact
  file targets, a two-sided/quantitative GATE (PASS/fail), the paid-for pitfalls, foreground
  polling for long runs (do NOT rely on background re-invocation), and commit-local-don't-push.

## Loop
1. **Scope** — Opus frames the move; consults Fable on hard/ambiguous ones.
2. **Critique** — Fable stress-tests the plan/hypothesis (kills bad ones before dispatch).
3. **Execute** — Opus builds directly, or dispatches Opus executors in parallel.
4. **Review** — Opus grades results honestly (gate PASS/fail); re-consults Fable on surprises.
5. **Land** — gate green + suite green → commit → push.

## Guardrails (unchanged)
Never touch viennaps-accel/ or plasma_sim/. Every claim = script + numbers + PASS/fail. Refuted
hypotheses documented, not deleted. Keep partner/Resona specifics out of the repo. 26-test suite
must stay green. Commit trailer + Claude-Session line.
