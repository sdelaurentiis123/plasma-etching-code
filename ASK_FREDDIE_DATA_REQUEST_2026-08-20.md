# Data request to Freddie — draft, 2026-08-20

Status: DRAFT for the owner to send. Nothing here has been sent. The point
of every item is to condition the boundary/surface deck on independent
same-run observables — never to fit the target profile. The SEMs stay the
sealed answer key until the blind board is frozen against them.

---

Subject: NPG80 TiO2 run — a short list of run observables that would let us
grade the blind prediction properly

Hey Freddie,

We have the full blind pillar board running through the model now (all seven
widths, the ion-energy/angle corners, and both selectivity witnesses). Before
we freeze the prediction against your SEMs, a handful of observables from the
actual run would let us pin the boundary without touching the answer key.
In rough priority order:

1. **Achieved DC self-bias** during the etch (the reading on the tool, not
   the setpoint), and anything the tool logs about forward/reflected power.
   Forward watts alone don't determine the ion energy; the self-bias does.

2. **Blanket TiO2 loss**: if any unpatterned area (or a witness piece) went
   through the same run, the measured TiO2 thickness before/after. This is
   the single most valuable number — it separates the reactor dose from the
   surface yield.

3. **Remaining Cr thickness** after the etch, anywhere measurable —
   especially interesting because our board currently says the 45 nm cap is
   marginal: parts of the witness-rate interval consume it entirely before
   20 min, corners first.

4. **Actual GDS dimensions** (drawn widths and pitch) and roughly where the
   sample sat on the electrode (center/edge, orientation).

5. **The rinse/dry procedure** after etch: DI rinse? N2 blowdry? IPA? Time
   between etch and imaging. (This matters because capillary forces at
   drying are one of the two candidate mechanisms for the clustered
   collapse.)

6. **SEMs, when you're ready to unseal them** — most diagnostic single shot:
   a tilted view **across the boundary of a collapse cluster**, where fallen
   and standing pillars are adjacent. Four things discriminate the collapse
   mechanism:
   - tops of fallen pillars: Cr cap still present vs eroded/rounded bare TiO2;
   - fall geometry: pairs/bundles leaning toward each other (capillary) vs
     all one direction (flow/handling) vs random;
   - cluster boundary shape: straight/rectangular (litho writefield) vs
     hugging array edges (loading) vs irregular blobs (film defects);
   - pillar bases: clean break vs undercut/notched feet.

Items 1–5 we can use immediately; item 6 stays sealed until we hand you the
frozen prediction.

---

## Why each item, in one line (internal)

| Item | Conditions | Board use |
| --- | --- | --- |
| self-bias | ion impact energy | replaces the cross-family 276 V transfer |
| blanket TiO2 loss | absolute rate | splits reactor dose from surface yield |
| remaining Cr | mask-survival call | direct test of the exhaustion knife-edge |
| GDS + position | geometry prior | replaces the assumed 400 nm/80–320 board |
| rinse/dry | collapse mechanism | capillary criterion inputs |
| SEM cluster boundary | mechanism fingerprints | discriminates all five candidates |
