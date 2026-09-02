# Draft email to Arun — 2026-09-02

Subject: CLKGATE X1 etch — what we found, what the sim says, what we need from you

Hi Arun,

Quick update on the clock-gate STL you sent. Short version: the geometry is fully pinned, we have a first-shot etch recipe with a predicted depth, and the one thing we can't compute from the file alone is how your printed polymer walls handle fluorine. Five measurements from a single sacrificial run would close that.

**What we did with the file**

- Recovered the source exactly. The STL is the Nangate/FreePDK45 CLKGATE_X1 Metal-1 cell, enlarged 40x: 98.8 x 62.8 um in plan, 30 um tall, narrowest gap 2.6 um, tracks 2.8 um. Every polygon matches to under 0.2 nm after scaling, so there is no ambiguity left in scale or footprint. (The "25x" in the filename is the objective preset, not a geometry scale.)
- Repaired nine zero-area export facets; the mesh is now closed and manifold with no geometric change.
- Ran the exact ten-polygon footprint through our deterministic transport operator. No proxy mask, no redrawn rectangles.

**What the physics says**

- The 30 um print is a tall entrance mask (aspect ratio about 11.5 at the narrowest gap). A narrow ion beam mostly clears it: about 92% of direct ion flux reaches the floor on average, 73% to 100% depending on position in the layout.
- Thermal radicals do not. Only about 7% of direct F/O reaches the floor. So the etch inside your openings is set by what F and O do after they hit the polymer sidewalls, which is a material property of your resist that no published number gives us.
- First-shot recipe we'd propose (Miao et al., Oxford Plasmalab 100 cryo ICP): -110 C, 8 mTorr, 1000 W ICP, 10 W platen, 52/8 sccm SF6/O2, 2 min. Open-field silicon depth 3.5 to 3.9 um; inside the clock-gate openings roughly 2.5 to 3.0 um under example wall-return assumptions. Polymer loss should be small, about 0.25 um, if your print behaves like the polymers in that paper. This is a conditional prediction, not a fitted one: nothing from a target result went in.

The attached HTML is a self-contained explorer: exact-mask transport map, cross-sections, time evolution, and the wall-return sensitivity, with a short physics explainer at the top. It runs offline in any browser.

**What we need from you (one run, five numbers)**

1. Polarity and target: raised silicon replica of Metal-1, or Metal-1-shaped trenches? And the silicon depth you want.
2. Which etcher, and whether it can hold -110 C. If not, tell us its range and we re-select the recipe.
3. From the run: DC self-bias or platen voltage, helium backside pressure, and sample position. Whatever the tool logs.
4. Polymer height before and after, plus one blanket or wide-open silicon depth on the same wafer. Those two numbers give us your tool's absolute rate and polymer selectivity without touching the patterned result.
5. A pre-etch image of the print and one post-etch cross-section, which we will only open after we freeze the profile prediction.

Point 5 is how we keep ourselves honest. We just did exactly this with another group's pillar run: prediction frozen from the recipe alone, SEMs opened afterwards, and the failure mechanism came out as predicted. We'd like to do the same on your geometry.

Happy to jump on a call to walk through the explorer.

Stan
