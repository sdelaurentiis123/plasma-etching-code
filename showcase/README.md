# showcase — reactor-to-feature visualization page

One self-contained dark HTML page (`index.html`) showing the engine across
scales: a live canvas animation of the Turner–Chabert moving RF sheath
(closed-form field re-implemented in JS, constants receipt-checked against the
frozen audit), an interactive knob panel over the three solved Oxford
PlasmaPro 80 power nodes, the feature-evolution movies, and validation charts
drawn from frozen curated boards.

Rebuild:

```bash
python showcase/build_data.py     # curated boards -> data/showcase_data.json
                                  # (fails closed if the sheath replay does
                                  #  not match the frozen audit to 1e-9)
python showcase/render_movies.py  # -> media/{reactor_scene,sheath_scene,
                                  #           flyby}.mp4 (~3 min, CPU)
python showcase/build_page.py     # template + data + media -> index.html
```

The movie series: `reactor_scene.mp4` (chamber scale, solved Oxford 105 W
state), `sheath_scene.mp4` (tracer ions through the replayed moving RF
field), and `flyby.mp4` (title → reactor → sheath → feature, one continuous
dive; embedded at the top of the page).

Media inlined from `viz/etch_hole.gif`, `viz/etch_trench.gif`, and
`../bench/out/hero3d.mp4`. Every number on the page traces to
`results/curated/` through `build_data.py`; nothing is typed into the template
except section prose. The race movie is deliberately not included: its frames
carry the retired 47x figure (honest wall-clock is ~14x, quoted as a stat
card instead).
