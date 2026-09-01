#!/usr/bin/env python3
"""Assemble the client-facing Oxford page: predicted pillar vs SEM, one-to-one."""
import base64, json
from pathlib import Path
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
def b64(p): return "data:image/png;base64," + base64.b64encode(Path(p).read_bytes()).decode()
t = (HERE / "template.html").read_text()
for k, f in {"__PRED__": HERE / "pred_intact_w185.png", "__SEM_INTACT__": HERE / "sem_intact_dose9.png",
             "__SEM_LOSS__": HERE / "sem_loss_dose9.png"}.items():
    t = t.replace(k, b64(f))
out = HERE / "oxford_tio2_prediction_vs_sem.html"; out.write_text(t)
print(out.relative_to(ROOT), f"{out.stat().st_size/1e6:.2f} MB")
