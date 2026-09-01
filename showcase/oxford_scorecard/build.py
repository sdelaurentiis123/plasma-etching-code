#!/usr/bin/env python3
"""Build the Oxford NPG80 TiO2 blind-prediction scorecard page for Freddie.

Reads only committed artifacts: the frozen prediction envelopes, the
digitized SEM measurements, and the original SEM frames (downscaled and
inlined).  The loss classification per dose is the whole-set visual review
recorded below, with the frames it rests on named.
"""
from __future__ import annotations

import base64
import io
import json
from pathlib import Path

from PIL import Image

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SEM = ROOT / "data" / "experimental" / "zhu_2026_tio2_npg80" / "sem_0817"
GDS = ROOT / "results" / "curated" / "zhu_npg80_gds_square_profiles_v1"

# Whole-set visual review (all 80 frames), block-level loss classification.
LOSS = {
    1: ("intact where imaged", "dose1_06–10"),
    2: ("intact where imaged", "dose2_02–04"),
    3: ("isolated missing sites, one partial pillar", "dose3_08–10"),
    4: ("intact where imaged", "dose4_02–07"),
    5: ("intact where imaged", "dose5_03–08"),
    6: ("intact where imaged", "dose6_02–08"),
    7: ("loss patch with a straight boundary; stubs in the patch", "dose7_15d_01–02"),
    8: ("catastrophic: most of the block removed, irregular islands survive", "dose8_02–06"),
    9: ("large loss zones with angular boundaries; stubs; intact zones adjacent", "dose9_15d_01–08"),
}

GALLERY = [
    ("dose3_09.tif", "Dose 3, plan view (2.6 nm/px). Standing pillars are 4-lobed clovers ~240 nm across; two lattice sites are empty and one pillar is partially gone."),
    ("dose9_15d_02.tif", "Dose 9, 15° tilt (5.6 nm/px). Left: intact flat-topped pillars. Right: every site holds a short stub — pillars consumed from the top, not lying down. Sharp boundary."),
    ("dose8_02.tif", "Dose 8, plan view (13.8 nm/px). The most fragile block: most of the array removed, irregular surviving islands."),
    ("dose7_15d_01.tif", "Dose 7, 15° tilt (9.5 nm/px). A loss patch with a near-straight boundary against an intact field."),
    ("dose4_04.tif", "Dose 4, plan view (3.4 nm/px). Intact block, ~222 nm pillars at ~340 nm pitch."),
    ("dose9_15d_01.tif", "Dose 9, 15° tilt, low mag (27.7 nm/px). Whole-block view of the loss zone geometry."),
]


def b64png(path: Path, width: int = 820) -> str:
    im = Image.open(path).convert("L")
    im = im.crop((0, 0, im.width, 683)) if im.height >= 683 else im
    im = im.resize((width, int(im.height * width / im.width)), Image.LANCZOS)
    buf = io.BytesIO(); im.save(buf, "PNG", optimize=True)
    return "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode()


def b64file(path: Path, mime="image/png") -> str:
    return f"data:{mime};base64," + base64.b64encode(path.read_bytes()).decode()


def main():
    summary = json.loads((SEM / "dose_summary.json").read_text())
    env = json.loads((GDS / "profile_envelopes.json").read_text())
    template = (HERE / "template.html").read_text()
    rows = []
    for d in range(1, 10):
        s = summary.get(str(d))
        loss, frames = LOSS[d]
        if s:
            cd = s["cd_equivalent_square_nm"]; xe = s["cd_x_extent_nm"]
            rows.append(f"<tr><td>{d}</td><td class=mono>{cd[0]:.0f} <span class=dim>({cd[1]:.0f}–{cd[2]:.0f})</span></td>"
                        f"<td class=mono>{xe[0]:.0f}</td><td class=mono>{s['pitch_nm_median']:.0f}</td>"
                        f"<td>{loss}</td><td class=dim>{frames}</td></tr>")
        else:
            rows.append(f"<tr><td>{d}</td><td>—</td><td>—</td><td>—</td><td>{loss}</td><td class=dim>{frames}</td></tr>")
    gallery = "".join(
        f'<figure><img src="{b64png(SEM / "originals" / f)}" alt="{f}"><figcaption><b>{f[:-4]}</b> — {cap}</figcaption></figure>'
        for f, cap in GALLERY)
    html = (template
            .replace("__DOSE_ROWS__", "\n".join(rows))
            .replace("__GALLERY__", gallery)
            .replace("__ENVELOPES_PNG__", b64file(GDS / "profile_envelopes.png"))
            .replace("__CONTACT_1__", b64file(SEM / "contact_sheet_1.png")))
    out = HERE / "oxford_tio2_blind_scorecard.html"
    out.write_text(html)
    print(out.relative_to(ROOT), f"{out.stat().st_size/1e6:.2f} MB")


if __name__ == "__main__":
    main()
