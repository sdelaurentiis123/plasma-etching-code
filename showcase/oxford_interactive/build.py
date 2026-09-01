#!/usr/bin/env python3
"""Interactive whole-board page for Freddie: per dose block, his SEM beside our
predicted pillar as a rotatable 3-D mesh, with etch-time scrub and numbers.
All geometry from the frozen exact-GDS board; all measurements from the
digitized SEM set.  Self-contained, no network."""
import base64, io, json
from pathlib import Path
import numpy as np
from PIL import Image, ImageDraw

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
AUD = json.loads((ROOT / "results/curated/zhu_npg80_gds_square_profiles_v1/audit.json").read_text())
SEMD = ROOT / "data/experimental/zhu_2026_tio2_npg80/sem_0817"
SUMM = json.loads((SEMD / "dose_summary.json").read_text())
FRAMES = {f["file"]: f for f in json.loads((SEMD / "manifest.json").read_text())["frames"]}
FILM_NM = 700.0

# dose block -> (nearest board width, representative frames [(file, crop box, caption)], loss label)
BLOCKS = {
    1: (195, [("dose1_07.tif", (300, 60, 1000, 585), "plan view, 10°")], "intact where imaged"),
    2: (195, [("dose2_02.tif", (520, 220, 1000, 580), "plan view, 10°")], "intact where imaged"),
    3: (235, [("dose3_09.tif", (0, 0, 1024, 683), "plan view, 10° — note the empty sites")], "isolated missing pillars"),
    4: (225, [("dose4_04.tif", (0, 0, 1024, 683), "plan view, 10°")], "intact where imaged"),
    5: (195, [("dose5_04.tif", (0, 0, 1024, 683), "plan view, 10°")], "intact where imaged"),
    6: (225, [("dose6_05.tif", (0, 0, 1024, 683), "plan view, 10°")], "intact where imaged"),
    7: (195, [("dose7_07.tif", (0, 0, 1024, 683), "plan view, 10°"), ("dose7_15d_01.tif", (0, 0, 1024, 683), "15° tilt — loss patch with a straight boundary")], "loss patch, stubs inside"),
    8: (185, [("dose8_06.tif", (0, 0, 1024, 683), "plan view, 10° — surviving islands"), ("dose8_02.tif", (0, 0, 1024, 683), "low mag — most of the block removed")], "catastrophic loss"),
    9: (185, [("dose9_15d_02.tif", (0, 40, 400, 340), "15° tilt — intact zone"), ("dose9_15d_02.tif", (560, 40, 960, 340), "15° tilt — loss zone: stubs"), ("dose9_15d_01.tif", (0, 0, 1024, 683), "15° tilt, low mag — zone geometry")], "large loss zones + intact zones"),
}


def crop64(fname, box, width=640):
    im = Image.open(SEMD / "originals" / fname).convert("L").crop(box)
    px = FRAMES[fname]["meta"]["pixel_size_nm"]
    scale = width / im.width
    im = im.resize((width, int(im.height * scale)), Image.LANCZOS)
    d = ImageDraw.Draw(im)
    bar_nm = 200 if px < 8 else 1000
    L = int(bar_nm / px * scale)
    d.rectangle([12, im.height - 18, 12 + L, im.height - 13], fill=255)
    d.text((12, im.height - 34), f"{bar_nm} nm" if bar_nm < 1000 else "1 µm", fill=255)
    buf = io.BytesIO(); im.save(buf, "PNG", optimize=True)
    return "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode()


def profile(width):
    ps = [p for p in AUD["profiles"] if p.get("width_nm") == float(width)
          and p["transport_scenario"]["name"] == "ion_low_tail_0p0" and p["tio2_to_cr_selectivity"] == 14.0]
    p = ps[0]; cs = p["profile"]["cross_section"]; fl = p["profile"]["floor_height_nm"]
    h = [c["height_um"] * 1000 - fl for c in cs]; w = [c["mean_width_nm"] for c in cs]
    pairs = sorted(zip(h, w))
    hs = [a for a, b in pairs if b > 2]; ws = [b for a, b in pairs if b > 2]
    # envelope across all 16 endpoints for this width
    deps = [q["profile"]["etched_depth_nm"] for q in AUD["profiles"] if q.get("width_nm") == float(width)]
    mids = [q["profile"]["middle_cd_nm"] for q in AUD["profiles"] if q.get("width_nm") == float(width)]
    return {"h": hs, "w": ws, "top_h": max(hs) if hs else 0, "depth": p["profile"]["etched_depth_nm"],
            "depth_range": [min(deps), max(deps)], "mid_cd_range": [min(mids), max(mids)],
            "cr_gone": bool(p["cr_mask"]["mask_exhausted_at_center"]), "sidewall": p["profile"]["sidewall_angle_from_wafer_deg"]}


data = {"film_nm": FILM_NM, "pitch_nm": 350.0, "rates_nm_min": [34.125, 43.4667], "selectivity": [14.0, 18.0167],
        "blocks": [], "profiles": {}}
for d, (wmodel, frames, loss) in BLOCKS.items():
    s = SUMM[str(d)]
    data["blocks"].append({"dose": d, "width_model": wmodel, "loss": loss,
                           "cd_meas": s["cd_equivalent_square_nm"], "tip_meas": s["cd_x_extent_nm"][0],
                           "pitch_meas": s["pitch_nm_median"],
                           "frames": [{"img": crop64(f, box), "caption": cap, "file": f[:-4]} for f, box, cap in frames]})
    if str(wmodel) not in data["profiles"]:
        data["profiles"][str(wmodel)] = profile(wmodel)
t = (HERE / "template.html").read_text().replace("__DATA__", json.dumps(data))
out = HERE / "oxford_tio2_interactive.html"; out.write_text(t)
print(out.relative_to(ROOT), f"{out.stat().st_size/1e6:.2f} MB")
