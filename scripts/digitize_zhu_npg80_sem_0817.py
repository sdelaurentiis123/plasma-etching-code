#!/usr/bin/env python3
"""Ingest and digitize Freddie's 17 Aug 2026 Columbia SEM set (Oxford NPG80
TiO2 pillar dose series).

- copies the original Zeiss TIFFs into data/experimental/.../sem_0817/originals
  with sha256 receipts;
- reads the exact Zeiss SmartSEM metadata (Image Pixel Size, Mag, stage tilt,
  WD, EHT, timestamp) from TIFF tag 34118 - no scale-bar eyeballing;
- segments pillars on plan-view frames (InLens: pillars bright on dark floor)
  and reports per-frame CD / pitch statistics in nm, plus an elongation-based
  fallen-pillar fraction; tilted frames are catalogued, not measured;
- writes manifest.json, measurements.json, and contact sheets.

No model quantity enters this script.  It is the answer key, digitized.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw
from scipy import ndimage

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "experimental" / "zhu_2026_tio2_npg80" / "sem_0817"
FOOTER_ROW = 683          # Zeiss info bar starts here on 1024x768 frames


def zeiss_meta(im: Image.Image) -> dict:
    raw = im.tag_v2.get(34118, b"")
    text = raw.decode("latin-1", "ignore") if isinstance(raw, bytes) else str(raw)
    def grab(key, pattern=r"([-+0-9.]+)\s*([A-Za-zµ°]*)"):
        m = re.search(re.escape(key) + r"\s*=\s*" + pattern, text)
        return (float(m.group(1)), m.group(2)) if m else (None, None)
    px, unit = grab("Image Pixel Size")
    scale = {"nm": 1.0, "µm": 1000.0, "um": 1000.0, "pm": 1e-3}.get(unit, 1.0)
    mag = re.search(r"Mag\s*=\s*([0-9.]+)\s*(K?)\s*X", text)
    tilt = re.search(r"Stage at T\s*=\s*([-0-9.]+)", text)
    wd = re.search(r"WD\s*=\s*([0-9.]+)\s*mm", text)
    eht = re.search(r"EHT\s*=\s*([0-9.]+)\s*kV", text)
    date = re.search(r"Date\s*:(\d{1,2} \w{3} \d{4})", text)
    time = re.search(r"Time\s*:(\d{2}:\d{2}:\d{2})", text)
    return {
        "pixel_size_nm": None if px is None else px * scale,
        "magnification": None if not mag else float(mag.group(1)) * (1000.0 if mag.group(2) else 1.0),
        "stage_tilt_deg": None if not tilt else float(tilt.group(1)),
        "working_distance_mm": None if not wd else float(wd.group(1)),
        "eht_kv": None if not eht else float(eht.group(1)),
        "date": date.group(1) if date else None,
        "time": time.group(1) if time else None,
    }


def segment_pillars(gray: np.ndarray, pixel_nm: float) -> dict:
    """Plan-view pillar statistics from an InLens frame (footer removed)."""
    img = gray[:FOOTER_ROW].astype(float)
    img = ndimage.gaussian_filter(img, 1.2)
    hist, edges = np.histogram(img, bins=256)
    centers = 0.5 * (edges[:-1] + edges[1:])
    w1 = np.cumsum(hist); w2 = w1[-1] - w1
    m1 = np.cumsum(hist * centers) / np.maximum(w1, 1)
    m2 = (np.cumsum((hist * centers)[::-1])[::-1]) / np.maximum(w2, 1)
    var = w1[:-1] * w2[1:] * (m1[:-1] - m2[1:]) ** 2
    thr = centers[int(np.argmax(var))]
    mask = img > thr
    mask = ndimage.binary_opening(mask, iterations=2)
    labels, n = ndimage.label(mask)
    if n == 0:
        return {"pillar_count": 0}
    objs = ndimage.find_objects(labels)
    rows = []
    h, w = mask.shape
    for i, sl in enumerate(objs, start=1):
        if sl is None:
            continue
        y0, y1, x0, x1 = sl[0].start, sl[0].stop, sl[1].start, sl[1].stop
        if x0 <= 1 or y0 <= 1 or x1 >= w - 1 or y1 >= h - 1:
            continue
        area = int(np.sum(labels[sl] == i))
        if area < 40:
            continue
        cy, cx = ndimage.center_of_mass(labels[sl] == i)
        rows.append({"area_px": area, "bbox_w_px": x1 - x0, "bbox_h_px": y1 - y0,
                     "cx": x0 + cx, "cy": y0 + cy})
    if len(rows) < 4:
        return {"pillar_count": len(rows)}
    areas = np.array([r["area_px"] for r in rows], float)
    bw = np.array([r["bbox_w_px"] for r in rows], float)
    bh = np.array([r["bbox_h_px"] for r in rows], float)
    aspect = np.maximum(bw, bh) / np.maximum(np.minimum(bw, bh), 1)
    elongated = aspect > 1.8
    keep = ~elongated
    eq_width_nm = np.sqrt(areas[keep]) * pixel_nm
    bbox_w_nm = bw[keep] * pixel_nm
    pts = np.array([[r["cx"], r["cy"]] for r in rows])[keep]
    pitch_nm = None
    if len(pts) >= 3:
        d = np.sqrt(((pts[:, None, :] - pts[None, :, :]) ** 2).sum(-1))
        np.fill_diagonal(d, np.inf)
        pitch_nm = float(np.median(np.sort(d, axis=1)[:, 0]) * pixel_nm)
    return {
        "pillar_count": int(len(rows)),
        "compact_count": int(keep.sum()),
        "elongated_fraction": float(elongated.mean()),
        "equivalent_square_width_nm": {
            "median": float(np.median(eq_width_nm)),
            "p10": float(np.percentile(eq_width_nm, 10)),
            "p90": float(np.percentile(eq_width_nm, 90)),
        } if keep.any() else None,
        "bbox_width_nm_median": float(np.median(bbox_w_nm)) if keep.any() else None,
        "nearest_neighbour_pitch_nm": pitch_nm,
        "otsu_threshold": float(thr),
    }


def contact_sheets(records, src_dir: Path, out_dir: Path, per_sheet=20, cols=5):
    thumbs = []
    for r in records:
        im = Image.open(src_dir / r["file"]).convert("L")
        im = im.crop((0, 0, 1024, FOOTER_ROW)).resize((256, 171))
        canvas = Image.new("L", (256, 190), 0)
        canvas.paste(im, (0, 0))
        d = ImageDraw.Draw(canvas)
        tilt = r['meta']['stage_tilt_deg']
        tag = f"{r['file'][:-4]}  {r['meta']['pixel_size_nm']:.1f}nm/px" + (f"  T{tilt:.0f}" if tilt else "")
        d.text((3, 173), tag, fill=255)
        thumbs.append(canvas)
    sheets = []
    for s in range(0, len(thumbs), per_sheet):
        chunk = thumbs[s:s + per_sheet]
        rows = (len(chunk) + cols - 1) // cols
        sheet = Image.new("L", (cols * 256, rows * 190), 0)
        for k, t in enumerate(chunk):
            sheet.paste(t, ((k % cols) * 256, (k // cols) * 190))
        p = out_dir / f"contact_sheet_{s//per_sheet+1}.png"
        sheet.save(p)
        sheets.append(p.name)
    return sheets


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", default=str(Path.home() / "Downloads" / "0817"))
    ap.add_argument("--copy-originals", action="store_true")
    args = ap.parse_args()
    src = Path(args.source)
    OUT.mkdir(parents=True, exist_ok=True)
    orig = OUT / "originals"
    records = []
    for p in sorted(src.glob("*.tif")):
        sha = hashlib.sha256(p.read_bytes()).hexdigest()
        if args.copy_originals:
            orig.mkdir(exist_ok=True)
            if not (orig / p.name).exists():
                shutil.copy2(p, orig / p.name)
        im = Image.open(p)
        meta = zeiss_meta(im)
        gray = np.asarray(im.convert("L"))
        m = re.match(r"dose(\d+)(?:_(\d+)d)?_(\d+)", p.stem)
        rec = {"file": p.name, "sha256": sha,
               "dose_index": int(m.group(1)) if m else None,
               "series_tilt_label_deg": int(m.group(2)) if m and m.group(2) else 0,
               "frame_index": int(m.group(3)) if m else None, "meta": meta}
        tilted = (meta["stage_tilt_deg"] or 0) > 12 or rec["series_tilt_label_deg"] > 0
        if not tilted and meta["pixel_size_nm"] and meta["pixel_size_nm"] <= 12.0:
            rec["measurement"] = segment_pillars(gray, meta["pixel_size_nm"])
        else:
            rec["measurement"] = {"skipped": "tilted or low magnification"}
        records.append(rec)
    sheets = contact_sheets(records, src, OUT)
    (OUT / "manifest.json").write_text(json.dumps({
        "source_folder": str(src), "frame_count": len(records),
        "instrument": "Zeiss SEM, Columbia University (from frame metadata)",
        "acquired": sorted({r["meta"]["date"] for r in records if r["meta"]["date"]}),
        "contact_sheets": sheets,
        "frames": [{k: v for k, v in r.items() if k != "measurement"} for r in records],
    }, indent=1))
    (OUT / "measurements.json").write_text(json.dumps(
        [{"file": r["file"], "dose_index": r["dose_index"],
          "pixel_size_nm": r["meta"]["pixel_size_nm"],
          "stage_tilt_deg": r["meta"]["stage_tilt_deg"],
          **r["measurement"]} for r in records], indent=1))
    print(f"{len(records)} frames -> {OUT}")


if __name__ == "__main__":
    main()
