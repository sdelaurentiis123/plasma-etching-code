#!/usr/bin/env python3
"""Assemble showcase/index.html: inline the data JSON and media into the
template so the page is a single self-contained file."""
from __future__ import annotations

import base64
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent

MEDIA = {
    "__FLYBY_MP4__": (HERE / "media" / "flyby.mp4", "video/mp4"),
    "__HERO3D_MP4__": (ROOT.parent / "bench" / "out" / "hero3d.mp4",
                       "video/mp4"),
    "__ETCH_HOLE_GIF__": (ROOT / "viz" / "etch_hole.gif", "image/gif"),
    "__ETCH_TRENCH_GIF__": (ROOT / "viz" / "etch_trench.gif", "image/gif"),
}


def main() -> None:
    html = (HERE / "template.html").read_text()
    data = (HERE / "data" / "showcase_data.json").read_text()
    html = html.replace("__DATA_JSON__", data)
    for token, (path, mime) in MEDIA.items():
        b64 = base64.b64encode(path.read_bytes()).decode()
        html = html.replace(token, f"data:{mime};base64,{b64}")
    out = HERE / "index.html"
    out.write_text(html)
    print(f"wrote {out} ({out.stat().st_size/1e6:.2f} MB)")


if __name__ == "__main__":
    main()
