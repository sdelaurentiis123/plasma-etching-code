#!/usr/bin/env python3
"""Inline the audited Arun explorer payload into the HTML fragment."""
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent


def main():
    template = (HERE / "template.fragment.html").read_text(encoding="utf-8")
    data = (ROOT / "results" / "etch_explorer_data.json").read_text(
        encoding="utf-8").strip()
    output = template.replace("__ARUN_EXPLORER_DATA__", data)
    if "__ARUN_EXPLORER_DATA__" in output:
        raise RuntimeError("explorer data token was not replaced")
    target = HERE / "arun_etch_explorer.html"
    target.write_text(output, encoding="utf-8")
    print(target.relative_to(ROOT))


if __name__ == "__main__":
    main()
