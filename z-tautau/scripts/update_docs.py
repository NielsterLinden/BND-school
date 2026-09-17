#!/usr/bin/env python
"""Inject the current result into the documents that quote it, so no number is ever typed by hand.

    python scripts/update_docs.py

`scripts/step6_report.py` writes `output/result_block.md` from the fit results. This script copies it into
every file that carries the markers

    <!-- RESULT:BEGIN --> ... <!-- RESULT:END -->

(README.md, handoff.md, docs/00-overview.md, docs/10-v4-plan.md). Run it after step 6; `run_all.py` does.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from ztautau import config  # noqa: E402

TARGETS = ["README.md", "handoff.md", "docs/00-overview.md", "docs/10-v4-plan.md"]
RX = re.compile(r"<!-- RESULT:BEGIN -->.*?<!-- RESULT:END -->", re.S)


def main():
    block = (config.OUTPUT_DIR_V4 / "result_block.md").read_text().strip()
    new = f"<!-- RESULT:BEGIN -->\n{block}\n<!-- RESULT:END -->"
    for name in TARGETS:
        p = config.CHANNEL_DIR / name
        if not p.exists():
            continue
        t = p.read_text()
        if not RX.search(t):
            print(f"[skip] {name}: no RESULT markers")
            continue
        p.write_text(RX.sub(lambda _: new, t))
        print(f"[ok]   {name}")


if __name__ == "__main__":
    main()
