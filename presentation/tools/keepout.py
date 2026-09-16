#!/usr/bin/env python
"""Scan rendered MP4s for ink in the user's PowerPoint overlays: the title band (y > 2.70,
top 175 px at 1080p) and the top-left block (3 cm x 9 cm of the slide: 170 x 510 px at 1080p).

    python tools/keepout.py clips/*/*.mp4            # every 15th frame of each clip
    python tools/keepout.py work/pipe_d1_tag_probe/pipe_d1_tag_probe.mp4 --step 1

A clip passes when at most 0.1 % of either region is darker than the background (H.264 noise
on white stays far below that). ``EXEMPT`` lists the clips whose full-frame endings cover the
overlays by design (they covered the title band before the block existed).
Anchors: style/bnd_style.py TITLE_BAND_Y, CORNER_X_MAX, CORNER_Y_MIN.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np

TITLE_PX, CORNER_W_PX, CORNER_H_PX = 175, 170, 510      # at 1920 x 1080
INK, TOL = 225, 0.001
EXEMPT = {
    "ee_detector": "full-frame zoom ending",
    "tautau_detector": "full-frame zoom ending",
    "tautau_a1_event": "opens on the tautau_detector zoom (first 1.4 s)",
    "mic_drop": "outro rays",
}


def clip_name(mp4: Path) -> str:
    m = re.match(r"\d-\d\d_(.+)_v\d+$", mp4.stem)
    return m.group(1) if m else mp4.stem


def scan(mp4: Path, step: int = 15) -> dict:
    import av
    worst = {"title": (0.0, 0.0), "corner": (0.0, 0.0)}
    with av.open(str(mp4)) as c:
        s = c.streams.video[0]
        fps = float(s.average_rate)
        for i, fr in enumerate(c.decode(s)):
            if i % step:
                continue
            im = fr.to_ndarray(format="gray")
            k = im.shape[1] / 1920
            regions = {"title": im[:round(TITLE_PX * k), :],
                       "corner": im[:round(CORNER_H_PX * k), :round(CORNER_W_PX * k)]}
            for key, r in regions.items():
                f = float((r < INK).mean())
                if f > worst[key][0]:
                    worst[key] = (f, i / fps)
    return worst


def report(mp4: Path, step: int = 15) -> tuple[bool, str]:
    w = scan(mp4, step)
    name = clip_name(mp4)
    bad = [k for k, (f, _) in w.items() if f > TOL]
    txt = "  ".join(f"{k} {f * 100:5.2f}% @ {t:4.1f}s" for k, (f, t) in w.items())
    if not bad:
        return True, f"ok      {name:28s} {txt}"
    if name in EXEMPT:
        return True, f"exempt  {name:28s} {txt}  ({EXEMPT[name]})"
    return False, f"WARN    {name:28s} {txt}  <- keep-out: {', '.join(bad)}"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("mp4", nargs="+")
    ap.add_argument("--step", type=int, default=15, help="scan every N-th frame")
    a = ap.parse_args()
    fails = 0
    for p in a.mp4:
        ok, line = report(Path(p), a.step)
        fails += not ok
        print(line, flush=True)
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
