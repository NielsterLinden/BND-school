#!/usr/bin/env python
"""Reserved-area check on rendered clips (06-chapter-anchors.md A1).

    python tools/zonecheck.py <clip.mp4 | directory> [...] [--fps 4] [--save work/zonecheck]

Nothing may be drawn where the user's PowerPoint overlays sit:
  * the title band       y > TITLE_BAND_Y (2.7)
  * the chapter identifier  x < CHAPTER_ID_X (-5.85) and y > CHAPTER_ID_Y (0.22)
    (the top-left 3 cm x 9 cm of the 33.867 x 19.05 cm slide)
(constants in style/bnd_style.py). Every clip is sampled at --fps plus its last frame, at 480p;
a pixel counts as ink when any channel is darker than --ink (white background = 255). A clip
FAILs when a sampled frame has more than --min-px ink pixels inside a zone; the worst frame is
saved to --save with the zones outlined. Directories are searched for *.mp4 (00_archive skipped).
Exit code 1 when any clip fails.
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from style.bnd_style import CHAPTER_ID_X, CHAPTER_ID_Y, TITLE_BAND_Y  # noqa: E402

W, H = 854, 480
FRAME_W, FRAME_H = 8.0 * 16 / 9, 8.0


def px(x: float, y: float) -> tuple[int, int]:
    return int(round((x + FRAME_W / 2) / FRAME_W * W)), int(round((FRAME_H / 2 - y) / FRAME_H * H))


EDGE = 2      # px at 480p (0.033 units) left out along the zone boundary: the slice's top edge sits exactly at y = 2.70
ZONES = {
    "title band": (0, 0, W, px(0.0, TITLE_BAND_Y)[1] - EDGE),
    "chapter id": (0, 0, px(CHAPTER_ID_X, 0.0)[0] - EDGE, px(0.0, CHAPTER_ID_Y)[1] - EDGE),
}


def frames(mp4: Path, fps: float) -> tuple[np.ndarray, np.ndarray]:
    """(times, frames[n, H, W, 3]) at `fps`, plus the last frame of the clip."""
    def run(args):
        r = subprocess.run(["ffmpeg", "-loglevel", "error", *args, "-vf", f"scale={W}:{H}", "-f", "rawvideo",
                            "-pix_fmt", "rgb24", "-"], capture_output=True, check=True)
        return np.frombuffer(r.stdout, dtype=np.uint8).reshape(-1, H, W, 3)
    body = run(["-i", str(mp4), "-r", str(fps)])
    tail = run(["-sseof", "-0.5", "-i", str(mp4)])[-1:]
    t = np.append(np.arange(len(body)) / fps, np.nan)
    return t, np.concatenate([body, tail])


def check(mp4: Path, args) -> bool:
    t, fr = frames(mp4, args.fps)
    ink = fr.min(axis=3) < args.ink
    worst, report = (0, None, None), []
    for zone, (x0, y0, x1, y1) in ZONES.items():
        n = ink[:, y0:y1, x0:x1].sum(axis=(1, 2))
        bad = np.flatnonzero(n > args.min_px)
        if bad.size:
            i = int(bad[np.argmax(n[bad])])
            ts = "end" if np.isnan(t[i]) else f"{t[i]:.2f}s"
            first = "end" if np.isnan(t[bad[0]]) else f"{t[bad[0]]:.2f}s"
            report.append(f"{zone}: {bad.size}/{len(n)} frames, from {first}, worst {int(n[i])} px at {ts}")
            if n[i] > worst[0]:
                worst = (int(n[i]), i, zone)
    ok = not report
    print(f"{'OK  ' if ok else 'FAIL'} {mp4.relative_to(ROOT) if mp4.is_relative_to(ROOT) else mp4}"
          + ("" if ok else "\n       " + "\n       ".join(report)))
    if not ok and args.save:
        out = Path(args.save)
        out.mkdir(parents=True, exist_ok=True)
        img = Image.fromarray(fr[worst[1]])
        d = ImageDraw.Draw(img)
        for x0, y0, x1, y1 in ZONES.values():
            d.rectangle([x0, y0, x1 - 1, y1 - 1], outline=(255, 0, 0), width=2)
        img.save(out / f"{mp4.stem}_zone.png")
    return ok


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("paths", nargs="+")
    ap.add_argument("--fps", type=float, default=4.0)
    ap.add_argument("--ink", type=int, default=235, help="a channel below this is ink (default 235)")
    ap.add_argument("--min-px", type=int, default=12, help="ink pixels tolerated per zone (H.264 noise)")
    ap.add_argument("--save", default=str(ROOT / "work" / "zonecheck"))
    args = ap.parse_args()
    clips = []
    for p in map(Path, args.paths):
        clips += sorted(q for q in p.rglob("*.mp4") if "00_archive" not in q.parts) if p.is_dir() else [p]
    n_bad = sum(not check(c.resolve(), args) for c in clips)
    print(f"[zonecheck] {len(clips) - n_bad}/{len(clips)} clips clear of the reserved areas")
    sys.exit(1 if n_bad else 0)


if __name__ == "__main__":
    main()
