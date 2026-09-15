#!/usr/bin/env python
"""Pixel diff between two frames (the joint check for chained clips).

    python tools/framediff.py <last_frame_of_clip_N.png> <first_frame_of_clip_N+1.png>

Reports mean|d|, max|d|, fraction of pixels over threshold, and a blob count
(pixels whose whole 5x5 neighbourhood is over threshold). Calibration at 1080p:
mean 0.2-0.6, max up to ~60 on hairlines, fraction < 0.1 %, blob 0.
FAIL at mean > 1.0, fraction > 0.5 %, or blob > 50.
"""
import sys

import numpy as np
from PIL import Image

a = np.asarray(Image.open(sys.argv[1]).convert("RGB"), dtype=float)
b = np.asarray(Image.open(sys.argv[2]).convert("RGB"), dtype=float)
if a.shape != b.shape:
    sys.exit(f"shape mismatch {a.shape} vs {b.shape}")
d = np.abs(a - b).max(axis=2)
thr = 24.0
over = d > thr
k = 5
er = over.copy()
for _ in range(k // 2):
    e = er.copy()
    e[1:] &= er[:-1]; e[:-1] &= er[1:]; e[:, 1:] &= er[:, :-1]; e[:, :-1] &= er[:, 1:]
    er = e
mean, mx, frac, blob = d.mean(), d.max(), over.mean() * 100, int(er.sum())
ok = mean <= 1.0 and frac <= 0.5 and blob <= 50
print(f"mean|d|={mean:.3f}  max|d|={mx:.0f}  over{thr:.0f}={frac:.3f}%  blob={blob}  -> {'OK' if ok else 'FAIL'}")
sys.exit(0 if ok else 1)
