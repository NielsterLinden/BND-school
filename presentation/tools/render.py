#!/usr/bin/env python
"""Render one scene and (at quality h) deliver it as a numbered, versioned MP4.

    python tools/render.py <section> <name> <scene_file> <SceneClass> [--quality l|h]
                           [--version N] [--frames t1,t2,...] [--no-deliver]

section  : 1..7 (see SECTIONS) or the directory name (e.g. 3_zee)
name     : short snake_case clip name (the registry key). Numbers are assigned
           in creation order per section from clips/CLIPLIST.tsv and never reused.
quality  : l (480p15, iterate) | h (1080p60, deliver). Default l.
version  : delivered version; default = 1 + highest existing _vN of this clip.

Outputs
  work/<name>/<name>.mp4 + frames         always (the thing you inspect)
  clips/<S>_<section>/<S>-<NN>_<name>_v<K>.mp4   quality h only (the deliverable)
The user only wants MP4s in clips/. Never write PNGs there.
"""
from __future__ import annotations

import argparse
import csv
import re
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SECTIONS = {
    "1": "1_theory", "2": "2_cms_methods", "3": "3_zee",
    "4": "4_zmumu", "5": "5_ztautau", "6": "6_combination", "7": "7_outro",
}
REGISTRY = ROOT / "clips" / "CLIPLIST.tsv"
FIELDS = ["number", "section", "name", "scene_file", "scene_class"]


def section_dir(s: str) -> str:
    if s in SECTIONS:
        return SECTIONS[s]
    if s in SECTIONS.values():
        return s
    sys.exit(f"unknown section {s!r}; use one of {list(SECTIONS)} or {list(SECTIONS.values())}")


def load_registry() -> list[dict]:
    if not REGISTRY.exists():
        return []
    with REGISTRY.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def save_registry(rows: list[dict]) -> None:
    REGISTRY.parent.mkdir(parents=True, exist_ok=True)
    with REGISTRY.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS, delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def clip_number(sec: str, name: str, scene_file: str, scene_class: str) -> int:
    rows = load_registry()
    for r in rows:
        if r["section"] == sec and r["name"] == name:
            if (r["scene_file"], r["scene_class"]) != (scene_file, scene_class):
                r["scene_file"], r["scene_class"] = scene_file, scene_class
                save_registry(rows)
            return int(r["number"])
    n = 1 + max((int(r["number"]) for r in rows if r["section"] == sec), default=0)
    rows.append(dict(number=n, section=sec, name=name, scene_file=scene_file,
                     scene_class=scene_class))
    save_registry(rows)
    return n


def next_version(sec_dir: Path, stem: str) -> int:
    pat = re.compile(re.escape(stem) + r"_v(\d+)\.mp4$")
    vs = [int(m.group(1)) for p in sec_dir.glob(f"{stem}_v*.mp4") if (m := pat.search(p.name))]
    return 1 + max(vs, default=0)


def grab(mp4: Path, out: Path, t: str) -> None:
    if t == "end":
        # decode the last 0.5 s and keep overwriting one PNG: the file ends as the true last frame
        cmd = ["ffmpeg", "-y", "-loglevel", "error", "-sseof", "-0.5", "-i", str(mp4),
               "-update", "1", str(out)]
    else:
        cmd = ["ffmpeg", "-y", "-loglevel", "error", "-ss", t, "-i", str(mp4),
               "-frames:v", "1", "-update", "1", str(out)]
    subprocess.run(cmd, check=True)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("section"); ap.add_argument("name")
    ap.add_argument("scene_file"); ap.add_argument("scene_class")
    ap.add_argument("--quality", "-q", default="l", choices=["l", "m", "h"])
    ap.add_argument("--version", "-v", type=int, default=None)
    ap.add_argument("--frames", default="0,end", help="comma list of times in s, or 'end'")
    ap.add_argument("--no-deliver", action="store_true")
    a = ap.parse_args()

    sec = section_dir(a.section)
    s_no = sec.split("_")[0]
    work = ROOT / "work" / a.name
    work.mkdir(parents=True, exist_ok=True)
    scene_file = str(Path(a.scene_file).as_posix())

    cmd = ["manim", "render", "--quality", a.quality, "--media_dir", str(work / "media"),
           "-o", a.name, scene_file, a.scene_class]
    print(">>", " ".join(cmd), flush=True)
    log = subprocess.run(cmd, cwd=ROOT, text=True, capture_output=True)
    (work / "render.log").write_text(log.stdout + log.stderr)
    tail = (log.stdout + log.stderr).strip().splitlines()[-25:]
    print("\n".join(tail))
    if log.returncode != 0 or "Rendered" not in (log.stdout + log.stderr):
        sys.exit(f"render FAILED (see {work/'render.log'})")

    mp4s = sorted((work / "media" / "videos").rglob(f"{a.name}.mp4"), key=lambda p: p.stat().st_mtime)
    if not mp4s:
        sys.exit("no mp4 produced")
    src = mp4s[-1]
    out = work / f"{a.name}.mp4"
    shutil.copy2(src, out)
    for t in [t.strip() for t in a.frames.split(",") if t.strip()]:
        png = work / f"{a.name}_{'final' if t == 'end' else 't' + t.replace('.', 'p')}.png"
        grab(out, png, t)
        print("frame:", png)
    print("work clip:", out)

    if a.quality == "h" and not a.no_deliver:
        n = clip_number(sec, a.name, scene_file, a.scene_class)
        sec_dir = ROOT / "clips" / sec
        sec_dir.mkdir(parents=True, exist_ok=True)
        stem = f"{s_no}-{n:02d}_{a.name}"
        v = a.version or next_version(sec_dir, stem)
        dst = sec_dir / f"{stem}_v{v}.mp4"
        shutil.copy2(out, dst)
        print("DELIVERED:", dst.relative_to(ROOT))


if __name__ == "__main__":
    main()
