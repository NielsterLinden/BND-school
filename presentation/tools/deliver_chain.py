#!/usr/bin/env python
"""Render a chain of scenes and deliver every section as its own numbered clip.

    python tools/deliver_chain.py <section> <scene_file> <Scene1> [<Scene2> ...]
                                  [-q l|h] [--no-deliver] [--only name1,name2] [--seam-before MP4]

The scenes mark their clips with ``clip_open(scene, name)`` / ``clip_cut(scene, name)``
(style/bnd_style.py = manim ``next_section``). Every scene is rendered in the given order
with ``--save_sections``. Sections carrying the same name at consecutive positions of the
chain (e.g. a zoom-out at the end of one scene and the zoom-in at the start of the next) are
joined into one clip (ffmpeg concat, stream copy: no re-encode, exact frames).

Outputs per clip (always): work/<name>/<name>.mp4, <name>_t0.png, <name>_final.png.
Quality h without --no-deliver: numbered through clips/CLIPLIST.tsv in chain order and copied to
clips/<S>_<section>/<S>-<NN>_<name>_v<K>.mp4 (same rules as tools/render.py).
Seams: tools/framediff.py on every (final of clip k, t0 of clip k+1); ``--seam-before`` adds
the seam from an already delivered clip into the first clip of the chain.
"""
from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from render import ROOT, clip_number, grab, next_version, section_dir  # noqa: E402
from keepout import report as keepout_report  # noqa: E402

QDIR = {"l": "480p15", "m": "720p30", "h": "1080p60"}


def render_scene(scene_file: str, cls: str, quality: str, media: Path) -> list[tuple[str, Path]]:
    cmd = ["manim", "render", "--quality", quality, "--save_sections", "--media_dir", str(media),
           "-o", cls, scene_file, cls]
    print(">>", " ".join(cmd), flush=True)
    log = subprocess.run(cmd, cwd=ROOT, text=True, capture_output=True)
    logdir = ROOT / "work" / "chains"
    logdir.mkdir(parents=True, exist_ok=True)
    (logdir / f"{cls}_{quality}.log").write_text(log.stdout + log.stderr)
    if log.returncode != 0 or "Rendered" not in (log.stdout + log.stderr):
        print("\n".join((log.stdout + log.stderr).strip().splitlines()[-30:]))
        sys.exit(f"render FAILED: {cls} (see {logdir / f'{cls}_{quality}.log'})")
    index = sorted((media / "videos").rglob(f"{QDIR[quality]}/sections/{cls}.json"),
                   key=lambda p: p.stat().st_mtime)
    if not index:
        sys.exit(f"no sections index for {cls}")
    entries = json.loads(index[-1].read_text())
    out = [(e["name"], index[-1].parent / e["video"]) for e in entries]
    print(f"   {cls}: " + ", ".join(n for n, _ in out), flush=True)
    return out


def duration(mp4: Path) -> float:
    r = subprocess.run(["ffprobe", "-v", "error", "-show_entries", "format=duration", "-of",
                        "default=nw=1:nk=1", str(mp4)], text=True, capture_output=True)
    return float(r.stdout.strip() or 0.0)


def join(parts: list[Path], out: Path) -> None:
    if len(parts) == 1:
        shutil.copy2(parts[0], out)
        return
    lst = out.with_suffix(".txt")
    lst.write_text("".join(f"file '{p}'\n" for p in parts))
    subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-f", "concat", "-safe", "0", "-i", str(lst),
                    "-c", "copy", str(out)], check=True)
    lst.unlink()


def framediff(a: Path, b: Path) -> str:
    r = subprocess.run([sys.executable, str(ROOT / "tools" / "framediff.py"), str(a), str(b)],
                       text=True, capture_output=True)
    return (r.stdout + r.stderr).strip()


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("section")
    ap.add_argument("scene_file")
    ap.add_argument("scenes", nargs="+")
    ap.add_argument("--quality", "-q", default="l", choices=["l", "m", "h"])
    ap.add_argument("--no-deliver", action="store_true")
    ap.add_argument("--only", default="", help="comma list of clip names to deliver (default: all)")
    ap.add_argument("--seam-before", default=None, help="delivered MP4 whose last frame opens the chain")
    a = ap.parse_args()

    sec = section_dir(a.section)
    s_no = sec.split("_")[0]
    scene_file = str(Path(a.scene_file).as_posix())
    media = ROOT / "work" / "chains" / Path(scene_file).stem / "media"

    sections: list[tuple[str, str, Path]] = []          # (clip name, scene class, section video)
    for cls in a.scenes:
        sections += [(n, cls, p) for n, p in render_scene(scene_file, cls, a.quality, media)]

    clips: list[dict] = []                               # contiguous equal names -> one clip
    for name, cls, p in sections:
        if clips and clips[-1]["name"] == name:
            clips[-1]["parts"].append(p)
            if cls not in clips[-1]["classes"]:
                clips[-1]["classes"].append(cls)
        else:
            if any(c["name"] == name for c in clips):
                sys.exit(f"section name {name!r} used twice, not contiguously")
            clips.append({"name": name, "classes": [cls], "parts": [p]})

    only = {s for s in a.only.split(",") if s}
    rows = []
    for c in clips:
        work = ROOT / "work" / c["name"]
        work.mkdir(parents=True, exist_ok=True)
        out = work / f"{c['name']}.mp4"
        join(c["parts"], out)
        c["mp4"], c["t0"], c["final"] = out, work / f"{c['name']}_t0.png", work / f"{c['name']}_final.png"
        grab(out, c["t0"], "0")
        grab(out, c["final"], "end")
        c["dur"] = duration(out)
        c["keepout"] = keepout_report(out, step=5 if a.quality == "l" else 15)[1]
        dst = ""
        if a.quality == "h" and not a.no_deliver and (not only or c["name"] in only):
            n = clip_number(sec, c["name"], scene_file, "+".join(c["classes"]))
            sec_dir = ROOT / "clips" / sec
            sec_dir.mkdir(parents=True, exist_ok=True)
            stem = f"{s_no}-{n:02d}_{c['name']}"
            dst_p = sec_dir / f"{stem}_v{next_version(sec_dir, stem)}.mp4"
            shutil.copy2(out, dst_p)
            dst = str(dst_p.relative_to(ROOT))
        rows.append((c["name"], "+".join(c["classes"]), c["dur"], dst))

    print("\nclips:")
    for name, cls, dur, dst in rows:
        print(f"  {name:28s} {dur:5.2f} s  {cls:28s} {dst}")
    print(f"  total {sum(r[2] for r in rows):.1f} s in {len(rows)} clips")

    print("\nkeep-out (title band y > 2.7, top-left 3 x 9 cm):")
    for c in clips:
        print("  " + c["keepout"])

    print("\nseams (final of k vs t0 of k+1):")
    if a.seam_before:
        prev = ROOT / "work" / "chains" / "seam_before_final.png"
        grab(Path(a.seam_before), prev, "end")
        print(f"  {Path(a.seam_before).name} -> {clips[0]['name']}: {framediff(prev, clips[0]['t0'])}")
    for k in range(len(clips) - 1):
        print(f"  {clips[k]['name']} -> {clips[k + 1]['name']}: {framediff(clips[k]['final'], clips[k + 1]['t0'])}")


if __name__ == "__main__":
    main()
