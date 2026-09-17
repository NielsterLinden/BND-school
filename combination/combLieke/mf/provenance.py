"""What a combination was built from: the input files of every channel with their checksum and date, and the last
commit of the channel directory. `run.py prepare` writes it to work/status.json, `results` copies it to result.json,
so that "the final inputs" is something one can check (the ROOT inputs are not in git)."""

from __future__ import annotations

import hashlib
import re
import subprocess
import time
from pathlib import Path

from .paths import REPO, repo_path


def _git(*args: str) -> str | None:
    r = subprocess.run(["git", "-C", str(REPO), *args], capture_output=True, text=True)
    return r.stdout.strip() or None if r.returncode == 0 else None


def _file(path: Path) -> dict:
    if not path.exists():
        return {"file": str(path.relative_to(REPO)), "missing": True}
    return {"file": str(path.relative_to(REPO)), "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "modified": time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(path.stat().st_mtime)), "bytes": path.stat().st_size}


def channel_files(spec: dict) -> list[Path]:
    """The channel's TRExFitter config, the histogram files it reads and the result the channel published."""
    config = repo_path(spec["config"])
    files = [config]
    if "inputs" in spec:                                 # ee: histograms and config come out of one tarball
        files.append(repo_path(spec["inputs"]["tarball"]))
    if "histo_path" in spec:
        names = dict.fromkeys(re.findall(r'^\s*HistoFile:\s*"?([^"\s]+)"?\s*$', config.read_text(), flags=re.M))
        files += [repo_path(spec["histo_path"]) / f"{n}.root" for n in names]
    if isinstance(spec.get("signal"), dict):
        files.append(repo_path(spec["signal"]["meta"]))
    source = (spec.get("published") or {}).get("source", "")
    if source.endswith(".json"):
        files.append(repo_path(source))
    return list(dict.fromkeys(files))


def record(m: dict) -> dict:
    out = {"repository_head": _git("log", "-1", "--format=%h %cI"), "channels": {}}
    for key, spec in m["channels"].items():
        top = Path(spec["config"]).parts[0]
        out["channels"][key] = {"directory": top, "last_commit": _git("log", "-1", "--format=%h %cI %s", "--", top),
                                "uncommitted": bool(_git("status", "--porcelain", "--", top)),
                                "files": [_file(p) for p in channel_files(spec)]}
    return out
