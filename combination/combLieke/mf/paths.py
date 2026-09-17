"""Locations. Everything the combination writes is under combination/combLieke/."""

from __future__ import annotations

import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]          # combination/combLieke
REPO = HERE.parents[1]                              # repository root
CONFIG = HERE / "config"
WORK = HERE / "work"                                # TRExFitter jobs, logs (git-ignored)
OUTPUT = HERE / "output"
PLOTS = OUTPUT / "plots"

if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))


def manifest() -> dict:
    return json.loads((CONFIG / "channels.json").read_text())


def references() -> dict:
    return json.loads((CONFIG / "references.json").read_text())


def repo_path(p: str | Path) -> Path:
    p = Path(p)
    return p if p.is_absolute() else REPO / p
