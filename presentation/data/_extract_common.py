"""Shared helpers of the presentation/data extractors (run in the LCG env, `source fitting/setup.sh`).

Every extractor freezes channel numbers into presentation/data/<name>.json with a provenance
block and a list of assertions against the channel anchors (z-mumu/handoff.md,
output/v2/RESULTS_v2.md, docs/). Mismatches are collected, printed and raised -- never
adjusted. Scenes read the JSON with `style.bnd_style.load_data(name)` and never recompute.
"""

from __future__ import annotations

import datetime as _dt
import json
import subprocess
import sys
from pathlib import Path

import numpy as np

DATA_DIR = Path(__file__).resolve().parent                      # presentation/data
PRESENTATION = DATA_DIR.parent
REPO = PRESENTATION.parent                                      # BND-school
ZMUMU = REPO / "z-mumu"
V2 = ZMUMU / "output" / "v2"
FIT_RESULTS = ZMUMU / "fit" / "results"
WORK = PRESENTATION / "work" / "extract"                        # part files (git-ignored)

LUMI_PB = 16393.381            # normtag, CMS Open Data record 1059 (docs/CONVENTIONS.md)
DATASET = "CMS 2016 Open Data, SingleMuon Run2016G+H, NanoAODv9 (records 30530, 30563)"

# ---- Z -> tautau (z-tautau v3: DeepTau Tight on both legs, MC-subtracted fake factor) -----------------------
ZTAUTAU = REPO / "z-tautau"
TAUTAU_NTUPLES = Path("/data/atlas/users/sjankovy/BND-school-cache/ztautau/ntuples_v1")   # flat TTree 'ntuple', one file per sample
DATASET_TAUTAU = "CMS 2016 Open Data, Tau Run2016G+H, NanoAODv9 (records 30532, 30565)"
VERSION_TAUTAU = "v3 DeepTau Tight, MC-subtracted FF"


def add_zmumu_path() -> None:
    """Make `import zmumu` and the v2 scripts importable (read-only use of the channel code)."""
    for p in (str(ZMUMU), str(ZMUMU / "scripts")):
        if p not in sys.path:
            sys.path.insert(0, p)


def add_ztautau_path() -> None:
    """Make `import ztautau` (z-tautau/ztautau: config, mass, fakes, analysis, corrections) importable, read-only."""
    for p in (str(ZTAUTAU), str(ZTAUTAU / "scripts")):
        if p not in sys.path:
            sys.path.insert(0, p)


def git_commit() -> dict:
    """Commit of the repository the numbers were frozen from (+ dirty flag)."""
    def run(*args):
        return subprocess.run(["git", "-C", str(REPO), *args], capture_output=True, text=True, check=False).stdout.strip()
    full = run("rev-parse", "HEAD")
    dirty = bool(run("status", "--porcelain", "--untracked-files=no"))
    return {"commit": full, "short": full[:8], "dirty_worktree": dirty}


def provenance(extractor: str, sources, **extra) -> dict:
    """Provenance block: extractor, date, git commit, luminosity, dataset, source files (+ anything else)."""
    srcs, links = [], {}
    for s in sources:
        if not str(s).startswith(str(REPO)):
            srcs.append(str(s))
            continue
        try:
            srcs.append(str(Path(s).resolve().relative_to(REPO)))
        except ValueError:                      # a git-ignored link inside the worktree that points to the main checkout
            rel = str(Path(s).relative_to(REPO))
            srcs.append(rel)
            links[rel] = str(Path(s).resolve())
    out = {"extractor": f"presentation/data/{extractor}", "generated": _dt.datetime.now().isoformat(timespec="seconds"),
           "git": git_commit(), "lumi_pb": LUMI_PB, "dataset": DATASET, "sources": srcs, **extra}
    if links:
        out["source_links"] = links
    return out


class Checker:
    """Collect assertions against anchors; `finish()` prints all of them and raises on any failure."""

    def __init__(self, label: str):
        self.label = label
        self.rows: list[dict] = []

    def check(self, name, got, want, tol=0.0, source="", rel=False, soft=False):
        """Scalar or array comparison. `tol` absolute (or relative with `rel=True`); arrays use allclose.
        `soft=True`: a documentation-level comparison -- a mismatch is printed as SOFT-FAIL and stored, but does
        not make `finish()` raise (used for numbers quoted in docs from an earlier pass / rounded tables)."""
        g, w = np.asarray(got, dtype=float), np.asarray(want, dtype=float)
        if g.shape != w.shape:
            ok, maxdiff = False, float("nan")
        else:
            diff = np.abs(g - w)
            scale = np.maximum(np.abs(w), 1e-300) if rel else 1.0
            ok = bool(np.all(diff <= tol * scale)) if g.size else True
            maxdiff = float(np.max(diff / scale)) if g.size else 0.0
        row = {"name": name, "ok": ok, "tol": float(tol), "rel": bool(rel), "max_diff": maxdiff, "source": source}
        if soft:
            row["soft"] = True
        if g.ndim == 0:
            row["got"], row["want"] = float(g), float(w)
        else:
            row["shape"] = list(g.shape)
        self.rows.append(row)
        return ok

    def check_true(self, name, cond, source="", detail=""):
        row = {"name": name, "ok": bool(cond), "source": source}
        if detail:
            row["detail"] = detail
        self.rows.append(row)
        return bool(cond)

    def finish(self, log=print) -> list[dict]:
        n_bad = 0
        for r in self.rows:
            flag = "PASS" if r["ok"] else ("SOFT-FAIL" if r.get("soft") else "FAIL")
            n_bad += (not r["ok"]) and not r.get("soft")
            if "got" in r:
                extra = f"got {r['got']:.10g}  want {r['want']:.10g}  tol {r['tol']:g}{' rel' if r['rel'] else ''}"
            elif "shape" in r:
                extra = f"shape {r['shape']}  max|d| {r['max_diff']:.3g}  tol {r['tol']:g}{' rel' if r['rel'] else ''}"
            else:
                extra = r.get("detail", "")
            log(f"  [{flag}] {r['name']}: {extra}" + (f"   ({r['source']})" if r["source"] else ""))
        n_soft = sum(1 for r in self.rows if (not r["ok"]) and r.get("soft"))
        log(f"[{self.label}] {len(self.rows) - n_bad - n_soft}/{len(self.rows)} checks passed"
            + (f", {n_soft} soft check(s) reported (not raised)" if n_soft else ""))
        if n_bad:
            raise AssertionError(f"{self.label}: {n_bad} check(s) failed -- numbers are NOT adjusted; see the list above")
        return self.rows


def _default(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.bool_,)):
        return bool(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, Path):
        return str(o)
    raise TypeError(f"not JSON serialisable: {type(o)}")


def dump_json(path, obj, indent=1) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".json.tmp")
    with open(tmp, "w") as fh:
        json.dump(obj, fh, indent=indent, default=_default, allow_nan=False, separators=None if indent else (",", ":"))
        fh.write("\n")
    tmp.rename(path)
    print(f"[{path.name}] wrote {path} ({path.stat().st_size / 1024:.1f} kB)")
    return path


def load_json(path):
    with open(path) as fh:
        return json.load(fh)


def rebin(values, factor: int):
    """Sum `factor` adjacent bins on the last axis (same as z-mumu/scripts/v2_5_fit.py:rebin)."""
    v = np.asarray(values, dtype=float)
    return v if factor == 1 else v.reshape(*v.shape[:-1], -1, factor).sum(axis=-1)


def standard_args(parser, default_json: str):
    """`--json PATH` (never `--out`: reserved by zmumu/batch.py) and `--check-only`."""
    parser.add_argument("--json", type=Path, default=DATA_DIR / default_json, help=f"output JSON (default {default_json})")
    parser.add_argument("--check-only", action="store_true", help="re-run the assertions on the existing JSON, write nothing")
    return parser


def finalize(args, build, verify, label, indent=1):
    """Common main: build (or load with --check-only), verify (prints + raises), write with the checks embedded."""
    if args.check_only:
        d = load_json(args.json)
        print(f"[{label}] --check-only on {args.json}")
        verify(d, Checker(label)).finish()
        return d
    d = build()
    rows = verify(d, Checker(label)).finish()
    d["checks"] = rows
    dump_json(args.json, d, indent=indent)
    return d
