"""Run trex-fitter and parse its text outputs.

    from fitting import run_trex
    run_trex.run("zmumu.config", "hwf", cwd="z-mumu/fit", log="results/zmumu/logs/hwf.log")
    res = run_trex.summarise("z-mumu/fit/results/zmumu", "zmumu", poi="mu_Z")

File formats (TRExFitter v1.8.0, verified on test/reference/*):
  Fits/<job>.txt                      NUISANCE_PARAMETERS block: "name  value +up -down"
  Fits/<job>_errDecomp_<poi>.txt      "KEY<TAB>sym<TAB>up<TAB>down" (TOT_ERROR, STAT_ERROR, ...)
  Fits/<job>_group_errDecomp_<poi>.txt / Fits/GroupedImpact_<poi>.txt
                                      "<Category>    <unc>  ( +<up>, -<down> )"
  Fits/NPRanking_<poi>.txt            "name pull +up -down dmu_up_post dmu_down_post dmu_up_pre dmu_down_pre"
  Tables/Table_postfit.yaml           yields per region/sample
The goodness-of-fit probability is only printed to the log of the `f` step.
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
from pathlib import Path


def trex_binary() -> str:
    exe = shutil.which("trex-fitter")
    if exe is None:
        home = os.environ.get("TREXFITTER_HOME")
        if home and (Path(home) / "build/bin/trex-fitter").exists():
            return str(Path(home) / "build/bin/trex-fitter")
        raise RuntimeError("trex-fitter not on PATH: source setup.sh and run fitting/build_trexfitter.sh")
    return exe


def run(config, actions, options: str | None = None, cwd=None, log=None, check=True) -> int:
    """Run `trex-fitter <actions> <config> [<options>]`, teeing stdout+stderr to `log`."""
    cmd = [trex_binary(), actions, str(config)]
    if options:
        cmd.append(options)
    cwd = str(cwd) if cwd else None
    proc = subprocess.Popen(cmd, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    fh = None
    if log:
        Path(cwd or ".", log).parent.mkdir(parents=True, exist_ok=True)
        fh = open(Path(cwd or ".", log), "w")
    for line in proc.stdout:
        if fh:
            fh.write(line)
        if re.search(r"ERROR|Error in|FATAL|GOODNESS|probability|fit status|Minimized", line):
            print("  [trex] " + line.rstrip())
    proc.wait()
    if fh:
        fh.close()
    if check and proc.returncode != 0:
        raise RuntimeError(f"trex-fitter {actions} failed with code {proc.returncode} (see {log})")
    return proc.returncode


# ------------------------------------------------------------------ parsers
def parse_fit_txt(path) -> dict:
    """{'nps': {name: (value, up, down)}, 'corr': {'names': [...], 'matrix': [[...]]}, 'nll': float}"""
    text = Path(path).read_text().splitlines()
    out = {"nps": {}, "corr": None, "nll": None}
    section = None
    corr_rows, corr_n = [], None
    for line in text:
        s = line.strip()
        if not s:
            continue
        if s in ("NUISANCE_PARAMETERS", "CORRELATION_MATRIX", "NLL"):
            section = s
            continue
        if section == "NUISANCE_PARAMETERS":
            m = re.match(r"^(\S+)\s+(\S+)\s+\+(\S+)\s+-(\S+)", s)
            if m:
                out["nps"][m.group(1)] = (float(m.group(2)), float(m.group(3)), float(m.group(4)))
        elif section == "CORRELATION_MATRIX":
            nums = s.split()
            if corr_n is None and len(nums) == 2:
                corr_n = int(nums[0])
            else:
                corr_rows.append([float(x) for x in nums])
        elif section == "NLL":
            out["nll"] = float(s)
    if corr_n:
        out["corr"] = {"n": corr_n, "matrix": corr_rows}
    return out


def parse_err_decomp(path) -> dict:
    """{KEY: (sym, up, down)} from the tab-separated error decomposition."""
    out = {}
    for line in Path(path).read_text().splitlines():
        parts = line.strip().split("\t")
        if len(parts) == 4:
            out[parts[0]] = tuple(float(x) for x in parts[1:])
    return out


def parse_grouped(path) -> dict:
    """{Category: (impact, up, down)} from GroupedImpact / group_errDecomp files."""
    out = {}
    rx = re.compile(r"^(\S.*?)\s{2,}(\S+)\s+\(\s*\+(\S+),\s*-(\S+)\s*\)")
    for line in Path(path).read_text().splitlines():
        m = rx.match(line.strip())
        if m:
            out[m.group(1)] = (float(m.group(2)), float(m.group(3)), float(m.group(4)))
    return out


def parse_ranking(path) -> list:
    """List of dicts sorted as in the file (largest post-fit impact first)."""
    rows = []
    for line in Path(path).read_text().splitlines():
        parts = line.split()
        if len(parts) >= 8:
            try:
                nums = [float(p.lstrip("+")) for p in parts[1:8]]
            except ValueError:
                continue
            rows.append({"name": parts[0], "pull": nums[0], "err_up": nums[1], "err_down": nums[2],
                         "dpoi_up_post": nums[3], "dpoi_down_post": nums[4],
                         "dpoi_up_pre": nums[5], "dpoi_down_pre": nums[6]})
    return rows


def parse_tables_yaml(path) -> dict:
    import yaml
    with open(path) as fh:
        return yaml.safe_load(fh)


def parse_gof(log_path):
    """Saturated-model goodness of fit from the `f` step log: (chi2/ndf string, probability)."""
    text = Path(log_path).read_text() if Path(log_path).exists() else ""
    prob = re.findall(r"probability\s*=\s*([0-9.eE+-]+)", text)
    chi2 = re.findall(r"chi2\s*/\s*ndf\s*=\s*([0-9.eE+-]+)\s*/\s*(\d+)", text)
    return {"gof_probability": float(prob[-1]) if prob else None,
            "gof_chi2_ndf": chi2[-1] if chi2 else None}


def summarise(job_dir, job, poi="mu_Z", stat_only_suffix="_statOnly") -> dict:
    """Collect the POI result and uncertainty breakdown of one TRExFitter job directory."""
    job_dir = Path(job_dir)
    fits = job_dir / "Fits"
    res = {"job": job, "poi": poi}
    full = parse_fit_txt(fits / f"{job}.txt")
    if poi not in full["nps"]:
        raise RuntimeError(f"{poi} not found in {fits / (job + '.txt')}")
    val, up, down = full["nps"][poi]
    res.update(poi_value=val, poi_err_up=up, poi_err_down=down, poi_err=0.5 * (up + down),
               nps={k: v for k, v in full["nps"].items() if k != poi}, nll=full["nll"])
    dec = fits / f"{job}_errDecomp_{poi}.txt"
    if dec.exists():
        d = parse_err_decomp(dec)
        res["err_decomp"] = {k: v[0] for k, v in d.items()}
    grouped = fits / f"GroupedImpact_{poi}.txt"
    if not grouped.exists():
        grouped = fits / f"{job}_group_errDecomp_{poi}.txt"
    if grouped.exists():
        res["grouped_impact"] = {k: v[0] for k, v in parse_grouped(grouped).items()}
        res["grouped_impact_file"] = grouped.name
    stat = fits / f"{job}{stat_only_suffix}.txt"
    if stat.exists():
        s = parse_fit_txt(stat)["nps"].get(poi)
        if s:
            res["poi_stat_only"] = {"value": s[0], "err_up": s[1], "err_down": s[2], "err": 0.5 * (s[1] + s[2])}
    rank = fits / f"NPRanking_{poi}.txt"
    if not rank.exists():
        rank = fits / "NPRanking.txt"
    if rank.exists():
        res["ranking"] = parse_ranking(rank)
    for kind in ("prefit", "postfit"):
        tab = job_dir / "Tables" / f"Table_{kind}.yaml"
        if tab.exists():
            try:
                res[f"table_{kind}"] = parse_tables_yaml(tab)
            except Exception as exc:  # pragma: no cover - yaml oddities
                res[f"table_{kind}_error"] = str(exc)
    return res
