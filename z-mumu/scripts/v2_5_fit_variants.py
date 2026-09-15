#!/usr/bin/env python
"""v2 step 5b -- stability of the fitted signal strength against the fit configuration.

    python scripts/v2_5_fit_variants.py

Re-runs scripts/v2_5_fit.py (quick mode: h w f i) with different signal-region binnings, without
the powheg SigModel template and with the template smoothing of the first v2 fit, and collects
mu_Z, its uncertainty, the goodness of fit and the notable pulls / constraints into
fit/results/stability.json (table in output/v2/RESULTS_v2.md). The nominal fit is the 5 GeV one;
REVIEW.md F3 explains why the 2 GeV configuration of the first v2 result was not robust.
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from zmumu import config

FIT = config.REPO_DIR / "fit"
VARIANTS = [
    ("zmumu", "nominal: 12 x 5 GeV, SigModel two-sided (LHE window), no smoothing", None),
    ("stab_2gev", "30 x 2 GeV", ["--rebin", "2"]),
    ("stab_10gev", "6 x 10 GeV", ["--rebin", "10"]),
    ("stab_1bin", "1 bin (counting)", ["--rebin", "60"]),
    ("stab_nosig", "12 x 5 GeV, no SigModel", ["--no-sigmodel"]),
    ("stab_smooth", "12 x 5 GeV, MuonScale/MuonRes smoothed", ["--smoothing"]),
    ("stab_2gev_nosig", "30 x 2 GeV, no SigModel", ["--rebin", "2", "--no-sigmodel"]),
]


def notable(nps, pull_min=1.0, constr_max=0.7):
    out = []
    for name, (val, up, down) in nps.items():
        if name.startswith("gamma") or name == "mu_Z":
            continue
        c = 0.5 * (up + down)
        if abs(val) >= pull_min or c <= constr_max:
            out.append(f"{name} {val:+.1f} ({c:.2f})")
    return out


def main():
    rows = []
    for tag, label, opts in VARIANTS:
        if opts is None:
            res_path = FIT / "results" / f"{tag}_fit_result.json"
            if not res_path.exists():
                sys.exit("run scripts/v2_5_fit.py first")
        else:
            print(f"[variants] {tag}: {label}", flush=True)
            cmd = [sys.executable, str(HERE / "v2_5_fit.py"), "--tag", tag, "--quick", *opts]
            if subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT).returncode != 0:
                sys.exit(f"{tag} failed")
            res_path = FIT / "results" / f"{tag}_fit_result.json"
        r = json.load(open(res_path))
        nps = {k: tuple(v) for k, v in r["nps"].items()}
        rows.append({"tag": tag, "label": label, "mu": r["mu"], "err_up": r["mu_err_up"], "err_down": r["mu_err_down"],
                     "sigma_fid_pb": r["sigma_fid_pb"], "gof_p": r["gof"].get("gof_probability"),
                     "syst_total": r["grouped_impacts_mu"].get("FullSyst"), "notable": notable(nps)})
    with open(FIT / "results" / "stability.json", "w") as fh:
        json.dump(rows, fh, indent=1)
    mus = [x["mu"] for x in rows]
    print(f"[variants] mu_Z spread: {min(mus):.4f} - {max(mus):.4f} (half-spread {50*(max(mus)-min(mus)):.2f}%)")
    for x in rows:
        print(f"  {x['label']:<55} mu = {x['mu']:.4f} +{x['err_up']:.4f} -{x['err_down']:.4f}  p = {x['gof_p']}  {'; '.join(x['notable'][:4])}")


if __name__ == "__main__":
    main()
