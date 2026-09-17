#!/usr/bin/env python
"""TRExFitter MultiFit of Z -> ee, Z -> mumu and Z -> tautau (tau_h tau_h, mu tau_h, e tau_h, e mu); CMS Open Data 2016 G+H.

    source ../../setup.sh
    python run.py prepare          # ee inputs from z-ee/Zee_fit.tar.gz; channel and MultiFit configs in work/
    python run.py condor --submit  # every fit below as one HTCondor DAG (mf/condor.py), ending with `results --interim`
    python run.py results          # -> output/result.json      (only once the channel results are final)
    python run.py plots            # -> output/plots/

The same chain on one machine, in order (hours with the four-channel tautau likelihood; the DAG is the normal route):

    python run.py all              # prepare workspaces fits variations impacts results plots

    python run.py workspaces       # trex-fitter hw for every channel config
    python run.py fits             # combined fit, three-POI fit, standalone channel fits, stat-only fit, likelihood scan
    python run.py variations       # the alternative likelihoods of config/channels.json "variations" and its "checks"
    python run.py impacts          # NP ranking of the combined fit (one refit job per parameter, parallel)

One job of the DAG (what condor/run_step.sh calls):

    python run.py likelihood <name>   # workspaces + fit of one likelihood (common, split, var_<name>, check_<name>)
    python run.py channelfit <key>    # standalone fit of one channel in the common likelihood + its stat-only workspace
    python run.py statonly            # the stat-only combined fit (needs every channelfit)
    python run.py scan                # the likelihood scan of the combined fit, in work/common/combination_scan/
    python run.py rank <NP> [<NP>..]  # ranking refits of single parameters
    python run.py mergeranking        # NPRanking_<NP>_mu_Z.txt -> NPRanking_mu_Z.txt

Likelihoods (work/<name>/):
    common       one POI mu_Z for all channels: the combined cross section sigma = mu_Z x poi.reference_pb
    split        mu_Z_ee, mu_Z_mumu, mu_Z_tautau with all shared nuisance parameters profiled together:
                 per-channel cross sections and the compatibility test -2 ln(L_common / L_split)
    var_<name>   one per entry of "variations", common POI
    check_<name> one per entry of "checks": a single channel, fitted on its own
TRExFitter output goes to work/ (git-ignored); what is kept is output/.
"""

from __future__ import annotations

import argparse
import copy
import json
import re
import subprocess
import time

from mf import ee_input, provenance, trexcfg
from mf.paths import WORK, manifest, repo_path

from fitting import run_trex  # noqa: E402  (mf.paths puts the repository root on sys.path)


def fit_options(m):
    return {k: v for k, v in m.get("fit", {}).items() if not k.startswith("_")}


def _public(d):
    return {k: v for k, v in (d or {}).items() if not k.startswith("_")}


def likelihoods(m):
    """{name: {"model": common|split|channel, "channels": {key: spec}}} for every likelihood this folder fits."""
    out = {"common": {"model": "common", "channels": m["channels"]},
           "split": {"model": "split", "channels": m["channels"]}}
    for name, var in _public(m.get("variations")).items():
        chans = {k: copy.deepcopy(v) for k, v in m["channels"].items() if k in var.get("channels", m["channels"])}
        for k, split in var.get("split", {}).items():
            chans[k]["split_shape_norm"] = split
        for k, over in var.get("overall", {}).items():
            chans[k]["overall"] = over
        for k, over in var.get("channel_overrides", {}).items():
            chans[k].update(copy.deepcopy(over))
        out[f"var_{name}"] = {"model": "common", "channels": chans, "label": var["label"],
                              "fit": {**fit_options(m), **_public(var.get("fit"))}}
    for name, chk in _public(m.get("checks")).items():
        spec = copy.deepcopy(m["channels"][chk["channel"]])
        spec.update(copy.deepcopy(chk.get("overrides", {})))
        out[f"check_{name}"] = {"model": "channel", "channels": {chk["channel"]: spec}, "label": chk["label"]}
    return out


def prepare(args):
    m = manifest()
    poi = m["poi"]
    status = {"ee_input": ee_input.extract(m["channels"]["ee"]["inputs"], WORK / "inputs" / "ee"),
              "inputs": provenance.record(m), "likelihoods": {}}
    if abs(status["ee_input"]["sigma_reference_pb"] - m["channels"]["ee"]["sigma_reference_pb"]) > 1e-6:
        raise ValueError("config/channels.json: ee sigma_reference_pb differs from the tarball's reference_cross_section.txt")
    for name, lk in likelihoods(m).items():
        wdir = WORK / name
        wdir.mkdir(parents=True, exist_ok=True)
        fits, pois, info = [], [], {}
        for key, spec in lk["channels"].items():
            pname = f'{poi["name"]}_{key}' if lk["model"] == "split" else poi["name"]
            histo_path = WORK / "inputs" / "ee" if key == "ee" else repo_path(spec["histo_path"])
            blocks, info[key] = trexcfg.adapt_channel(key, spec, poi, pname, wdir, histo_path, lk.get("fit", fit_options(m)))
            cfg = wdir / f"{key}.config"
            cfg.write_text(trexcfg.render(blocks, header=f"Generated by combination/combLieke/run.py prepare from {spec['config']}\n"
                                                         f"(likelihood '{name}'). Edit config/channels.json or mf/trexcfg.py, not this file."))
            fits.append({"job": key, "config": str(cfg), "directory": str(wdir / key), "label": spec["label"]})
            pois.append(pname)
        mpois = pois if lk["model"] == "split" else [poi["name"]]
        status["likelihoods"][name] = {"model": lk["model"], "pois": mpois, "channels": info, "label": lk.get("label", name)}
        if lk["model"] == "channel":
            continue

        def mf(job, options, scan=False):
            return trexcfg.render(trexcfg.multifit(job, fits, mpois, wdir, {**poi, "fit": options}, scan=scan),
                                  header=f"Generated by combination/combLieke/run.py prepare (likelihood '{name}', MultiFit '{job}').")
        (wdir / "multifit.config").write_text(mf("combination", lk.get("fit", fit_options(m))))
        if name == "common":
            # the scan is its own MultiFit (combination_scan/): 41 refits, which nothing else has to wait for
            (wdir / "multifit_scan.config").write_text(mf("combination_scan", fit_options(m), scan=True))
            # the stat-only fit (StatOnly=TRUE:Suffix=_statOnly on the command line)
            (wdir / "multifit_statonly.config").write_text(mf("combination", fit_options(m)))
            # the ranking refits read the nominal fit (Fits/combination.root) and use their own minimiser settings
            (wdir / "multifit_ranking.config").write_text(mf("combination", {**fit_options(m), **_public(m.get("ranking_fit"))}))
    (WORK / "status.json").write_text(json.dumps(status, indent=1) + "\n")
    print(f"prepared {len(status['likelihoods'])} likelihoods in {WORK}; ee inputs: {status['ee_input']['histograms']} histograms, "
          f"max. relative difference to {status['ee_input'].get('check_against')}: {status['ee_input'].get('max_relative_difference')}")


STAT_ONLY = "StatOnly=TRUE:Suffix=_statOnly"


def trex(config, actions, log, options=None):
    """trex-fitter <actions> <config> [<options>]; full output to `log`, the fit summary lines echoed."""
    t0 = time.time()
    print(f"trex-fitter {actions} {config.relative_to(WORK)}" + (f" {options}" if options else ""), flush=True)
    log.parent.mkdir(parents=True, exist_ok=True)
    with open(log, "w") as fh:
        code = subprocess.run([run_trex.trex_binary(), actions, str(config)] + ([options] if options else []),
                              cwd=config.parent, stdout=fh, stderr=subprocess.STDOUT).returncode
    keep = re.compile(r"minuit status|hess status|minos status|Edm =|probability =|Fit failure")
    summary = [re.sub(r"\x1b\[[0-9;]*m", "", line).split("|", 1)[-1].strip()
               for line in log.read_text(errors="replace").splitlines() if keep.search(line)]
    print("   " + "; ".join(summary[-5:]) + f"  [exit {code}, {time.time() - t0:.0f} s]", flush=True)
    if code != 0:
        raise RuntimeError(f"trex-fitter {actions} {config} failed (exit {code}); see {log}")


def _names(args, prefix=None):
    names = list(likelihoods(manifest()))
    if prefix is not None:
        names = [n for n in names if n.startswith(prefix)]
    return [n for n in names if not args.only or n in args.only]


def _workspaces(name):
    for cfg in sorted((WORK / name).glob("*.config")):
        if not cfg.name.startswith("multifit"):
            trex(cfg, "hw", WORK / name / "logs" / f"{cfg.stem}_hw.log")


def _fit(name):
    if (WORK / name / "multifit.config").exists():
        trex(WORK / name / "multifit.config", "mwf", WORK / name / "logs" / "multifit_mwf.log")
    else:                                                   # a single-channel check
        for cfg in sorted((WORK / name).glob("*.config")):
            trex(cfg, "f", WORK / name / "logs" / f"{cfg.stem}_f.log")


def workspaces(args):
    for name in _names(args):
        _workspaces(name)


def likelihood(args):
    """Workspaces and fit of the named likelihoods; after `common`, the list of parameters to rank."""
    known = likelihoods(manifest())
    for name in args.names:
        if name not in known:
            raise SystemExit(f"unknown likelihood {name}; known: {', '.join(known)}")
        _workspaces(name)
        _fit(name)
        if name == "common":
            write_rank_params()


def channelfit(args):
    """Standalone fit of a channel in the common likelihood, then its stat-only workspace."""
    for key in args.names or list(manifest()["channels"]):
        trex(WORK / "common" / f"{key}.config", "f", WORK / "common" / "logs" / f"{key}_f.log")
        trex(WORK / "common" / f"{key}.config", "hw", WORK / "common" / "logs" / f"{key}_hw_statOnly.log", STAT_ONLY)


def statonly(args):
    """Data-statistics-only uncertainty: the combination of the stat-only channel workspaces (no NPs, no gammas)."""
    trex(WORK / "common" / "multifit_statonly.config", "mwf", WORK / "common" / "logs" / "multifit_mwf_statOnly.log", STAT_ONLY)


def scan(args):
    trex(WORK / "common" / "multifit_scan.config", "mwf", WORK / "common" / "logs" / "multifit_scan_mwf.log")


def fits(args):
    for name in [n for n in ("common", "split") if not args.only or n in args.only]:
        _fit(name)
    if not args.only or "common" in args.only:
        args.names = []
        channelfit(args)
        statonly(args)
        scan(args)


def variations(args):
    for name in _names(args, "var_") + _names(args, "check_"):
        _fit(name)


def ranked_parameters():
    """Every fitted parameter of the combined fit except the POI and the MC-statistics gammas."""
    poi = manifest()["poi"]["name"]
    fit = run_trex.parse_fit_txt(WORK / "common" / "combination" / "Fits" / "combination.txt")
    return [n for n in fit["nps"] if not n.startswith(("gamma", poi))]


def write_rank_params():
    path = WORK / "condor" / "params_rank.txt"
    path.parent.mkdir(parents=True, exist_ok=True)
    names = ranked_parameters()
    path.write_text("".join(f"rank {n}\n" for n in names))
    print(f"{len(names)} parameters to rank -> {path}", flush=True)


def _rank_one(name):
    wdir = WORK / "common"
    (wdir / "logs").mkdir(parents=True, exist_ok=True)
    with open(wdir / "logs" / f"rank_{name}.log", "w") as fh:
        return name, subprocess.run([run_trex.trex_binary(), "mr", str(wdir / "multifit_ranking.config"), f"Ranking={name}"],
                                    cwd=wdir, stdout=fh, stderr=subprocess.STDOUT).returncode


def rank(args):
    """Ranking refits of the named parameters (one HTCondor job each). TRExFitter selects by substring, so
    Ranking=PDF also refits Acc_PDF and PDF_eeShape into the same file; mergeranking keeps each parameter's own row."""
    for name in args.names:
        t0 = time.time()
        _, code = _rank_one(name)
        # no exception: DAGMan would remove every other ranking job of the cluster; mergeranking lists what is missing
        print(f"ranking {name}: {'FAILED, ' if code else ''}exit {code}, {time.time() - t0:.0f} s; work/common/logs/rank_{name}.log", flush=True)


def mergeranking(args):
    poi = manifest()["poi"]["name"]
    fits_dir = WORK / "common" / "combination" / "Fits"
    names, rows = ranked_parameters(), {}
    for n in names:
        path = fits_dir / f"NPRanking_{n}_{poi}.txt"
        for line in (path.read_text().splitlines() if path.exists() else []):
            if line.split() and line.split()[0] == n:
                rows[n] = line
    missing = [n for n in names if n not in rows]
    (fits_dir / f"NPRanking_{poi}.txt").write_text("\n".join(rows[n] for n in names if n in rows) + "\n")
    print(f"ranking: {len(rows)}/{len(names)} parameters" + (f"; missing: {missing}" if missing else ""), flush=True)


def impacts(args):
    """Refit-based ranking of every nuisance parameter except the MC-statistics gammas, one
    `trex-fitter mr ... Ranking=<NP>` per parameter in parallel, merged into Fits/NPRanking_mu_Z.txt.
    The refits run on multifit_ranking.config (config/channels.json "ranking_fit"); the post-fit values and
    errors they fix each parameter at come from the nominal combined fit of `fits`.

    The uncertainty groups come from the covariance decomposition that mwf writes: TRExFitter's refit-based
    grouped impacts (`mi`) fail HESSE in this likelihood (README)."""
    from concurrent.futures import ThreadPoolExecutor
    poi = manifest()["poi"]["name"]
    for old in (WORK / "common" / "combination" / "Fits").glob("NPRanking*"):
        old.unlink()
    t0 = time.time()
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        codes = dict(pool.map(_rank_one, ranked_parameters()))
    failed = [n for n, c in codes.items() if c != 0]
    print(f"ranking refits of {poi}: {time.time() - t0:.0f} s" + (f"; failed: {failed}" if failed else ""), flush=True)
    mergeranking(args)


def final(args):
    """The last node of the DAG: merge the ranking and collect what has been fitted, without touching output/."""
    mergeranking(args)
    args.interim = True
    results(args)


def condor(args):
    from mf import condor as cd
    cd.write(likelihoods(manifest()), list(manifest()["channels"]), submit=args.submit)


def results(args):
    from mf import results as res
    res.collect(interim=args.interim)


def plots(args):
    from mf import plots as pl
    pl.make_all()


STEPS = {"prepare": prepare, "workspaces": workspaces, "fits": fits, "variations": variations,
         "impacts": impacts, "results": results, "plots": plots}
JOBS = {"likelihood": likelihood, "channelfit": channelfit, "statonly": statonly, "scan": scan, "rank": rank,
        "mergeranking": mergeranking, "final": final, "condor": condor}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("step", choices=list(STEPS) + list(JOBS) + ["all"])
    ap.add_argument("names", nargs="*", help="likelihoods (likelihood), channels (channelfit) or parameters (rank)")
    ap.add_argument("--only", nargs="*", default=None, help="restrict workspaces/fits/variations to these likelihoods")
    ap.add_argument("--workers", type=int, default=6, help="parallel trex-fitter processes for the local ranking (impacts)")
    ap.add_argument("--interim", action="store_true", help="results: write interim/result.json instead of output/result.json")
    ap.add_argument("--submit", action="store_true", help="condor: also condor_submit_dag")
    args = ap.parse_args()
    for step in (STEPS if args.step == "all" else [args.step]):
        {**STEPS, **JOBS}[step](args)


if __name__ == "__main__":
    main()
