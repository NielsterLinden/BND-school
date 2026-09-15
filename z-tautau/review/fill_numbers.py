#!/usr/bin/env python
"""Fill the @@TOKEN@@ placeholders of README.md, handoff.md, REVIEW.md, docs/06, docs/07 and slides/*.tex
from output/results.json (run after step 6). Idempotent: tokens that are gone stay gone."""
import json, re, sys
from pathlib import Path
HERE = Path(__file__).resolve().parents[1]
r = json.loads((HERE / "output/results.json").read_text())
f = r["fit"][r["nominal_variant"]]; fo = r["fit"]["nosub"]
b = r["bdt"]; y = r["yields_prefit"][r["nominal_variant"]]; yr = r["yields_prefit_per_region"]
s60, sf = f["sigma_60_120_pb"], f["sigma_fid_pb"]
gi = f["grouped_impact"]
def pm(v, up, dn, n=3): return f"{v:.{n}f} +{up:.{n}f} -{dn:.{n}f}"
def pmtex(v, up, dn, n=3): return f"{v:.{n}f}^{{+{up:.{n}f}}}_{{-{dn:.{n}f}}}"
cats = ["tautau_SR0", "tautau_SR1", "tautau_SR2"]
sb2 = yr["tautau_SR2"]["DYtautau"]["value"] / max(yr["tautau_SR2"]["Fakes"]["value"], 1)
closure = ", ".join(f"{k.split('_')[-2]}/{k.split('_')[-1]} {100*v['delta']:.0f}\\,\\%" for k, v in r["closure_nps"].items())
closure_md = ", ".join(f"{k.split('_')[-2]}/{k.split('_')[-1]} {100*v['delta']:.0f} %" for k, v in r["closure_nps"].items())
catrows = " \\\\\n".join(f"SR{i} ({e0:.2f}--{e1:.2f}) & {yr[c]['Data']['value']:.0f} & {yr[c]['Fakes']['value']:.0f} & {yr[c]['DYtautau']['value']:.0f} & {yr[c]['DYtautau_nonfid']['value']:.0f} & {yr[c]['DYtautau']['value']/max(yr[c]['Fakes']['value'],1):.2f}"
                         for i, (c, e0, e1) in enumerate(zip(cats, r["category_edges"][:-1], r["category_edges"][1:]))) + " \\\\"
impacts_md = "| source | impact |\n|---|---:|\n" + "\n".join(f"| {k} | {100*v:.1f} % |" for k, v in sorted(gi.items(), key=lambda kv: -kv[1]) if k != "FullSyst") + f"\n| **data statistics** | **{100*f['mu_stat']:.1f} %** |"
n_bkg = sum(v["value"] for k, v in y.items() if k not in ("Data", "Total", "DYtautau"))
tok = {
    "MU": pm(f["mu"], f["mu_err_up"], f["mu_err_down"]), "MUNOSUB": pm(fo["mu"], fo["mu_err_up"], fo["mu_err_down"]) if fo else "n/a",
    "SIG": f"{s60['value']:.0f}\\pm{s60['stat']:.0f}\\,^{{+{s60['err_up']:.0f}}}_{{-{s60['err_down']:.0f}}}\\pm{s60['acceptance']:.0f}",
    "SIGFULL": f"{s60['value']:.0f} ± {s60['stat']:.0f} (stat) +{s60['err_up']:.0f}/−{s60['err_down']:.0f} (syst+stat) ± {s60['acceptance']:.0f} (acc)",
    "SIGFID": f"{sf['value']:.2f}\\pm{sf['stat']:.2f}\\pm{sf['syst']:.2f}", "SIGFIDFULL": f"{sf['value']:.2f} ± {sf['stat']:.2f} (stat) ± {sf['syst']:.2f} (syst)",
    "SIG60": f"{s60['value']:.0f} +{s60['err_up']:.0f} −{s60['err_down']:.0f}", "SIG60NOSUB": f"{fo['sigma_60_120_pb']['value']:.0f}" if fo else "n/a",
    "SIGFIDV": f"{sf['value']:.2f} ± {sf['stat']:.2f} ± {sf['syst']:.2f}",
    "EXP": f"+{f['mu_expected_asimov']['err_up']:.3f} −{f['mu_expected_asimov']['err_down']:.3f}" if f.get("mu_expected_asimov") else "n/a",
    "OBS": f"+{f['mu_err_up']:.3f} −{f['mu_err_down']:.3f}", "STAT": f"{f['mu_stat']:.3f}",
    "GOF": f"{f['gof_probability']:.2f}" if f.get("gof_probability") is not None else "n/a", "GOFNOSUB": f"{fo['gof_probability']:.2f}" if fo and fo.get("gof_probability") is not None else "n/a",
    "AUC": f"{sum(b['training']['auc_test'])/len(b['training']['auc_test']):.3f}", "SB2": f"{sb2:.1f}",
    "NSIGTRAIN": f"{b['training']['n_sig']:,}", "NBKGTRAIN": f"{b['training']['n_bkg']:,}",
    "CLO": f"{r['sigmodel']['C_LO_over_NLO_fiducial']:.3f}" if r.get("sigmodel") else "n/a",
    "C": f"{r['for_combination']['C']:.4f}", "ACCEFF": f"{r['for_combination']['acc_eff']:.4g}", "NBKG": f"{n_bkg:.0f}",
    "NFAKE": f"{y['Fakes']['value']:.0f}", "NSIG": f"{y['DYtautau']['value']:.0f}", "NNONFID": f"{y['DYtautau_nonfid']['value']:.0f}",
    "NTT": f"{y['TTbar']['value']:.0f}", "NWJ": f"{y['WJets']['value']:.0f}", "NEE": f"{y['DYee']['value']:.0f}",
    "NVV": f"{y['WW']['value']+y['WZ']['value']+y['ZZ']['value']:.0f}",
    "NNP": str(len(f["pulls"])), "GAMMAS": f"{100*gi.get('Gammas', 0):.1f}\\,\\%", "TAUID": f"{100*gi.get('Tau ID', 0):.1f}\\,\\%",
    "SIGMODELIMPACT": "0 (removed)", "CATROWS": catrows, "CLOSURENPS": closure,
    "FAKEPULL": f"\\texttt{{FakeOSSS}} pulled {f['pulls']['FakeOSSS_tautau'][0]:+.2f}$\\sigma$",
    "FINDING1": f"Consistent with NNLO and with $\\mu\\mu$ ($\\mu_Z = 0.990$) within $1\\sigma$; data statistics {100*f['mu_stat']:.1f}\\,\\%, total $^{{+{100*f['mu_err_up']:.0f}}}_{{-{100*f['mu_err_down']:.0f}}}\\,\\%$.",
    "IMPACTS": impacts_md, "CLOSURE": closure_md,
    "COSSS": f"{r['fake_factors']['mcsub']['C_OS_SS_inclusive']:.3f}",
    "CTABLE": "\n".join(f"| Run2016{e} | " + " | ".join(f"{c:.3f}" for c in row) + " |" for e, row in zip("GH", r['fake_factors']['mcsub']['C_OS_SS_per_njet'])),
}
tok["MU"] = tok["MU"]; mu_tex = pmtex(f["mu"], f["mu_err_up"], f["mu_err_down"])
munosub_tex = pmtex(fo["mu"], fo["mu_err_up"], fo["mu_err_down"]) if fo else "n/a"
for path in ["README.md", "handoff.md", "REVIEW.md", "docs/05-fake-factors.md", "docs/06-cross-section.md", "docs/07-corrections-and-systematics.md", "docs/08-fit-and-results.md", "docs/00-overview.md", "docs/02-selection.md", "slides/ztautau_slides.tex"]:
    p = HERE / path
    if not p.exists():
        continue
    s = p.read_text()
    t = dict(tok)
    if path.endswith(".tex"):
        t["MU"], t["MUNOSUB"] = mu_tex, munosub_tex
    else:
        t["CLOSURENPS"] = closure_md; t["GAMMAS"] = f"{100*gi.get('Gammas', 0):.1f} %"; t["TAUID"] = f"{100*gi.get('Tau ID', 0):.1f} %"
        t["SIG"] = tok["SIGFULL"]
    if path == "README.md":
        t["RESULT"] = (f"σ(pp → Z/γ* → ττ, 60 < m < 120 GeV) = {tok['SIGFULL']} pb     NNLO: 1945 pb\n"
                       f"σ_fid(τhτh, vis. pT > 40 GeV, |η| < 2.1, 60 < m < 120 GeV) = {tok['SIGFIDFULL']} pb   pred.: 4.50 pb\n"
                       f"μ_Z = {tok['MU']}   (without MC subtraction in the FF regions: {tok['MUNOSUB']})")
    n = 0
    for k, v in t.items():
        s, c = re.subn(f"@@{k}@@", lambda m, v=v: v, s); n += c
    left = re.findall(r"@@[A-Z0-9]+@@", s)
    if "--dry-run" in sys.argv:
        print(f"{path}: would fill {n}" + (f", LEFT: {sorted(set(left))}" if left else ""))
        continue
    p.write_text(s)
    print(f"{path}: {n} filled" + (f", LEFT: {sorted(set(left))}" if left else ""))
if "--dry-run" in sys.argv:
    for k, v in tok.items():
        print(f"  {k:14s} = {str(v)[:110]!r}")
