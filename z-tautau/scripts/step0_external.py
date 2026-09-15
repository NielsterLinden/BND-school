#!/usr/bin/env python
"""Step 0 -- file lists and external inputs (docs/07-corrections-and-systematics.md).

    python scripts/step0_external.py

1. File lists of every sample (dCache manifest or CERN Open Data API) -> filelists/*.tsv.
2. TauPOG scale factors for UL2016 postVFP (DeepTau2017v2p1), from the public GitHub repositories
   (the CMS cvmfs area and the CERN GitLab jsonpog-integration need a CERN login, GitHub does not):
     cms-tau-pog/TauIDSFs      TauID_SF_dm_DeepTau2017v2p1VSjet_UL2016_postVFP.root   (tau ID, per DM)
                               TauID_SF_pt_DeepTau2017v2p1VSjet_UL2016_postVFP.root   (tau ID, per pT; cross-check)
                               TauES_dm_DeepTau2017v2p1VSjet_UL2016_postVFP.root      (tau energy scale)
                               TauID_SF_eta_DeepTau2017v2p1VSe_UL2016_postVFP.root    (e -> tau_h fakes)
                               TauID_SF_eta_DeepTau2017v2p1VSmu_UL2016_postVFP.root   (mu -> tau_h fakes)
     cms-tau-pog/TauTriggerSFs 2016ULpostVFP_tauTriggerEff_DeepTau2017v2p1.root       (di-tau leg efficiencies)
   The raw files go to $BND_TAUTAU_CACHE/downloads; the numbers the analysis uses are written to the
   small, committed external/tau_pog_UL2016postVFP.json, so no ROOT and no network are needed later.
3. Checks the luminosity-by-lumisection table used for the pileup profile.
"""

from __future__ import annotations

import json
import subprocess
import sys
import urllib.request
from datetime import date
from pathlib import Path

import numpy as np
import uproot

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from ztautau import config, samples  # noqa: E402

TAUID = "https://raw.githubusercontent.com/cms-tau-pog/TauIDSFs/master/data/{}"
TRIG = "https://raw.githubusercontent.com/cms-tau-pog/TauTriggerSFs/run2_SFs/data/{}"
FILES = {
    "id_dm": TAUID.format("TauID_SF_dm_DeepTau2017v2p1VSjet_UL2016_postVFP.root"),
    "id_pt": TAUID.format("TauID_SF_pt_DeepTau2017v2p1VSjet_UL2016_postVFP.root"),
    "tes": TAUID.format("TauES_dm_DeepTau2017v2p1VSjet_UL2016_postVFP.root"),
    "vse": TAUID.format("TauID_SF_eta_DeepTau2017v2p1VSe_UL2016_postVFP.root"),
    "vsmu": TAUID.format("TauID_SF_eta_DeepTau2017v2p1VSmu_UL2016_postVFP.root"),
    "trigger": TRIG.format("2016ULpostVFP_tauTriggerEff_DeepTau2017v2p1.root"),
}
TRIG_GRID = np.arange(20.0, 300.5, 0.5)       # the fitted curves are flat above ~150 GeV


def download(url: str) -> Path:
    dest = config.DOWNLOAD_DIR / url.rsplit("/", 1)[1]
    if not dest.exists():
        dest.parent.mkdir(parents=True, exist_ok=True)
        with urllib.request.urlopen(url, timeout=120) as resp:
            dest.write_bytes(resp.read())
    return dest


def th1(obj):
    return {"edges": obj.axis().edges().tolist(), "values": obj.values().tolist(), "errors": obj.errors().tolist()}


def pt_binned_id(path: Path) -> dict:
    """The pT-binned SFs are TF1 formulas; evaluate them with ROOT at a few pT values (cross-check only)."""
    code = (f"import ROOT, json; f = ROOT.TFile.Open('{path}'); out = {{}}\n"
            "for wp in ('VVVLoose', 'Medium'):\n"
            "    out[wp] = {v: [f.Get(wp + '_' + v).Eval(x) for x in (45., 100., 300.)] for v in ('cent', 'up', 'down')}\n"
            "print(json.dumps(out))")
    res = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
    if res.returncode != 0:
        print("  (ROOT evaluation of the pT-binned SFs failed; skipped)")
        return {}
    return json.loads(res.stdout.strip().splitlines()[-1])


def main():
    for key in samples.SAMPLES:
        samples.write_filelist(key)
    print(f"file lists written to {config.FILELIST_DIR}")

    paths = {k: download(u) for k, u in FILES.items()}
    out = {"retrieved": str(date.today()), "sources": FILES, "era": "UL2016 postVFP (Run2016G+H)",
           "tau_id": "DeepTau2017v2p1"}
    f = uproot.open(paths["id_dm"])
    out["id_vsjet_dm"] = {}
    for wp in ("VVVLoose", "VVLoose", "VLoose", "Loose", "Medium", "Tight", "VTight", "VVTight"):
        h = f[wp]
        out["id_vsjet_dm"][wp] = {str(int(dm)): [float(h.values()[int(dm)]), float(h.errors()[int(dm)])]
                                  for dm in config.TAU_DMS}
    out["id_vsjet_pt_check"] = {"pt": [45.0, 100.0, 300.0], **pt_binned_id(paths["id_pt"])}
    h = uproot.open(paths["tes"])["tes"]
    out["tes_dm"] = {str(dm): [float(h.values()[dm]), float(h.errors()[dm])] for dm in config.TAU_DMS}
    out["vse"] = {wp: th1(uproot.open(paths["vse"])[wp]) for wp in ("VVLoose",)}
    out["vsmu"] = {wp: th1(uproot.open(paths["vsmu"])[wp]) for wp in ("VLoose",)}
    t = uproot.open(paths["trigger"])
    for wp in ("Medium", "Tight", "VTight"):          # offline VSjet working point of the legs
        trig = {"pt": TRIG_GRID.tolist()}
        for dm in config.TAU_DMS:
            entry = {}
            for kind in ("sf", "data", "mc"):
                hist = t[f"{kind}_ditau_{wp}_dm{dm}_fitted"]
                idx = np.clip(np.searchsorted(hist.axis().edges(), TRIG_GRID, side="right") - 1, 0, len(hist.values()) - 1)
                entry[kind] = hist.values()[idx].round(5).tolist()
                entry[kind + "_err"] = hist.errors()[idx].round(5).tolist()
            trig[str(dm)] = entry
        out[f"trigger_ditau_{wp}"] = trig
    trig = out["trigger_ditau_Medium"]
    config.EXTERNAL_DIR.mkdir(parents=True, exist_ok=True)
    dest = config.EXTERNAL_DIR / "tau_pog_UL2016postVFP.json"
    dest.write_text(json.dumps(out, indent=1))
    print(f"TauPOG inputs -> {dest} ({dest.stat().st_size / 1e3:.0f} kB)")
    for dm in config.TAU_DMS:
        sf, err = out["id_vsjet_dm"][config.TAU_WP][str(dm)]
        tes, terr = out["tes_dm"][str(dm)]
        i = int(np.searchsorted(TRIG_GRID, 50.0))
        print(f"  DM{dm:>2d}: ID SF {sf:.3f} +- {err:.3f}   TES {tes:.3f} +- {terr:.3f}   "
              f"trigger SF(50 GeV) {trig[str(dm)]['sf'][i]:.3f} +- {trig[str(dm)]['sf_err'][i]:.3f}")
    ok = config.LUMIBYLS_CSV.exists()
    print(f"lumi-by-LS table {config.LUMIBYLS_CSV}: {'present' if ok else 'MISSING -> xrdcp ' + config.LUMIBYLS_URL}")


if __name__ == "__main__":
    main()
