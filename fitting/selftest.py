#!/usr/bin/env python
"""End-to-end self-test of the fitting helpers: toy histograms -> fit inputs -> config ->
trex-fitter h/w/f -> parsed result. Run after building TRExFitter:

    source setup.sh && python fitting/selftest.py [workdir]
"""

from __future__ import annotations

import sys
from pathlib import Path

import hist
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from fitting import run_trex, trexconfig as tc, trexhist


def main():
    work = Path(sys.argv[1] if len(sys.argv) > 1 else "/tmp/trex_selftest").resolve()
    (work / "fitinputs").mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(1)
    ax = hist.axis.Regular(30, 60, 120, name="mll")

    def make(n, mean, sigma, w=1.0):
        h = hist.Hist(ax, storage=hist.storage.Weight())
        h.fill(rng.normal(mean, sigma, n), weight=w)
        return h

    sig = make(200000, 91.0, 3.0, w=0.5)          # 100k expected
    bkg = make(20000, 90.0, 30.0, w=0.1)          # 2k expected, flat-ish
    data = hist.Hist(ax, storage=hist.storage.Double())
    data.fill(rng.normal(91.0, 3.0, 100500)); data.fill(rng.normal(90.0, 30.0, 2000))
    sig_up = make(200000, 91.0, 3.0, w=0.5 * 1.02)
    sig_dn = make(200000, 91.0, 3.0, w=0.5 * 0.98)
    hists = {"toy_SR__Data": data, "toy_SR__DYmumu": sig, "toy_SR__TTbar": bkg,
             "toy_SR__DYmumu__MuonIDUp": sig_up, "toy_SR__DYmumu__MuonIDDown": sig_dn}
    rep = trexhist.write_fitinputs(work / "fitinputs/toy.root", hists, meta={"lumi_pb": 1.0})
    trexhist.check_fitinputs(work / "fitinputs/toy.root", list(hists))
    print("fit inputs:", rep)

    blocks = [
        tc.job("toy", POI="mu_Z", ReadFrom="HIST", HistoPath="fitinputs", HistoFile="toy",
               OutputDir="results", MCstatThreshold=0, UseGammaPulls=True, DebugLevel=1,
               ImageFormat="png"),
        tc.fit("fit", FitType="SPLUSB", FitRegion="CRSR", UseMinos="mu_Z", NumCPU=2),
        tc.region("toy_SR", Type="SIGNAL", HistoName="toy_SR", VariableTitle="m [GeV]", Label="toy SR"),
        tc.sample("Data", Type="DATA", Title="Data", HistoNameSuff="__Data"),
        tc.sample("DYmumu", Type="SIGNAL", Title="signal", HistoNameSuff="__DYmumu", FillColor=800),
        tc.sample("TTbar", Type="BACKGROUND", Title="bkg", HistoNameSuff="__TTbar", FillColor=632),
        tc.normfactor("mu_Z", Title="#mu_{Z}", Nominal=1, Min=0.5, Max=1.5, Samples="DYmumu"),
        tc.overall_syst("Lumi", 0.012, -0.012, "DYmumu,TTbar", "Luminosity", title="Luminosity"),
        tc.overall_syst("XS_TTbar", 0.06, -0.06, "TTbar", "Background normalisation"),
        tc.histo_syst("MuonID", "DYmumu", "Muon efficiency", title="Muon ID SF"),
    ]
    (work / "toy.config").write_text(tc.render(blocks, header="fitting/ self-test"))
    for actions in ("h", "w", "f", "dp", "i"):
        run_trex.run("toy.config", actions, cwd=work, log=f"results/toy/logs/{actions}.log")
    run_trex.run("toy.config", "wf", options="StatOnly=TRUE:Suffix=_statOnly", cwd=work,
                 log="results/toy/logs/statonly.log")
    res = run_trex.summarise(work / "results/toy", "toy", poi="mu_Z")
    print(f"mu_Z = {res['poi_value']:.4f} +{res['poi_err_up']:.4f} -{res['poi_err_down']:.4f}"
          f"   stat-only {res.get('poi_stat_only', {}).get('err', float('nan')):.4f}")
    print("grouped impacts:", res.get("grouped_impact"))
    print("err decomposition:", res.get("err_decomp"))
    print("GoF:", run_trex.parse_gof(work / "results/toy/logs/f.log"))
    assert abs(res["poi_value"] - 1.005) < 0.02, "toy fit did not recover mu ~ 1.005"
    print("SELFTEST OK")


if __name__ == "__main__":
    main()
