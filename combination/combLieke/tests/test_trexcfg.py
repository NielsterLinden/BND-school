"""Unit checks of the config adapter (no TRExFitter needed; the fit inputs of the channels must exist).

    source ../../setup.sh && python -m pytest -q tests      (or: python tests/test_trexcfg.py)
"""

from __future__ import annotations

import copy
import json
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))

from mf import trexcfg  # noqa: E402
from mf.paths import manifest, repo_path  # noqa: E402


def _adapt(key, name="mu_Z", **override):
    m = manifest()
    spec = copy.deepcopy(m["channels"][key])
    spec.update(override)
    histo = repo_path(spec["histo_path"]) if "histo_path" in spec else None
    return trexcfg.adapt_channel(key, spec, m["poi"], name, Path(tempfile.mkdtemp()), histo)


def test_roundtrip_every_channel_config():
    for spec in manifest()["channels"].values():
        blocks = trexcfg.read(repo_path(spec["config"]))
        with tempfile.NamedTemporaryFile("w", suffix=".config", delete=False) as fh:
            fh.write(trexcfg.render(blocks))
        again = trexcfg.read(Path(fh.name))
        assert [(b.kind, b.name, b.opts) for b in blocks] == [(b.kind, b.name, b.opts) for b in again]


def test_poi_and_reference_factor():
    m = manifest()
    blocks, info = _adapt("tautau", name="mu_Z_tautau")
    job = trexcfg.first(blocks, "Job")
    assert job.opts["POI"] == "mu_Z_tautau"
    nf = {b.name: b for b in blocks if b.kind == "NormFactor"}
    signal = trexcfg.listopt(nf["mu_Z_tautau"].opts["Samples"])
    assert "mu_Z" not in nf and len(signal) == 15 and all(s.startswith("DYtautau_tDM") for s in signal)
    assert nf["xsref_tautau"].opts["Samples"] == nf["mu_Z_tautau"].opts["Samples"]      # never DYtautau_out
    assert info["acceptance"] == [] and "mu_ttbar" in nf                                # acceptance inside the tautau fit
    factor = float(nf["xsref_tautau"].opts["Nominal"])
    assert abs(factor * m["channels"]["tautau"]["sigma_reference_pb"] - m["poi"]["reference_pb"]) < 1e-9
    assert nf["xsref_tautau"].opts["Constant"] == "TRUE"


def test_acceptance_nps_follow_the_channel_metadata():
    _, info = _adapt("mumu")
    acc = {a["name"]: a["rel"] for a in info["acceptance"]}
    assert set(acc) == {"Acc_PDF", "Acc_AlphaS", "Acc_QCDScale", "Acc_PTZ_mumu", "Acc_Generator_mumu", "Acc_PS_FSR",
                        "Acc_QEDFSR_mumu", "AccStat_mumu"}
    meta = json.loads(repo_path("z-mumu/fit/fitinputs/zmumu.root.meta.json").read_text())
    total = sum(v * v for v in acc.values()) ** 0.5
    assert abs(total - meta["acceptance"]["A_rel_unc"]) < 1e-9     # every row of the frozen acceptance block (17 Sep 2026)
    assert abs(total - 0.0070078) < 1e-6


def test_empty_bin_drop_is_checked():
    blocks, _ = _adapt("tautau")
    regions = {b.name: b for b in blocks if b.kind == "Region"}
    assert regions["tautau_SR2"].opts["DropBins"] == "1"           # z-tautau's own config drops its empty bin
    assert len(regions) == 13 and not any("SRlo" in r or "SRhi" in r for r in regions)   # never the pT-split copies
    try:
        _adapt("tautau", drop_empty_bins={"tautau_SR1": [5]})
    except ValueError as err:
        assert "not empty" in str(err)
    else:
        raise AssertionError("a non-empty bin was dropped")


def test_tautau_lepton_parameters_are_decorrelated():
    blocks, _ = _adapt("tautau")
    nps = {b.opts.get("NuisanceParameter", b.name) for b in blocks if b.kind == "Systematic"}
    shared_with_mumu_or_ee = {"MuonID", "MuonIso", "MuonTrigger", "MuonScale", "ElectronID", "ElectronReco", "ElectronScale"}
    assert not nps & shared_with_mumu_or_ee and {f"{n}_tautau" for n in shared_with_mumu_or_ee} <= nps
    assert {"Lumi", "Pileup", "L1Prefiring", "PDF", "QCDScale", "PS_ISR", "PS_FSR"} <= nps and "XS_TTbar" not in nps


def test_mumu_emu_region_is_dropped():
    blocks, info = _adapt("mumu")
    assert [b.name for b in blocks if b.kind == "Region"] == ["mumu_SR"]
    assert info["dropped_regions"]["regions"] == ["mumu_CRemu"]
    names = {(b.kind, b.name) for b in blocks}
    assert not names & {("Sample", "WJets"), ("Systematic", "XS_WJets"), ("Systematic", "ElectronEff_mumu")}
    for b in blocks:
        assert "mumu_CRemu" not in b.opts.get("Regions", "") and "WJets" not in trexcfg.listopt(b.opts.get("Samples", ""))
    # as delivered it is a validation region: not part of the likelihood either way
    delivered = trexcfg.read(repo_path(manifest()["channels"]["mumu"]["config"]))
    assert [b.opts["Type"] for b in delivered if b.kind == "Region" and b.name == "mumu_CRemu"] == ["VALIDATION"]


def test_ttbar_variation_replaces_the_free_normalisation():
    blocks, _ = _adapt("tautau", normfactor_to_overall={"mu_ttbar": {"name": "XS_TTbar", "rel": 0.06}})
    assert not [b for b in blocks if b.kind == "NormFactor" and b.name == "mu_ttbar"]
    xs = [b for b in blocks if b.kind == "Systematic" and b.name == "XS_TTbar"]
    assert len(xs) == 1 and xs[0].opts["OverallUp"] == "0.06" and "TTbar_tDMnone" in xs[0].opts["Samples"]


def test_split_shape_norm():
    blocks = [trexcfg.Block("Region", "r"),
              trexcfg.Block("Systematic", "Pileup", {"Type": "HISTO", "Samples": "s", "Category": "Pileup"}),
              trexcfg.Block("Systematic", "Lumi", {"Type": "OVERALL", "OverallUp": "0.012", "OverallDown": "-0.012"})]
    done = trexcfg.split_shape_norm("ee", blocks, True)
    assert done == ["Pileup"]
    norm, shape = blocks[1], blocks[2]
    assert norm.opts["DropShapeIn"] == "r" and "NuisanceParameter" not in norm.opts
    assert shape.opts["DropNorm"] == "r" and shape.opts["NuisanceParameter"] == "Pileup_eeShape"
    assert "shape_" not in shape.opts["NuisanceParameter"]          # reserved by TRExFitter


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn()
            print("ok", name)
