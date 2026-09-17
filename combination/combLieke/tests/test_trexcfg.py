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
    assert "mu_Z" not in nf and nf["mu_Z_tautau"].opts["Samples"] == "DYtautau"
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
    sr2 = [b for b in blocks if b.kind == "Region" and b.name == "tautau_SR2"][0]
    assert sr2.opts["DropBins"] == "1"
    try:
        _adapt("tautau", drop_empty_bins={"tautau_SR2": [5]})
    except ValueError as err:
        assert "not empty" in str(err)
    else:
        raise AssertionError("a non-empty bin was dropped")


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
