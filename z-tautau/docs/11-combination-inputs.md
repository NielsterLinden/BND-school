# 11 — What the combination gets from Z → ττ, and how to use it

Everything in this document is committed in the repository, so a combination can be built from a fresh
checkout without running anything here. It is the channel-side counterpart of `fitting/CONVENTIONS.md`
and of `z-mumu/docs/15-combination-inputs.md`.

Read section 7 before quoting a number: the four-channel result is a likelihood compromise between two
sub-measurements that do not agree, and which of them a combination should use is a choice the
combination has to make consciously.

---

## 1. Deliverables

| what | path | note |
|---|---|---|
| fit inputs, all channels | `fit/fitinputs/ztautau.root` (+ `.root.meta.json`) | `TH1D` with `Sumw2`; 13 fitted regions, plus 16 pT-split copies of the ℓτh ones that are **not** part of the fit (§2) |
| fit inputs, one channel | `fit/fitinputs/ztautau_<ch>.root` (+ `.meta.json`) | `<ch>` = `tautau`, `mutau`, `etau`, `emu`; histograms copied bit for bit out of the file above |
| TRExFitter config, all channels | `fit/ztautau.config` | job `ztautau`, POI `mu_Z` |
| TRExFitter config, one channel | `fit/ztautau_<ch>.config` | τh ID scale factors **free** |
| workspace, all channels | `fit/results/ztautau/RooStats/ztautau_combined_ztautau_model.root` | copy of TRExFitter's `ztautau_allBinsFitRegions_combined_ztautau_model.root` (a region uses `DropBins`) |
| workspace, one channel | `fit/results/ztautau_<ch>/RooStats/ztautau_<ch>_combined_ztautau_<ch>_model.root` | what a MultiFit reads |
| fit results | `fit/results/ztautau_fit_result.json`, `fit/results/ztautau_<ch>_fit_result.json` | POI, NPs, grouped impacts, ranking |
| MultiFit of our four channels | `fit/comb.config`, `fit/results/comb_fit_result.json` | the closure check of section 6 |
| everything in one place | `output/results.json` → `for_combination` | the channel result as a flat block: μ, σ^pred, grouped impacts, ranking |
| the numbers in prose | `output/RESULTS.md`, `handoff.md` | |

Regenerate with `python run_all.py --from 4` (needs `fit/fitinputs/tautau_base.root`, itself built by
`python run_tautau_base.py --from 3`).

## 2. Regions, samples, histograms

Histogram names follow `fitting/CONVENTIONS.md` §1: `<region>__<sample>`, `<region>__Data`,
`<region>__<sample>__<syst>Up|Down`. Every histogram is a `TH1D` of the di-τ mass `m_tt` with the sum of
squared weights stored; all samples of a region share its binning; there are no dots in any name.

| region | channel | what it is |
|---|---|---|
| `tautau_SR0`, `tautau_SR1`, `tautau_SR2` | τhτh | BDT categories; **`tautau_SR0` is fitted only above 110 GeV** (`DropBins` in the config): it is the fake sideband |
| `mutau_SR_dm0/1/10/11` | μτh | one region per τh decay mode |
| `etau_SR_dm0/1/10/11` | eτh | one region per τh decay mode |
| `emu_SR` | eμ | signal region (no τh) |
| `emu_CRtt` | eμ | tt̄ control region (D_ζ < −40, MET > 80), fitted with `mu_ttbar` |

The files also contain `mutau_SRlo/hi_dm*` and `etau_SRlo/hi_dm*`: the **same events** as
`*_SR_dm*`, split at pT(τh) = 40 GeV, used only by the cross-check fit `ztautau_ptsplit` (§7). They are
marked `"ptsplit"` in `meta["region_sets"]`. **A combination must not fit them together with the
`*_SR_dm*` regions — it would double count the data.** Take the regions whose `region_sets` entry is
`"nominal"`, or simply use `fit/ztautau.config`, which already does.

Samples are the conventional names with a suffix that carries the τh decay mode of the *genuine* τh legs:
`<sample>_tDM<key>`, `key` = `none`, a decay mode (`0`, `1`, `10`, `11`) or a pair (`0_10`) in τhτh. The
split exists so that one τh ID scale factor per decay mode can multiply the right templates (and its square
the τhτh ones, via TRExFitter `Expression`). Base samples: `DYtautau` (the signal), `DYtautau_out`, `DYee`,
`DYmumu`, `DYlowmass`, `WJets`, `TTbar`, `SingleTop`, `WW`, `WZ`, `ZZ`, `Fakes`, `Data`.

> **`mu_Z` must scale every `DYtautau_tDM*` template and nothing else.** `DYtautau_out` (the Z/γ*→ττ
> simulation outside 60–120 GeV) is a theory-normalised background with a 5 % normalisation NP.

The histograms are already in events (luminosity applied): never set `Lumi:` in a Job block that reads them.

## 3. The measured quantity and its reference

* Signal: **Z/γ* → ττ with 60 < m_LHE < 120 GeV, all τ decays**, one `mu_Z` for all four channels.
* σ^pred(60–120) = **1944.9 pb** = 6077.22 pb × Σw(LHE ττ, 60 < m_LHE < 120) / Σw of the aMC@NLO
  `DYJetsToLL_M-50` sample. This is *our* reference; z-mumu uses 1953.9 pb and z-ee 1954.1 pb.
* σ̂ = μ̂ × 1944.9 pb. `CONVENTIONS.md` §6: **combine the cross sections, never `mu_Z`** — a single shared
  `mu_Z` across channels would fit one parameter against three different references.
* Luminosity 16393.381 pb⁻¹ (normtag), ±1.2 %, the same number in all channels.
* There is **no external acceptance uncertainty to add** (`acc = {}`). The theory templates (`QCDScale`,
  `PDF`, `PS_ISR`, `PS_FSR`) are renormalised to a constant σ(60–120) per variation member, exactly as
  `CONVENTIONS.md` §3 demands, so they vary A × ε only and are already inside the fit and inside the
  `Signal modelling` group. The normalisation effect of every theory NP on the signal template is ≤ 6 %
  in τhτh (the double 40 GeV cut selects the Z pT tail), ≤ 3 % in ℓτh and 0.7 % in eμ; κ_theory ≡ 1.
  This is the test z-ee failed (`combination/docs/01-inputs.md`); ττ passes it.

## 4. Nuisance parameters: what may be correlated by name

TRExFitter correlates by name. Ours, and what they are:

**Correlate with z-mumu / z-ee (same quantity, same source):** `Lumi` (normtag, 1.2 %), `Pileup`
(Open Data lumi table, 69.2 mb), `L1Prefiring` (NanoAOD weights), `PDF`, `QCDScale`, `PS_ISR`, `PS_FSR`
(the same aMC@NLO DY sample and the same members), `XS_TTbar`, `XS_SingleTop`, `XS_WW`, `XS_WZ`, `XS_ZZ`,
`XS_WJets`, `XS_DYll`, `XS_DYlowmass`.

**Ours alone — never correlate:** everything with a channel suffix (`Fake*_mutau`, `Fake*_etau`,
`Fake*_tautau`, `QCDOSSS_emu`, `MCStatNorm_WJets_tautau`), all `gamma_*`, `TauIDSF_DM*` (free
NormFactors), `TauES_DM*`, `TauTrigger_DM*`, `TauFakeEle`, `TauFakeMu`, `XS_DYtautau_out`, `BTag`, `JES`,
`TopPt`, `MET_Unclustered`, and `mu_ttbar` (see below).

**Decide explicitly, do not let the name decide:**

| NP | ours is | z-mumu / z-ee is | recommendation |
|---|---|---|---|
| `MuonID`, `MuonIso`, `MuonScale` | POG UL2016 json (`ztautau/pog.py`) | z-mumu measures its own with tag-and-probe | **decorrelate** (`DecorrSysts`) unless both groups agree they are the same measurement |
| `MuonTrigger` | POG IsoMu24 json | z-mumu tag-and-probe | **decorrelate** |
| `ElectronID`, `ElectronReco`, `ElectronScale` | EGM POG json | z-ee's `ElectronID`/`ElectronRECO`/`ElectronES` are its own | **decorrelate** (the names already differ from z-ee's, which is enough) |
| `ElectronTrigger`, `ElectronTrigger_lowpt`, `EmuTrigger` | **our own in-situ measurement** (`external/trigger_insitu_v4.json`) | not measured elsewhere | uncorrelated, ρ = 0 |
| `mu_ttbar` | a free normalisation of tt̄, determined by `emu_CRtt` | z-mumu constrains tt̄ with `XS_TTbar` (6 %) | decide: either drop `XS_TTbar` from the correlation or fix `mu_ttbar`; they describe the same yield twice |
| `SigModel_tautau` | madgraph LO vs aMC@NLO, reported not fitted | z-mumu's `SigModel` is powheg vs aMC@NLO | different quantities; ours is not in the fit at all |

`DecorrSysts` / `DecorrSuff` must be set **before** the workspace is built, i.e. in each channel's own Job
block, not in the MultiFit.

## 5. Overlap with the other channels

* **τhτh, μτh, eτh are orthogonal** to z-mumu (≥ 2 tight isolated muons) and z-ee (≥ 2 medium electrons):
  every v4 channel vetoes any additional muon (loose ID, pT > 10, I_rel < 0.3) or electron (MVA noIso 90 %,
  pT > 10, I_rel < 0.3).
* **The eμ channel overlaps z-mumu's tt̄ control region `mumu_CRemu`**, which selects exactly one muon and
  one electron. A joint fit must **drop `mumu_CRemu`** (or drop our `emu_CRtt` and `emu_SR`). This is the
  one hard orthogonality condition of a ττ + μμ combination.
* Data statistics are uncorrelated between the channels once `mumu_CRemu` is dropped.
* Our eμ channel would also constrain the shared `QCDScale` to ≈ 0.6 of its prior through the
  τhτh/ℓτh acceptance ratio, which then tightens the μμ and ee acceptance uncertainties. That is only
  legitimate if the constraint is believed — see §7.

## 6. How the combination reads this channel

The live combination is **`combination/combLieke/`**: one TRExFitter v1.8.0 MultiFit of ee, μμ and ττ with a
common POI (`combination/combLieke/README.md`). It takes each channel's own config and fit inputs and
adapts them with `mf/trexcfg.adapt_channel`, driven by `config/channels.json`. Its ττ entry points at
`z-tautau/fit/ztautau.config` and `z-tautau/fit/fitinputs/ztautau.root` — the paths this measurement now
owns — so it will pick up the four-channel model automatically.

**Three entries of `config/channels.json` must be changed with it, and the fit will not run until they are.**

1. **`"signal": ["DYtautau"]` no longer matches anything.** The signal is split by the decay mode of the
   genuine τh legs, so the samples are `DYtautau_tDMnone`, `DYtautau_tDM0`, `DYtautau_tDM1`,
   `DYtautau_tDM10`, `DYtautau_tDM11`, `DYtautau_tDM0_0`, `DYtautau_tDM0_1`, … The adapter renames the
   signal NormFactor to `mu_Z` and puts the constant `xsref_tautau` on those samples, so it needs all of
   them. The exact list for the current inputs is in the metadata as `signal_samples`
   (`fit/fitinputs/ztautau.root.meta.json`), so the entry can read that key instead of being hard-coded.

2. **`"acceptance"` must become `null`.** It currently reads `signal_prediction.A_unc.*` and
   `signal_prediction.A_mc_stat` and adds `Acc_PDF`, `Acc_AlphaS`, `Acc_QCDScale`, `Acc_PS_ISR`,
   `Acc_PS_FSR`, `AccStat_tautau` as OVERALL parameters on the signal — a 3.7 % acceptance uncertainty
   *outside* the fit. That was right for the τhτh-only measurement, which quoted σ = σ_fid / A with A from
   aMC@NLO. It is **wrong here and would double count**: the four-channel fit measures σ(60–120) directly,
   and its `PDF`, `QCDScale`, `PS_ISR` and `PS_FSR` templates are renormalised per variation member to a
   constant σ(60–120) (§3), so they already vary A × ε inside the likelihood, per region. Those keys are
   therefore absent from the metadata and the adapter will raise a `KeyError` rather than guess; the
   metadata says so explicitly in `acceptance_in_fit` and `acceptance_note`. The right entry is

   ```json
   "acceptance": null,
   "acceptance_note": "none outside the fit: the theory templates are renormalised to a constant
                       sigma(60-120), so PDF/QCDScale/PS_* vary A x epsilon inside the likelihood
                       (z-tautau/docs/11-combination-inputs.md section 3)"
   ```

   which is exactly what the ee entry already carries, for the same reason.

3. **`drop_empty_bins` has to be re-derived.** `tautau_SR2` bin 1 is still empty, but the ℓτh and eμ regions
   are new. The adapter refuses to drop a bin that is not exactly zero, so the safe procedure is to run the
   combination once, read which γ parameters are unconstrained at zero, and list those bins.

Two further decisions are the combination's to make, not ours (§4 has the details): whether our POG-based
`MuonID/Iso/Trigger/Scale` and `ElectronID/Reco/Scale` are the same quantity as z-mumu's tag-and-probe and
z-ee's own measurements — by default the MultiFit correlates them **by name**, which is a claim, not a
default — and what to do about tt̄, which we normalise with a free `mu_ttbar` from `emu_CRtt` while z-mumu
constrains it with `XS_TTbar`. Both channels describing the same tt̄ yield twice is not harmless.

`σ_reference_pb` (1944.88 pb) and `sigma_reference_source` are unchanged and still correct.

### Our own MultiFit, as a closure test

`fit/comb.config` is the MultiFit of our four per-channel workspaces and reproduces the single-file fit;
that is the check that the per-channel exports are usable at all. To add the other groups, append `Fit:`
blocks pointing at `../../z-mumu/fit/zmumu.config` and the z-ee config and run `trex-fitter mwf comb.config`
from `z-tautau/fit/`. Remember §3: a single shared `mu_Z` across channels with different σ^pred fits one
parameter against three references, which is precisely why `combLieke` adds the constant `xsref_<channel>`.

### A covariance / BLUE combination, if one is ever wanted again

`output/results.json` carries a ready-made `for_combination.channel_result` block with the fields the
earlier BLUE combination used (`combination/comb/inputs.ChannelResult`, removed in commit `b6aeeb2`):
`mu`, its MINOS and stat errors, `sigma_pred`, the grouped impacts, the ranking, `acc = {}` (§3),
`sigmodel = 0`, `residual = 0` and `gof_p`. One caveat, and it is the reason the block also carries
`groups_rescaled`:

> **The grouped impacts overlap and their quadrature sum over-shoots the MINOS total.** `mu_ttbar`, the eμ
> trigger efficiency and `mu_Z` form one chain (the control region fixes μ_tt̄ × ε, the signal region
> μ_Z × ε), so the categories `NormFactors` and `Emu trigger` contain the same degeneracy. Use
> `groups_rescaled` (every row multiplied by `groups_rescale_factor`, so the quadrature sum equals the
> MINOS total) with `residual = 0`, or the MINOS total directly. The raw rows inflate our uncertainty by
> more than half and silently change our weight.

Categories and the recommended ρ with μμ/ee are in `for_combination.category_rho_recommended` and in the
grouped-impact table of `output/RESULTS.md`. v4 adds categories the old model did not know
(`NormFactors`, `Emu trigger`, `Electron trigger`, `Electron energy`, `Jets`, `b tagging`,
`Background modelling`, `Tau ID (fitted)`); all are ρ = 0 against μμ and ee.

## 7. What a combiner must know before quoting our number

1. **The four-channel result is a compromise between two sub-measurements.** The eμ channel measures
   `mu_Z` directly (no τh, no ID scale factor). The three τ channels with free scale factors measure it
   only through the ratio τhτh / (ℓτh)², because the scale factor enters the τhτh templates squared and the
   ℓτh templates once. The two are quoted separately in `output/RESULTS.md` with their difference. Which
   number a combination uses changes the ττ line by ±10 %; the combined one is the likelihood answer and
   the default, but it is not independent of the eμ trigger efficiency.
2. **The eμ trigger prior sets the scale of the eμ sub-measurement.** It is our own in-situ measurement
   plus the 2 % of CMS arXiv:1801.03535, applied **once per event** (before the review it was applied once
   per leg, which doubled it and moved `mu_Z` by a full standard deviation). `ztautau_emutrig2x` is the
   fit with the old, doubled prior and is quoted in `RESULTS.md` as the size of that choice.
3. **The τ-channel lever assumes one τh ID scale factor per decay mode over pT(τh) > 30 GeV.** The
   cross-check fit `ztautau_ptsplit` gives the ℓτh regions below 40 GeV their own scale factors and
   measures the difference; its table is in `RESULTS.md`.
4. `mu_ttbar` is free and comes out above 1; it is the eμ trigger–tt̄ degeneracy, not a statement about
   the tt̄ cross section (§4).
5. The τh ID scale factors and energy scales are *results* of this channel (`RESULTS.md`), measured in
   situ at 4–5 % and 0.6–2.2 %. In a full combination where μμ and ee fix the cross section, that is what
   this channel mainly adds.
