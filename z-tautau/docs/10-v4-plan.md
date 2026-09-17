# 10 — v4: the four-channel measurement (τhτh + μτh + eτh + eμ)

v4 extends the τhτh measurement (v3, docs 00–09) by the three lepton channels of the CMS 13 TeV
measurement (arXiv:1801.03535, 2.3 fb⁻¹, 2018): **μτh**, **eτh** and **eμ**. The **μμ channel is deliberately
absent** and every channel vetoes a second muon or electron, so v4 is orthogonal to the Z→μμ and Z→ee
selections of the other two groups. As in the paper, the τh identification scale factors and the τh energy
scale are constrained by the data in the combined fit, and every systematic source of the paper's Table 2
is either implemented from an official CMS correction, measured in situ, or replaced by a stated estimate
(section 9).

Status of every item is kept in this file; the numbers of the result are in `output_v4/RESULTS.md`.

## 1. Data streams and channels

| channel | stream (Run2016G+H, NanoAODv9) | trigger | region names |
|---|---|---|---|
| τhτh | Tau (records 30532, 30565) | di-τ 35 GeV (`config.DITAU_TRIGGERS`) | `tautau_SR0/1/2` (v3 BDT categories, unchanged) |
| μτh | SingleMuon (30530, 30563) | `HLT_IsoMu24 ‖ HLT_IsoTkMu24` | `mutau_SR` |
| eτh | SingleElectron (30529, 30562; EOS only, no dCache mirror) | `HLT_Ele27_WPTight_Gsf` | `etau_SR` |
| eμ | MuonEG (30528, 30561) | `HLT_Mu8_TrkIsoVVL_Ele23_…_IsoVL(_DZ)` ‖ `HLT_Mu23_TrkIsoVVL_Ele12_…_IsoVL(_DZ)` | `emu_SR`, `emu_CRtt` |

Luminosity 16393.381 pb⁻¹ (normtag, record 1059) for every stream; the same golden JSON and run range.
The DZ versions of the eμ paths exist only in Run2016H (the non-DZ ones were prescaled there): the OR of
all four is used in data and in the simulation, and the efficiency of the OR is measured in situ (section 7).

Same-stream side samples: SingleMuon events with a muon and an electron (`skim_cat` bit EM_MU) and
SingleElectron events with an electron and a muon (bit EM_EL) are kept for the trigger-efficiency
measurements only; they are never fitted.

## 2. Objects

| object | definition | source of the scale factor |
|---|---|---|
| muon (μτh) | pT > 26, \|η\| < 2.4, `tightId`, \|dxy\| < 0.045, \|dz\| < 0.2, I_rel(ΔR 0.4, Δβ) < 0.15, matched to an IsoMu24/IsoTkMu24 object (ΔR < 0.3) | Muon POG UL16postVFP `muon_Z.json`: `NUM_TightID_DEN_TrackerMuons`, `NUM_TightRelIso_DEN_TightIDandIPCut`, `NUM_IsoMu24_or_IsoTkMu24_DEN_CutBasedIdTight_and_PFIsoTight` |
| electron (eτh) | pT > 29, \|η\| < 2.1, not 1.4442 < \|η_SC\| < 1.566, `mvaFall17V2noIso_WP90`, \|dxy\| < 0.045, \|dz\| < 0.2, I_rel(0.3) < 0.10, conversion veto, ≤ 1 lost hit, matched to an Ele27_WPTight object | EGM `electron.json`: `RecoAbove20`, `wp90noiso`; trigger: in situ (section 7) |
| eμ leptons | μ pT > 10 (tight ID, iso < 0.15), e pT > 13 (WP90, iso < 0.15); the leading lepton pT > 24 (the 23 GeV leg), the other above its 8/12 GeV leg; \|η\| < 2.4 / 2.5 | same ID/iso SFs; cross trigger in situ |
| τh (μτh, eτh) | pT > 30, \|η\| < 2.3, DM 0/1/10/11, \|dz\| < 0.2, VSjet ≥ Tight (`config.NOMINAL_WP`), μτh: VSe ≥ VVLoose & VSmu ≥ Tight, eτh: VSe ≥ Tight & VSmu ≥ VLoose, ΔR(ℓ, τh) > 0.5, opposite charge | ID SF: **fitted** (section 8); e→τh / μ→τh: TAU `tau.json` `DeepTau2017v2p1VSe/VSmu` at the channel's working point; energy scale: `tau_energy_scale` central value, ±3 % prior |
| τh (τhτh) | v3 (pT > 40, \|η\| < 2.1, Tight, VVLoose VSe, VLoose VSmu, both trigger-matched) | trigger SF TauPOG (v3); ID SF fitted |
| jets | anti-kT 0.4 CHS, pT > 30, \|η\| < 4.7, tight ID, ΔR > 0.5 to the selected leptons/τh; b jets: pT > 20, \|η\| < 2.4, DeepJet > 0.2489 (medium) | b-tag: BTV `btagging.json` `deepJet_comb` (b, c) / `deepJet_incl` (light) with efficiencies from the tt̄ simulation |
| MET | PF MET with the τh / lepton energy-scale shifts propagated | unclustered energy from NanoAOD |

## 3. Event selection per channel

Common: MET filters, ≥ 1 good vertex, the veto of section 4.

* **μτh**: exactly the μ and τh above, OS, **m_T(μ, MET) < 40 GeV**. AR: τh VVVLoose-not-Tight. Same-sign,
  m_T > 70 (W DR) and anti-isolated-lepton (0.15 < I_rel < 0.5) sidebands as in section 5.
* **eτh**: the same with the electron.
* **eμ**: e and μ above, OS, **D_ζ = P_ζ^miss − 0.85 P_ζ^vis > −20 GeV**, no b jet. tt̄ control region
  `emu_CRtt`: D_ζ < −40 GeV and MET > 80 GeV (the paper's CR, b-inclusive).
* **τhτh**: v3.

Fit variable everywhere: the MET-likelihood mass m_ττ of docs/04 (`mass.likelihood_mass`), extended with the
leptonic decay phase space (dN/dx = (5 − 9x² + 4x³)/3 on x ∈ [0, 1]) for the e/μ legs. Binning
`config.FIT_BINS_LTAU/EMU`; the τhτh categories keep `config.FIT_BINS`.

## 4. Orthogonality to Z→μμ and Z→ee (and between the channels)

Every channel rejects the event if it contains **any additional muon** (`looseId`, pT > 10, \|η\| < 2.4,
I_rel < 0.3, IP cuts) **or electron** (WP90 noIso, pT > 10, \|η\| < 2.5, I_rel < 0.3, conversion veto)
beyond the channel's own leptons (`config.VETO_MU_V4`, `VETO_EL_V4`):

* the z-mumu signal region needs two tight muons with pT > 20 and I_rel < 0.15 — any such second muon
  fails our veto, so no μτh / eμ / eτh / τhτh event is in the z-mumu SR;
* the z-ee signal region needs two medium cut-based electrons with pT > 20 — a medium cut-based electron
  is isolated (the cut-based ID contains an isolation cut) and passes WP90 in > 99 % of the cases; the residual
  overlap (cut-based-medium but MVA-fail electrons) is checked on the eτh skim and reported in
  `output_v4/data/orthogonality.json`;
* μτh has no electron, eτh no muon, eμ exactly one of each and τhτh none: the four channels are disjoint
  (τhτh additionally uses a different stream; an event in two streams would need two triggers *and* pass
  two disjoint selections, which cannot happen).

**Known overlap with a z-mumu *control* region:** z-mumu's tt̄ control region `mumu_CRemu` (SingleMuon,
IsoMu24, one tight muon pT > 26 + one electron, 60 < m_eμ < 120) shares events with our eμ signal region
(different stream, but the same collisions). A joint fit of z-mumu and v4 must drop `mumu_CRemu` (its tt̄
constraint is irrelevant for Z→μμ); the user asked for orthogonality to the μμ/ee *signal* selections, which
holds.

## 5. Backgrounds

| channel | jet→τh (or multijet) | other |
|---|---|---|
| τhτh | v3 fake-factor method (same-sign DR, era × DM × N_jets × pT, closure corrections, C_OS/SS, closure NPs) + the paper's **30 % on the W+jets/tt̄ part** of the AR (from simulation, `FakeWTT_tautau`) | v3 simulation; τh ID SF fitted |
| μτh, eτh | **fake factors as in the paper**: FF(DM, pT, N_jets) measured in a multijet DR (same sign, lepton isolated) and a W+jets DR (OS, m_T > 70, no b jet), tt̄ FF from the tt̄ simulation; applied to the AR (τh VVVLoose-not-Tight) weighted with the fractions R_p(m_T, N_jets) of multijet / W / tt̄ jet fakes in the AR (W and tt̄ from simulation, multijet = data − simulation); OS/SS extrapolation of the multijet FF from the anti-isolated-lepton sideband; m_T extrapolation of the W FF from the W simulation; genuine τh and ℓ→τh subtracted from every region. DY jet fakes use the W FF (paper). Same-sign validation region. | Z→ττ (genuine τh), Z→μμ/ee (μ/e→τh; VSmu/VSe SFs), tt̄, single top, dibosons with a genuine or lepton-faked τh from simulation |
| eμ | multijet = same-sign data − simulation × OS/SS from the isolation sideband SB1 (both I_rel < 0.6, ≥ 1 above 0.15), measured vs pT and ΔR(e, μ); SB2 (≥ 1 lepton I_rel > 0.6) gives the systematic | tt̄ (normalisation free, `mu_ttbar`, constrained by `emu_CRtt`; top-pT reweighting shape), single top, WW/WZ/ZZ, Z→ℓℓ, W+jets from simulation |

### 5.1 What the fake-factor measurement showed (findings while building v4)

* **The same-sign region with an isolated lepton is about half W+jets**, not multijet. A "multijet" FF measured
  there with only the genuine-τh simulation subtracted is the FF of that mixture; used together with a
  separate W fraction it double counts the W part (16 % over-prediction of the same-sign data). v4 subtracts the
  simulated W+jets / Drell-Yan / top **jet fakes** from the same-sign DR and from the anti-isolated sidebands as
  well, so `FF_qcd` is the FF of pure multijet events and the fractions R_p are consistent with it.
* **The W+jets FF depends on the relative charge of the lepton and the τh candidate**: the jet recoiling against
  a W is a quark jet whose charge is anti-correlated with the W charge, so opposite-sign W fakes are quark-like
  and same-sign ones gluon-like. Measured in data at m_T > 70 GeV: FF_W(OS) ≈ 0.08, FF_W(SS) ≈ 0.043 in μτh.
  The same-sign validation therefore uses a same-sign W FF (`w_ss` table); with it the same-sign data are
  reproduced to 0.97 ± 0.01 (μτh) and 1.03 ± 0.02 (eτh), and those numbers are the `FakeClosure_<ch>` priors.
* The e μ isolation sidebands are 0.15–0.5 (both leptons, ≥ 1 above 0.15) and 0.3–0.5 because the skims keep
  leptons with I_rel < 0.5; the OS/SS ratio falls from 2.9 (ΔR < 1) to 1.65 (ΔR > 4) and the two sidebands agree
  within 2 %.

## 6. Signal definition and parameter of interest

v4 measures **σ(pp → Z/γ* → ττ, 60 < m_LHE < 120 GeV)** with all τ decays as the signal
(`DYtautau`, one `mu_Z` for all channels). Z/γ*→ττ outside the window (`DYtautau_out`, mostly
m > 120 GeV) is a theory-normalised background (5 %). The theory variations (QCD scales, PDF, PS) keep
σ(60–120) fixed and vary only the acceptance × efficiency per channel. v3's τhτh visible fiducial volume is
still reported per channel (σ_fid = μ × σ_fid^pred) but is no longer the definition of the POI. The theory
reference is σ(60–120) of the aMC@NLO sample normalised to 6077.22 pb (m > 50).

## 7. Corrections, and what is measured in situ

| correction | applied to | source |
|---|---|---|
| pileup, L1 prefiring | all MC | v3 (lumi-by-LS profile, NanoAOD weights) |
| muon ID, iso, IsoMu24 trigger | genuine muons | Muon POG json (section 2), ±1σ each |
| electron reco, ID | genuine electrons | EGM json, ±1σ each |
| Ele27 trigger efficiency | eτh, EM_EL | **in situ**: eμ events of SingleMuon (IsoMu24, μ pT > 26): ε(Ele27 ‖ electron) in data and simulation vs pT, η → SF ± stat ⊕ 2 % above 35 GeV (`ElectronTrigger`); in the 29–35 GeV turn-on bin the flat 2 % is replaced by the step to the next pT bin of the measured table (4 % / 2 % / 13 % for \|η_SC\| < 0.8 / < 1.479 / < 2.1) as a separate parameter `ElectronTrigger_lowpt` |
| eμ cross-trigger efficiency | eμ | **in situ**: the electron legs from SingleMuon eμ events, the muon legs from SingleElectron eμ events; SF(pT_e) × SF(pT_μ), varied by the statistical error of **both** table bins and the 2 % of the paper **once per event** (2.7–2.9 % in total). Applying the 2 % per leg, as the first version did, doubles the prior and moves μ_Z by one standard deviation |
| di-τ trigger | τhτh | TauPOG (v3) |
| τh ID SF per DM | genuine τh, every channel | **free NormFactors** `TauIDSF_DM0/1/10/11` (nominal = POG value) |
| τh energy scale per DM | genuine τh, every channel | POG central value; **±3 % prior**, constrained in situ |
| e→τh, μ→τh | prompt e/μ reconstructed as τh | TAU json at the channel's VSe/VSmu working point |
| electron energy scale/resolution | eτh, eμ | NanoAOD `Electron_dEscaleUp/Down`, `dEsigmaUp/Down` |
| muon momentum scale | μτh, eμ | ±0.2 % (Rochester-corrections size, estimate) |
| b-tagging | eμ b veto, W DR b veto | BTV json + MC efficiencies; if the efficiency map is not available, ±1 % (signal) / ±6 % (tt̄) estimate |
| jet energy scale | N_jets (FF binning), b jets, MET | JME json `Total` uncertainty on the jet pT, propagated to N_jets, b jets and MET (type-1) |
| top pT reweighting | tt̄ | applied (CMS 2016 fit), shape NP = no / twice |
| MET unclustered | all MC | v3 |

## 8. Fit model

Regions: `tautau_SR0/1/2` (v3, SR0 above 110 GeV), `mutau_SR_dm{0,1,10,11}`, `etau_SR_dm{0,1,10,11}` (one region
per τh decay mode: with a single m_ττ distribution per channel the per-decay-mode scale factors were degenerate and
the first combined fit drove two of them to the boundary; per decay mode each region measures SF(DM) × μ_Z and the
eμ channel fixes μ_Z), `emu_SR`, `emu_CRtt`.
Samples per region: `Data`, the simulation split by fit sample **and by the decay modes of the genuine τh
legs** (`<sample>_tDM<key>`, key ∈ {none, 0, 1, 10, 11, 0_0, 0_1, …}), `Fakes` (data-driven).
NormFactors: `mu_Z` (POI, on every `DYtautau_*`), `TauIDSF_DM<d>` on every template with one genuine τh of
that decay mode, `Expression` products (`TauIDSF_DM0*TauIDSF_DM1`, `TauIDSF_DM1*TauIDSF_DM1`, …) on the
τhτh templates with two genuine legs, `mu_ttbar` on tt̄. Nuisance parameters: section 9. The τhτh templates
carry **no** POG ID SF (the NormFactor is the SF); their v3 `TauID_DM*` NPs are gone.

Outputs: `fit/ztautau.config`, `fit/fitinputs/ztautau.root`, `fit/results/ztautau_fit_result.json` — the
canonical names of `fitting/CONVENTIONS.md`, because this is the only measurement of the channel. The τhτh
Data and fake templates come from `fit/fitinputs/tautau_base.root` (`run_tautau_base.py`).

**Regions that exist in the file but are not in the fit.** `mutau_SRlo/hi_dm*` and `etau_SRlo/hi_dm*` hold the
same events as `mutau_SR_dm*` / `etau_SR_dm*`, split at pT(τh) = 40 GeV. They are filled by step 4 (tagged
`"ptsplit"` in `meta["region_sets"]`) and fitted only by the cross-check job `ztautau_ptsplit`, which gives
the 30–40 GeV part its own `TauIDSF_DM*_lowpt` scale factors. That is the test of the one model assumption
the τhτh / (ℓτh)² lever rests on: the τhτh legs are above 40 GeV, most of the ℓτh signal is not, and a
scale factor that varies between 30 and 40 GeV would show up in μ_Z as ≈ 1–2 times that variation.

**Cross-check jobs of step 5** (`run_all.py --help`, `REVIEW_v4_RESPONSE.md`): `ztautau_<ch>` (one channel,
scale factors free: the workspaces of the MultiFit), `ztautau_<ch>_fixedid` (POG scale factors fixed),
`ztautau_taulep` (the three τ channels without eμ), `ztautau_ptsplit`, `ztautau_emutrig2x`.

## 9. Systematic uncertainties: CMS Table 2 ↔ v4

| paper source (Table 2) | v4 nuisance parameter(s) | status |
|---|---|---|
| Integrated luminosity 2.3 % | `Lumi` 1.2 % (record 1059) | as v3 |
| Hadronic τ ID and trigger | `TauIDSF_DM*` free; `TauTrigger_DM*` (τhτh) | fitted / TauPOG |
| τh ES (3 % prior) | `TauES_DM*` ±3 % prior, m_ττ recomputed, migration in/out of the selection | implemented |
| Rate of e misidentified as τh | `TauFakeEle` (VSe at Tight in eτh, VVLoose elsewhere) | TAU json |
| Rate of μ misidentified as τh | `TauFakeMu` (VSmu Tight in μτh, VLoose elsewhere) | TAU json |
| Electron ID and trigger | `ElectronReco`, `ElectronID`, `ElectronTrigger`, `ElectronTrigger_lowpt` | EGM json / in situ (plateau and turn-on separately) |
| e ES | `ElectronScale`, `ElectronRes` | NanoAOD variations |
| Muon ID and trigger | `MuonID`, `MuonIso`, `MuonTrigger` | Muon POG json |
| μ ES | `MuonScale` ±0.2 % | estimate |
| MET response and resolution | `MET_Unclustered`, `JES` | implemented / JME json |
| Norm. Z→ee, μμ (unconstrained) | `XS_DYll` 5 % | fixed prior (the ee/μμ channels of this project measure it to 2 %) |
| Norm. and shape of false τh | `FakeStat_*` (γ), `FakeOSSS_<ch>`, `FakeWmT_<ch>`, `FakeFrac_<ch>`, `FakeClosure_<ch>`, `FakeMCSub_<ch>`; τhτh: v3 set + `FakeWTT_tautau` 30 % | implemented |
| Norm. and shape of multijet (eμ) 20 % | `QCDOSSS_emu` (SB1−SB2 ⊕ stat), `QCDStat_emu` (γ) | implemented |
| Norm. tt̄ 7 % | `mu_ttbar` free + `emu_CRtt` | implemented |
| Shape tt̄ | `TopPt` | implemented |
| Norm. SM H 30 % | — | not simulated; ≤ 0.2 % of the signal (H→ττ σ×B = 3.3 pb vs 1900 pb), noted |
| Norm. single top 15 %, diboson 15 %, W+jets 15 % | `XS_SingleTop`, `XS_WW/WZ/ZZ`, `XS_WJets` 10 % | `fitting/CONVENTIONS.md` values (shared with ee/μμ); the paper's 15 % noted |
| PDF, scale dependence, UE and PS | `PDF`, `QCDScale`, `PS_ISR`, `PS_FSR` (fixed σ(60–120), per channel A×ε) | as v3, extended to all channels |
| — (not in the paper) | `L1Prefiring`, `Pileup`, `BTag`, `MCStatNorm_WJets_tautau`, γ | as v3 |

## 10. Skims, ntuples, scripts

* `scripts/step1_skim.py --v4` → `$BND_TAUTAU_CACHE/skims_v4/<sample>/` (`skim.v4_categories`: bitmask
  `skim_cat` TT=1, MT=2, ET=4, EM=8, EM_MU=16, EM_EL=32; a data stream keeps only its own bits, the
  simulation all). Loose objects (muon medium ID, I_rel < 0.5; electron WP90, I_rel < 0.5; τh pT > 18,
  VVVLoose) so every sideband is inside. Trigger objects of ids 11/13/15. Per-file provenance with the
  category counts; measured rates: SingleMuon 2.5 %, MuonEG 1.4 %, SingleElectron 1.6 %, DY 2.7 %, tt̄ 25 %.
  The Tau stream is **not** re-skimmed: τhτh reuses `skims_v1` / `ntuples_v1`.
* `scripts/step2_ntuples_lepton.py` → `ntuples_v4/<sample>_<channel>.root` for `mutau`, `etau`, `emu`
  (+ `_trig` side samples), one flat tree with the pair, the region flags, MET, jets, masses and weights.
* `scripts/step3d_fakes_lepton.py` (lepton-channel fake factors, fractions, OS/SS, closure; `output/data/fakes_<ch>.json`),
  `scripts/step3c_trigger.py` (in-situ trigger efficiencies; `external/trigger_insitu_v4.json`),
  `scripts/step4_histograms.py` (all regions, all variations, DM-split templates → `fit/fitinputs/`),
  `scripts/step4b_export_channels.py` (one file per channel), `scripts/step5_fit.py` (TRExFitter),
  `scripts/step5b_multifit.py` (the MultiFit of the four channel workspaces), `scripts/step6_report.py`,
  `run_all.py`. The τhτh base (`fit/fitinputs/tautau_base.root`) comes from `run_tautau_base.py`
  (`scripts/step1_skim.py`, `step2_ntuples_tautau.py`, `step3_fakefactors.py`, `step3b_bdt.py`,
  `step4a_tautau_base.py`).

## 11. Limitations (v4)

* SM H→ττ not simulated (≤ 0.2 %). EWK Z→ττ absent (0.3 %, as v3).
* The eμ channel shares events with z-mumu's `mumu_CRemu` (section 4).
* **The τh ID scale factor is not flat in pT, and the measurement depends on that.** The fit uses one
  scale factor per decay mode, while the τhτh legs are above 40 GeV and most of the ℓτh signal is not.
  `ztautau_ptsplit` gives the 30–40 GeV ℓτh regions their own scale factors and finds them **higher** than
  the > 40 GeV ones by 6 % (DM0), 0 % (DM1), 6 % (DM10) and 14 % (DM11) — each 0.5–1 σ on its own, but
  coherent in sign for three of the four decay modes. Freeing them moves the result from
  μ_Z = 1.019 ± 0.037 to **0.954 ± 0.034**, i.e. by 1.8 times the total uncertainty, and it removes most of
  the disagreement between the eμ channel (0.960) and the τ channels. So the single-scale-factor
  assumption is not an academic caveat: it is worth more than every experimental systematic in the table.
  It is not adopted as the nominal model because the four ratios are individually compatible with one and
  the split doubles the number of free scale factors on the same data; but any use of this measurement
  should carry the 6 % spread between the two fits. The TauPOG pT-binned prescription, which is 7 % below
  the DM-binned one above 40 GeV (`docs/07`), points the same way.
* With the scale factors free, the eμ channel (no τh) and the three τ channels measure μ_Z in different ways
  and the combined number is their compromise; both are quoted (`docs/11-combination-inputs.md` §7).
* `MET_Unclustered` changes the ℓτh *yield* by −3.8 % / +2.9 % through the m_T < 40 GeV cut and the MET term
  of the likelihood mass. It is one nuisance parameter for both the shape and that normalisation, on purpose
  (`REVIEW_v4_RESPONSE.md` §4), and the fit constrains it to ≈ 0.3 of its prior.
* Muon momentum scale, b-tag (if no efficiency map) and W+jets/single-top/diboson normalisations are
  estimates rather than measurements; they are marked as such in the results tables.

## 12. Result

The full tables are `output/RESULTS.md` (regenerated by `scripts/step6_report.py`); the headline is
injected here by `scripts/update_docs.py`:

<!-- RESULT:BEGIN -->
**Result** (four channels, τh ID scale factors and energy scale fitted in situ)

```
σ(pp → Z/γ* → ττ, 60 < m < 120 GeV) = 1981 +73 −70 pb   (stat ±9;  prediction 1945 pb: aMC@NLO acceptance, NNLO normalisation)
μ_Z = 1.019 +0.037 −0.036   (stat ±0.005, syst ±0.036),   expected ±0.036,   GoF p = 0.15,   μ_tt̄ = 1.11 ± 0.04
τh ID SF (Tight):  DM0 0.989 ± 0.040 (POG 0.90 ± 0.13)  DM1 0.963 ± 0.034 (POG 0.89 ± 0.05)  DM10 0.896 ± 0.035 (POG 0.94 ± 0.15)  DM11 0.795 ± 0.050 (POG 0.81 ± 0.15)
τh energy scale:   DM0 -0.6 ± 0.9 %  DM1 -0.2 ± 0.6 %  DM10 +0.6 ± 1.0 %  DM11 +3.2 ± 2.2 %
μ_Z of the two sub-measurements the fit combines:  e mu alone 0.960 +0.042 −0.040;  τ channels alone (SF free) 1.203 +0.091 −0.084   (2.6 σ apart)
per channel alone (POG SFs fixed):  tautau 1.043 +0.079 −0.072  mutau 1.048 +0.055 −0.051  etau 1.016 +0.070 −0.065  emu 0.960 +0.042 −0.040
```

Largest grouped impacts on μ_Z: NormFactors 3.1 %, Emu trigger 2.5 %, Gammas 1.4 %, Tau trigger 1.2 %, Fakes 1.2 %, Electron energy 1.2 %; data statistics 0.5 %. The categories overlap, so their quadrature sum exceeds the MINOS total by a factor 1.46 (`output/RESULTS.md`).

Ranking: mu_ttbar 3.1 %, EmuTrigger 2.6 %, TauIDSF_DM1 2.2 %, TauIDSF_DM10 2.1 %, TauIDSF_DM0 1.9 %, TauIDSF_DM11 1.3 %.
<!-- RESULT:END -->

## 13. Channel readiness and the per-channel exports for a combination

All four channels are fitted and exported; nothing of the v4 chain is pending. What "ready" means per channel:

| channel | regions | data stream | fake estimate | status |
|---|---|---|---|---|
| τhτh | `tautau_SR0/1/2` | Tau | v3 same-sign fake factor | ready (v3 chain, templates rebuilt with the v4 signal definition and free τh ID SFs) |
| μτh | `mutau_SR_dm0/1/10/11` | SingleMuon | per-process FF (multijet / W / tt̄) | ready |
| eτh | `etau_SR_dm0/1/10/11` | SingleElectron | per-process FF | ready |
| eμ | `emu_SR`, `emu_CRtt` | MuonEG | same-sign × OS/SS | ready |
| μμ, ee | — | — | — | deliberately absent (the other groups' channels; v4 is orthogonal to them) |

`scripts/step4b_export_channels_v4.py` splits the combined fit inputs into **one file per channel**
(`fit_v4/fitinputs/ztautau_v4_<channel>.root` + `.meta.json`, histograms copied bit for bit, the systematic
registry restricted to the channel), `step5_fit_v4.py --channels <ch> --job ztautau_v4_<ch>` builds the channel's
config, workspace (`fit_v4/results/ztautau_v4_<ch>/RooStats/ztautau_v4_<ch>_combined_ztautau_v4_<ch>_model.root`)
and fit result from that file with the τh ID scale factors **free**, and `fit_v4/comb_v4.config` is the TRExFitter
MultiFit of the four workspaces (`trex-fitter mwf comb_v4.config` from `fit_v4/`). Nuisance parameters are
correlated by name across the channel workspaces (`TauIDSF_DM*`, `TauES_DM*`, `TauFakeEle/Mu`, `Lumi`, `Pileup`,
`L1Prefiring`, `Muon*`, `Electron*`, `EmuTrigger`, `BTag`, `JES`, `MET_Unclustered`, `TopPt`, `QCDScale`, `PDF`,
`PS_*`, `XS_*`, `mu_ttbar`); channel-specific ones carry the channel suffix (`Fake*_<ch>`, `QCDOSSS_emu`,
`MCStatNorm_WJets_tautau`). The MultiFit of the four workspaces therefore reproduces the single-file fit of step 5
(checked, section 12 / `output_v4/RESULTS.md`). The `_fixedid` jobs are the per-channel cross-checks with the
POG scale factors held fixed (a single ℓτh channel alone cannot separate μ_Z from its scale factors).

Adding the other groups' channels: append `Fit:` blocks to `comb_v4.config` pointing at `../../z-mumu/fit/zmumu.config`
and the z-ee config once it follows `fitting/CONVENTIONS.md`; drop z-mumu's `mumu_CRemu` (it shares events with
`emu_SR`); the shared names above then correlate automatically, and every channel's σ^pred(60–120) must be its own
(`CONVENTIONS.md` §6). The per-channel signal templates are split by τh decay mode (`DYtautau_tDM<key>`), so a
combination must scale all `DYtautau_*` templates with `mu_Z` (the exported configs do).
