# Handoff – Combination

## Team

- Lieke Gijsen
- Caspar van Eck

## What we did

We combined the three channels, `z-ee` (Z → e⁺e⁻), `z-mumu` (Z → μ⁺μ⁻) and `z-tautau` (Z → ττ: τhτh + μτh + eτh +
eμ), in **one TRExFitter v1.8.0 MultiFit**: one signal strength for all channels, and every shared nuisance
parameter fitted jointly. The fits run on HTCondor (17 Sep 2026).

> **σ(pp → Z/γ\* → ℓℓ, 60 < m_ℓℓ < 120 GeV) = 1949 ⁺³⁰₋₃₀ pb**
> = 1949.3 ± 0.5 (stat) ⁺³⁰·²₋₂₉.₆ (syst, of which luminosity 21.9 pb and acceptance 12.8 pb) pb,
> per lepton flavour, assuming lepton universality. aMC@NLO: 1954 ⁺⁵⁶₋₈₂ pb.

Per channel, same likelihood model: ee 2094 ⁺¹²²₋₁₁₄, μμ 1931 ⁺³¹₋₃₀, ττ 1981 ⁺¹⁵³₋₁₃₆ pb (z-tautau's own:
1981 ⁺⁷³₋₇₀; ours carries their τh ID p_T-dependence as an extra ±6.3 %). They are compatible:
−2Δln L (one σ against three) = 1.95 for 2 dof, p = 0.38.

This supersedes the combination **frozen on 17 Sep 2026** with τhτh only, 1945 ⁺³¹₋₃₀ pb (git tag
`zmumu-freeze-2026-09-17`, `../FREEZE.md`). The ee and μμ inputs are the frozen ones.

Everything is in `combLieke/`: `README.md` (method and result), `output/plots/`, `output/result.json`,
`docs/`. Rerun with `source ../setup.sh && cd combLieke && python run.py prepare && python run.py condor --submit`
(~15 min), then `python run.py results && python run.py plots`.

## Inputs

| channel | what we took |
|---|---|
| ee | `z-ee/fit.config` with the histograms of the TRExFitter job in `z-ee/Zee_fit.tar.gz` (16 Sep, 14:44 delivery; identical to `datasets/z-ee/zee.root`), σ_ref = 1954.10 pb from the tarball |
| μμ | `z-mumu/fit/zmumu.config` + `fit/fitinputs/zmumu.root` (v2 shape fit, 12 × 5 GeV), σ_ref = 1953.93 pb, measured reconstruction SF, acceptance uncertainty 0.70 % (eight parameters, `z-mumu/zmumu/acceptance.py`) |
| ττ | `z-tautau/fit/ztautau.config` + `fit/fitinputs/ztautau.root` (v4, four final states, 13 regions, τh ID scale factors free), σ_ref = 1944.88 pb, acceptance inside its fit; its published fits `fit/results/*_fit_result.json` for the modelling term and the `tautau_channels` figure |

## Decisions (details and reasons in `combLieke/README.md`)

1. σ = μ_Z × 1953.93 pb in every channel, via a constant per-channel reference factor.
2. The μμ acceptance uncertainties are nuisance parameters in the fit; ee and ττ carry A × ε in their theory
   templates.
3. Correlation by nuisance-parameter name, per `fitting/CONVENTIONS.md`. The ττ muon and electron parameters (POG)
   are decorrelated from μμ's (tag-and-probe) and ee's.
4. **ee template systematics: shape and normalisation are separate parameters.** Only the
   normalisation part is correlated with μμ/ττ. Without this the correlated fit has no minimum,
   because the ee peak shape (GoF p = 4 × 10⁻⁴²) pins its ±5.9 % electron-ID normalisation and
   pulls the shared parameters.
5. `mumu_CRemu` is removed from μμ: same data as the ττ eμ channel. It was a validation region, so μμ does not change.
6. **ττ carries `TauIDpT_tautau`, ±6.3 %**, the p_T dependence of the τh ID scale factor that z-tautau measured with
   its p_T-split fit and asked a conservative combination to carry; the size is read from their two result files.
7. tt̄: the free `mu_ttbar` of ττ (from its eμ control region) and the `XS_TTbar` prior of ee/μμ stay two parameters.
   Sharing the control-region factor with all channels gives +4.5 pb, one constrained `XS_TTbar` everywhere −5.0 pb.
8. MINUIT strategy 3 for every fit (strategy 2 returns status 0 with a covariance forced positive-definite). Two
   exceptions: `split_all_channels` uses strategy 2, and the ranking refits TRExFitter's default escalation
   (`combination/CLAUDE.md` item 9).

## Checked

- **Orthogonality, on data, event by event:** μμ ∩ ee ≤ 97 events, τhτh ∩ μμ = 0, τhτh ∩ ee ≤ 4 events
  (`combLieke/docs/orthogonality.md`). μτh, eτh and eμ veto a second lepton by construction; not re-measured for them.
- **Systematics against ATLAS (arXiv:1603.09222) and CMS (arXiv:2408.03744, 1801.03535):** μμ
  matches. ττ now measures its τh ID in situ (3–5 %). ee lacks trigger SF, charge-misID
  and multijet uncertainties, and its electron-ID uncertainty is inflated by the ECAL gap
  (`combLieke/docs/systematics.md`).
- **Each channel and the combination against the published measurements:** figures
  `channels_vs_published` and `combined_vs_published`.
- **Stability** (figure `variations`, 14 alternative likelihoods): without ee 1937 pb, without ττ 1944 pb; ττ as
  delivered +1.7 pb, with its p_T-split model −13.8 pb; ee electron-ID normalisation at the official ±1.2 % gives
  1939 pb. Every fit converges (MIGRAD, HESSE and MINOS status 0, no forced covariance), HESSE agrees with MINOS on
  μ_Z to 0.07 %, 102/102 ranking refits converge.

## Open, for the next round

1. **z-ee:** veto 1.444 < |η_SC| < 1.566, apply and vary an `HLT_Ele27_WPTight_Gsf` scale factor,
   add charge-misID and multijet uncertainties, check the goodness of fit (coarser bins if needed).
   Then rerun `python run.py all`, and try without `split_shape_norm` for ee
   (`config/channels.json`).
2. **z-tautau:** p_T-dependent τh ID scale factors in the nominal model (then `TauIDpT_tautau` can go), and the
   2.6σ between the eμ and the τh sub-measurements.
3. The event-level orthogonality measurement for μτh, eτh and eμ (`combLieke/docs/orthogonality.md`).
