# Handoff – Combination

## Team

- Lieke Gijsen
- Caspar van Eck

## What we did

We combined the three channels, `z-ee` (Z → e⁺e⁻), `z-mumu` (Z → μ⁺μ⁻) and `z-tautau` (Z → τhτh), in
**one TRExFitter v1.8.0 MultiFit**: one signal strength for all channels, and every shared nuisance
parameter fitted jointly.

> **σ(pp → Z/γ\* → ℓℓ, 60 < m_ℓℓ < 120 GeV) = 1948 ⁺³³₋₃₂ pb**
> = 1948 ± 0.5 (stat) ⁺³³₋₃₂ (syst, of which luminosity 20.7 pb and acceptance 10.6 pb) pb,
> per lepton flavour, assuming lepton universality. aMC@NLO: 1954 ⁺⁵⁶₋₈₂ pb.

Per channel, same likelihood model: ee 2094 ⁺¹²²₋₁₁₄, μμ 1931 ⁺³³₋₃₃, τhτh 2082 ⁺²³⁷₋₂₀₇ pb. They are
compatible: −2Δln L (one σ against three) = 2.41 for 2 dof, p = 0.30.

Everything is in `combLieke/`: `README.md` (method and result), `output/plots/`, `output/result.json`,
`docs/`. Rerun with `source ../setup.sh && cd combLieke && python run.py all` (~10 min).

## Inputs

| channel | what we took |
|---|---|
| ee | `z-ee/fit.config` with the histograms of the TRExFitter job in `z-ee/Zee_fit.tar.gz` (16 Sep, 14:44 delivery; identical to `datasets/z-ee/zee.root`), σ_ref = 1954.10 pb from the tarball |
| μμ | `z-mumu/fit/zmumu.config` + `fit/fitinputs/zmumu.root` (v2 shape fit, 12 × 5 GeV), σ_ref = 1953.93 pb, acceptance uncertainty 0.61 % |
| ττ | `z-tautau/fit/ztautau.config` + `fit/fitinputs/ztautau.root` (v3, DeepTau Tight, three BDT categories), σ_ref = 1944.88 pb, acceptance uncertainty 3.7 % |

## Decisions (details and reasons in `combLieke/README.md`)

1. σ = μ_Z × 1953.93 pb in every channel, via a constant per-channel reference factor.
2. μμ and ττ acceptance uncertainties are nuisance parameters in the fit (PDF, α_s, scale and PS
   correlated between μμ and ττ; MC statistics per channel).
3. Correlation by nuisance-parameter name, per `fitting/CONVENTIONS.md`.
4. **ee template systematics: shape and normalisation are separate parameters.** Only the
   normalisation part is correlated with μμ/ττ. Without this the correlated fit has no minimum,
   because the ee peak shape (GoF p = 7 × 10⁻⁴³) pins its ±5.9 % electron-ID normalisation and
   pulls the shared parameters.
5. The empty 0–40 GeV bin of the highest ττ BDT category (0 data, 0 prediction) is dropped.

## Checked

- **Orthogonality, on data, event by event:** μμ ∩ ee ≤ 97 events, ττ ∩ μμ = 0, ττ ∩ ee ≤ 4 events
  (`combLieke/docs/orthogonality.md`).
- **Systematics against ATLAS (arXiv:1603.09222) and CMS (arXiv:2408.03744, 1801.03535):** μμ
  matches. ττ is complete but its τh ID is not constrained in situ. ee lacks trigger SF, charge-misID
  and multijet uncertainties, and its electron-ID uncertainty is inflated by the ECAL gap
  (`combLieke/docs/systematics.md`).
- **Each channel and the combination against the published measurements:** figures
  `channels_vs_published` and `combined_vs_published`.
- **Stability:** without ee 1932 pb; ee electron-ID normalisation at the official ±1.2 % gives 1934 pb,
  with ee alone at 1918 ± 39 pb. Every fit converges (MIGRAD, HESSE and MINOS status 0).

## Open, for the next round

1. **z-ee:** veto 1.444 < |η_SC| < 1.566, apply and vary an `HLT_Ele27_WPTight_Gsf` scale factor,
   add charge-misID and multijet uncertainties, check the goodness of fit (coarser bins if needed).
   Then rerun `python run.py all`, and try without `split_shape_norm` for ee
   (`config/channels.json`).
2. **z-tautau:** add eτh/μτh to constrain the τh ID in situ, as CMS does.
3. A and C theory variations are separate parameters for μμ and ττ. Only matters for ττ, which carries a
   few per cent of the weight.
