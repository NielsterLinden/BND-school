# How the combination with the four-channel Z → ττ was built (17 Sep 2026)

**The result and its discussion are in `../README.md` (1949 ⁺³⁰₋₃₀ pb); this file is the record of the two steps
that led there.** The frozen combination (git tag `zmumu-freeze-2026-09-17`) is ee ⊕ μμ ⊕ τhτh.

1. **z-tautau's inputs arrive (b7ec3d9), its fits are not done yet.** The likelihood is rebuilt as described below
   and fitted: **1951.0 ⁺²⁹·⁹₋₂₉.₃ pb**, ττ alone 1981.3 ⁺⁷²·⁷₋₇₀.₀ pb, compatibility p = 0.35 — with ττ exactly as
   delivered. `output/` is left untouched. The ττ standalone number predicted z-tautau's own result, which arrived
   as 1981 ⁺⁷³₋₇₀ pb.
2. **z-tautau's results arrive (e0b0b1d)** with what a combiner has to decide (`z-tautau/docs/11` §7): the τh ID
   scale factor is not flat in p_T (−6.3 % on their μ_Z), two sub-measurements 2.6σ apart, the eμ trigger prior.
   The baseline gains `TauIDpT_tautau`, every item gets a variation, the tt̄ question gets two
   (`../README.md`): **1949.3 ⁺³⁰·²₋₂₉.₆ pb**, p = 0.38. Step 1 is the variation `tautau_as_delivered`.

The numbers in the tables below are those of **step 1** (no `TauIDpT_tautau`); the final ones are in `../README.md`.

## What changed in the likelihood

| | frozen | now | where |
|---|---|---|---|
| ττ regions | `tautau_SR0/1/2` | + `mutau_SR_dm0/1/10/11`, `etau_SR_dm0/1/10/11`, `emu_SR`, `emu_CRtt`: the 13 regions of `z-tautau/fit/ztautau.config`, one `Fit` of the MultiFit. Never the pT-split copies (same events) | `channels.json` → `tautau` |
| ττ signal | `DYtautau` | the 15 `DYtautau_tDM*` templates, read from the ττ metadata (`signal_samples`); `DYtautau_out` is a background | `tautau.signal` |
| ττ acceptance | 3.7 % block of `Acc_*` parameters | none: the ττ theory templates are renormalised to a constant σ(60–120), so `PDF`/`QCDScale`/`PS_*` carry A × ε inside the fit. The μμ `Acc_*` block stays and is now μμ-only | `tautau.acceptance_note` |
| empty bins | `tautau_SR2` bin 1 dropped by the combination | dropped by z-tautau's own config. Re-derived on the v4 inputs: it is the only bin of the 13 regions without any prediction | `tautau.drop_empty_bins_note` |
| `mumu_CRemu` | in the μμ config | **removed** (next section) | `mumu.drop_regions` |
| ττ `Muon*`, `Electron*` | — | decorrelated from μμ / ee (`MuonID_tautau`, …, `ElectronScale_tautau`) | `tautau.rename_systematics` |
| tt̄ | `XS_TTbar` everywhere | ττ: free `mu_ttbar` (from `emu_CRtt`); μμ, ee: `XS_TTbar` 6 % | `tautau.ttbar_note` |
| MINUIT | strategy 2 (escalated to 3 by TRExFitter) | **strategy 3**, `MaximumNumberFCNcalls` 2 × 10⁶ | `fit` |
| where the fits run | one machine, ~10 min | HTCondor DAG, ~15 min wall time including 101 ranking refits | `mf/condor.py` |

## `mumu_CRemu` is dropped, and it costs nothing

`mumu_CRemu` (exactly one muon and one electron, opposite sign) selects the same data as the ττ eμ channel
(`emu_SR`, `emu_CRtt`), so the two cannot be in one likelihood. It is removed from the μμ config by the adapter,
with the `WJets` sample and the `XS_WJets` / `ElectronEff_mumu` parameters that act only there.

* It is `Type: VALIDATION` in `z-mumu/fit/zmumu.config`: it **never was part of the μμ likelihood**. The μμ
  standalone fit with the region as delivered is identical to the one without it (`checks.mumu_cremu_validation`:
  1931.1 ⁺³⁰·⁸₋₃₀.₁ pb, Δ = 0.0), and `without_tautau` reproduces the frozen 1943.9 pb exactly.
* What μμ would gain by fitting it as a control region (`checks.mumu_cremu_control`): −0.4 pb (−0.02 %) and a
  2.6 % smaller uncertainty (30.0/29.3 instead of 30.8/30.1 pb). That is the whole price of leaving the eμ data to ττ.

## The two decisions z-tautau left to the combination

1. **Are the ττ muon and electron parameters the same quantity as μμ's and ee's? — No, decorrelated.** ττ takes
   them from the POG jsons, z-mumu measures its own with tag-and-probe, and z-ee's `ElectronID` is the ±5.9 % ECAL-gap
   artefact. Sharing a name would claim one measurement. Size of the choice (`tautau_leptons_correlated`, everything
   shared by name, `ElectronReco`/`ElectronScale` mapped onto ee's `ElectronRECO`/`ElectronES`): **−1.0 pb**, same
   uncertainty.
2. **`mu_ttbar` or `XS_TTbar`? — Both, each in its own channel.** ττ has no `XS_TTbar`; its tt̄ is scaled by the free
   `mu_ttbar` = 1.105 ± 0.026, which absorbs the eμ trigger efficiency and is therefore not the cross section that
   `XS_TTbar` describes. No yield is described twice. The alternative (`tautau_ttbar_constrained`, now `ttbar_xs_constrained`: one 6 % `XS_TTbar`
   for all channels, no free factor) gives **−6.7 pb** and a smaller uncertainty (⁺²⁹·¹₋₂₈.₅): with tt̄ tied to its
   theory cross section the 46 k events of `emu_CRtt` become a luminosity × efficiency monitor, and a 10 % tt̄ excess
   is pushed into `Lumi`, `EmuTrigger` and `BTag`. That is a different measurement; it is not the baseline.
   The same mechanism is why `mu_ttbar` is second in the refit-based ranking (±19 pb): fixing it at ±1σ forces
   `Lumi` to compensate in the control region. It is a correlation, not an uncertainty that is missing from the total.

The shared theory parameters are constrained by the four ττ final states through their acceptance ratios, as
z-tautau predicted (`QCDScale` 0.69 → 0.53, `PS_FSR` 0.56 → 0.36, `Pileup` 0.66 → 0.43). Because the μμ acceptance
sits in its own `Acc_*` parameters, this cannot shrink the μμ acceptance uncertainty; with the ττ theory parameters
decorrelated altogether (`tautau_theory_decorrelated`) the result moves by −0.4 pb.

| likelihood | σ [pb] | Δ |
|---|---|---|
| baseline | 1951.0 ⁺²⁹·⁹₋₂₉.₃ | |
| without ee | 1938.9 ⁺²⁹·⁹₋₂₉.₃ | −12.1 |
| without ττ | 1943.9 ⁺³⁰·⁸₋₃₀.₂ | −7.1 |
| shape/normalisation split in all channels | 1959.3 ⁺³¹·⁴₋₃₀.₈ | +8.3 |
| ee: split only the shared parameters | 1885.1 ⁺²⁶·⁴₋₂₆.₀ | −65.9 |
| ee: `ElectronID` normalisation ±1.2 % (diagnostic) | 1940.4 ⁺²⁸·⁶₋₂₈.₁ | −10.6 |
| ττ lepton parameters shared with μμ / ee | 1950.1 ⁺²⁹·⁹₋₂₉.₄ | −1.0 |
| ττ theory parameters not shared | 1950.7 ⁺³⁰·⁰₋₂₉.₄ | −0.4 |
| ττ: `XS_TTbar` instead of free `mu_ttbar` | 1944.3 ⁺²⁹·¹₋₂₈.₅ | −6.7 |

## Fit health, and why MINUIT strategy 3

Every fit (combined, three-POI, three standalone, stat-only, nine variations, two checks) ends with MIGRAD 0,
HESSE 0, MINOS 0 and **no covariance forced positive-definite**; HESSE/MINOS on μ_Z = 1.0014; no constrained
parameter has a post-fit error above 1.02; 101/101 ranking refits converged.

That needed strategy 3. With strategy 2 every fit also ended with status 0 and the same MINOS interval, but HESSE
ran on a matrix "forced pos-def by adding to diagonal": HESSE error on μ_Z 0.0018 instead of 0.0152, `Lumi`
post-fit error 0.15, μ_Z most correlated with a single ττ MC-statistics γ; the three-POI fit lost the MINOS upper
error of μ_Z(ee) and `without_ee` ended with status 1. The frozen fit asked for strategy 2 as well, but there MIGRAD
failed at 2 and TRExFitter escalated to 3 by itself; the four-channel ττ likelihood "converges" at 2, so nothing
escalates. Strategy 3 needs `MaximumNumberFCNcalls` raised (the ee + μμ fit otherwise exhausts the default in HESSE,
status 300, segfault). **For z-tautau:** the ττ standalone fit shows the same forced matrix at strategy 2 and is clean
with `FitStrategy: 3` + `MaximumNumberFCNcalls: 2000000` — likely what keeps their own fits from converging.

Also fixed on the way: TRExFitter writes the correlation matrix of `Fits/*.txt` with its rows in reverse parameter
order; `mf/results._corr` read other elements. The POI correlations of the three-POI fit are ee–μμ 0.39, ee–ττ 0.43,
μμ–ττ 0.67 (shared luminosity). The `compatibility.correlations` of the frozen `output/result.json` are affected
(nothing else used them).

## Still open after step 1 (see `../README.md` for the state now)

* **Final `results` and `plots`**, waiting for z-tautau's published fit result. Expect their own number near
  1981 pb (μ_Z ≈ 1.019 against their 1944.9 pb reference) if their inputs stay as they are.
* `docs/orthogonality.md` was measured for ee, μμ and τhτh. μτh, eτh and eμ veto a second muon or electron by
  construction, and the one known overlap (`mumu_CRemu`) is removed, but the event-level check on data has not been
  repeated for the three new final states (the ττ v4 data ntuples with run/lumi/event are in
  `$BND_TAUTAU_CACHE/ntuples_v4/`).
* `docs/systematics.md` and `checks/systematics.json` still describe the τhτh-only ττ.
* The ee limitations of `../README.md` ("The ee channel") are unchanged.
