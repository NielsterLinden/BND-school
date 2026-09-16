# Handoff – Combination

## Team

- Lieke Gijsen
- Caspar van Eck

## What we did

We combined all three channels — `z-mumu` (Z → μ⁺μ⁻), `z-tautau` (Z → τhτh) and `z-ee` (Z → e⁺e⁻)
— into one Z cross section.

> **σ(pp → Z/γ* → ℓℓ, 60 < m_ℓℓ < 120 GeV) = 1980 ± 32 pb**
> = 1980 ± 0.5 (stat) ± 22 (syst) ± 7 (acc) ± 21 (lumi) pb, per lepton flavour, assuming lepton
> universality.
> **χ²/ndf = 9.67/2 (p = 0.008)** — the three channels are *not* mutually consistent.
>
> **R = σ(Z→ττ)/σ(Z→μμ) = 1.08 ± 0.13** (+0.6σ) and
> **R = σ(Z→ee)/σ(Z→μμ) = 1.064 ± 0.021** (+3.0σ).

Updated on 16 Sep 2026: z-ee added, z-tautau moved to v3 (DeepTau Tight nominal). The previous
round (μμ ⊕ ττ) gave 1931 ± 35 pb.

Everything: [`result.md`](result.md) (generated), [`docs/`](docs/), [`README.md`](README.md),
[`CLAUDE.md`](CLAUDE.md). Deck: [`presentation/`](presentation/). What to ask z-ee for:
[`../ASK_Z_EE.md`](../ASK_Z_EE.md).

**Checked on 16 Sep 2026** (`docs/05-vs-published.md`): the three selections are orthogonal to
better than 0.002 % (μμ ∩ ee = 61 ZZ→4ℓ events; ττ vetoes e and μ); each channel is compared with
the published measurement of the same decay (ATLAS Z→ee and Z→μμ, CMS Z→ττ); and the five
systematics CMS/ATLAS carry that we do not are listed — four of them in z-ee.

## Read this before quoting the number

**The z-ee input is biased by +12 % and we know exactly why.** Its `PDF` and `QCDScale` templates
are raw LHE weight envelopes that were never renormalised to a constant yield
(`fitting/CONVENTIONS.md` §3 requires it), so `QCDScale` enters the fit as a ±5.87 %
*normalisation* of the signal — degenerate with the POI. The fit pulls it to −1.88σ, rescaling the
prediction by κ = 0.893, and `mu_signal` = 1.0509 is measured against that rescaled prediction.
The data-over-prediction ratio is 0.9700; μ̂ × (response of every NP) = 0.9706 closes with it to
0.05 %, so the diagnosis is arithmetic, not opinion.

* Baseline (the channel's own number): σ(ee) = 2054 ± 40 pb → combination 1980 ± 32 pb.
* `ee_normfix` (κ divided out): σ(ee) = 1833 pb → combination 1882 ± 30 pb.
* Without ee at all: 1931 ± 35 pb.

We kept the channel's own number as the baseline — this folder follows the channels, it does not
re-fit them — and flagged it everywhere. **The honest reading today is 1931 ± 35 pb from μμ ⊕ ττ,
with an ee channel that does not yet agree with it.**

### What z-ee needs to change — full list in [`../ASK_Z_EE.md`](../ASK_Z_EE.md)

1. **Renormalise the `PDF`/`QCDScale` envelopes to a constant yield** — divide each LHE replica by
   its own Σw, or rescale each replica histogram to the nominal integral. Only the shape may vary.
   *This is the whole 99 pb.*
2. Accumulate the envelopes per LHE flavour, not over the inclusive DY sample.
3. Conventions §1–§4: POI `mu_Z`, `DY_ee` as `Type: SIGNAL`, region `ee_SR`, the §3 sample and NP
   names, `z-ee/fit/fitinputs/zee.root`, and a published `z-ee/fit/results/zee_fit_result.json`
   from `fitting/run_trex.py` — then `tools/extract_zee.py` can be deleted.
4. Run the saturated-model **goodness of fit**. z-ee is the only channel without one.
5. Fit in coarser bins than 60 × 1 GeV. `Pileup` −2.9σ, `L1Prefiring` +3.6σ and `ElectronID`
   constrained to 0.08 are a data/MC shape mismatch being absorbed by calibration parameters;
   z-mumu hit the same thing (`REVIEW.md` F3/F4) and fixed it with 12 × 5 GeV.
6. `Wjets` has 24 % MC statistics and ~30 bins that were negative before being floored to +10⁻⁶.
   Flooring keeps TRExFitter happy; smooth the template or merge W+jets into the other
   backgrounds.
7. `HistoPath` in `fit.config` points at `/project/atlas/users/lvdurenw/...` — make it
   repo-relative so anyone can rerun the fit.

Already fixed by z-ee on 16 Sep 2026: the reference cross section (their 1954.1032 pb agrees with
our independent 1954.1110 pb to 4 × 10⁻⁶) and the negative bins.

## Method, in one paragraph

A covariance (BLUE) combination of the three channels' published TRExFitter v1.8.0 results, with
the shared systematics correlated through the `Category` strings of `../fitting/CONVENTIONS.md`
(ρ = 1 for luminosity, pileup, L1 prefiring, background cross sections, PDF/α_s/scale/PS; ρ = 0
for the object calibrations, MC statistics, fakes, the generator comparison and the z-ee fit
residual). Cross-checked against an explicit profile likelihood that carries the ττ channel's
asymmetric uncertainty — same answer to 0.3 pb. Combination is done on the cross sections, **not**
on `mu_Z`, because the three channels normalise to references that differ by up to 0.47 %
(`docs/01-inputs.md`).

## Inputs we used

| channel | what we took | value | weight |
|---|---|---|---|
| μμ | the **v2 shape fit** (12 × 5 GeV, two-sided `SigModel`, MINOS), the channel's baseline after its review fixes | μ_Z = 0.9881 ± 0.0159, σ = 1931 ± 35 pb | +0.596 |
| ττ | **v3**: DeepTau **Tight** on both legs, MC-subtracted fake factor (its `nominal_variant`) | μ_Z = 1.071 ⁺⁰·¹¹⁴₋₀.₁₀₀, σ = 2082 ± 253 pb | −0.0003 |
| ee | the published TRExFitter job `z-ee/Zee_fit.tar.gz` (16 Sep 2026), converted by `tools/extract_zee.py` | μ_signal = 1.0509 ± 0.0202, σ = 2054 ± 40 pb | +0.404 |

z-tautau v3 publishes only its nominal, so the ττ cross-check rows the previous round carried
(`nosub`, Tight-WP) are gone — Tight *is* the nominal now.

`run_combination.py --check` asserts all three against the channels' published files and fails if
anything changes. **If a channel re-publishes, run that first** — it will tell you exactly what
moved.

## The main messages

1. **Adding ee finally makes combining worth something**: 35.4 → 31.6 pb, a 10.7 % improvement,
   where μμ ⊕ ττ gained 0.00 %. It is limited to 10.7 % rather than 29 % because ρ(μμ, ee) = 0.44,
   almost all of it luminosity.
2. **But the χ² is 9.67/2.** See above. Fix z-ee before quoting 1980 ± 32 pb anywhere.
3. **The floor is still the luminosity**: 21 of the 32 pb, fully correlated, irreducible.
4. Data statistics are 0.5 pb — 0.02 %. The MC statistics (`Gammas`) are 10.6 pb.

## Not done

- **No TRExFitter MultiFit.** The channel workspaces are not all committed and were not rebuilt,
  so the shared NPs are correlated by hand instead of profiled jointly. The config is generated and
  syntax-checked ([`fit/comb.config`](fit/comb.config)) and
  [`fit/run_multifit.py`](fit/run_multifit.py) has the full recipe. This matters more than it used
  to: ee now carries 40 % of the weight, and a joint fit would expose its `QCDScale` problem
  structurally instead of as a variation.
- **No χ² inflation.** BLUE returns ±32 pb assuming the inputs are consistent. No PDG-style scale
  factor (√(χ²/ndf) = 2.2) has been applied — we prefer to fix the input.
- The caveats we inherited and did not fix are listed at the end of `result.md`.

## Next steps, in order of what they would change

1. **z-ee re-runs with renormalised theory templates** (list above). 99 pb, and it is the
   difference between a quotable number and one with a p = 0.008 caveat attached.
2. Run the MultiFit and compare.
3. Settle the μμ 60–80 GeV lineshape (`z-mumu` open issue 3): 11.1 pb between counting and fitting.
4. More simulated events: `Gammas` is 10.6 pb on the combination and 23.8 pb inside the ee channel.
5. Add eτh/μτh to ττ so the τh ID is constrained in situ — that improves the universality test,
   not the cross section.
