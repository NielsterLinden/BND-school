# Answer to `REVIEW_v4.md`

What was changed in the four-channel measurement in response to the review of 17 September 2026, where it
lives, and what was deliberately left undone. The review itself is `REVIEW_v4.md`; the numbers after the
changes are in `output/RESULTS.md`.

The repository was also reorganised at the same time: **the four-channel measurement is now the only
measurement**. It owns the canonical paths (`fit/`, `output/`, job `ztautau`, per-channel jobs
`ztautau_<ch>`) that `fitting/CONVENTIONS.md` and the combination expect, the τhτh-only v3 fit and its
outputs are gone, and the part of the v3 chain that is still needed — the τhτh Data and fake templates —
is built by `run_tautau_base.py` into `fit/fitinputs/tautau_base.root`. What the combination gets, and what
it must know before using it, is `docs/11-combination-inputs.md`.

## 1. Summary

| # | severity | finding | what was done |
|---|---|---|---|
| 1 | critical | the quoted expected uncertainty was the HESSE error of a failed Asimov MINOS | the Asimov fit now runs from its own config with `FitStrategy: 2`, the MINOS status is parsed, and **nothing is quoted if it fails**; the observed fit runs last so the correlation matrix and the NP plots are the observed ones | 
| 2 | critical | the four-channel number is a compromise between two sub-measurements that disagree | both are now fitted, quoted side by side with their difference, and the report says in so many words that the combined number is not independent of the eμ trigger prior |
| 3 | major | the eμ trigger prior was 5.4 %: the paper's 2 % applied once per leg | **fixed in the weights**: the flat 2 % is applied once per event (2.7–2.9 % total with the in-situ statistics). The old treatment is kept as the cross-check fit `ztautau_emutrig2x` |
| 4 | major | the τhτh / (ℓτh)² lever assumes one τh ID scale factor from 30 GeV upwards | the ℓτh regions are now also filled split at pT(τh) = 40 GeV, and the cross-check fit `ztautau_ptsplit` gives the low-pT part its own scale factors and measures the difference |
| 5 | major | the grouped impacts over-shoot the total by 1.6 | the over-shoot is computed, printed, and written to `results.json` together with the rescaled groups; `EmuTrigger` was moved out of `Electron efficiency` into its own category so the table is readable; the report states that the total is the MINOS one |
| 6 | minor | `FakeFrac_etau` pulled +2.6σ | the prior on the W+jets fraction of the application region was 20 %, an estimate; it is now 40 % (`config.LTAU_FF_FRAC_SYST`), which is what the data say it should have been |
| 7 | minor | eτh electrons sit on the Ele27 turn-on, covered by one flat 3 % | the turn-on bin (29–35 GeV) has its own nuisance parameter `ElectronTrigger_lowpt`, sized by the step to the next pT bin of the in-situ table (4 % / 2 % / 13 % in the three η bins used) |
| 8 | minor | `MET_Unclustered` acts as a normalisation knob | **not split**; explained instead (section 4) |
| 9 | minor | report hygiene | `mu_ttbar` is a labelled triplet, the statistical error is rounded, the "expected" row is honest, and the plots copied to `output/plots/` come from the observed fit |

## 2. Findings 1 and 5: what the fit reports

`scripts/step5_fit.py` now runs, in this order, `h`, `w`, the stat-only fit, the Asimov fit, and only then
`f`, `d`, `p`, `i`, `r`. TRExFitter writes `CorrelationMatrix.yaml`, `CorrMatrix.png` and the nuisance-parameter
plots without a suffix, so whichever fit ran last owned them; the observed fit now does.

The Asimov fit is run from `fit/ztautau_asimov.config`, a copy of the configuration with `FitStrategy: 2`.
`minos_status` is parsed out of the log of both fits and stored. If the Asimov MINOS does not converge, the
HESSE numbers are stored under `asimov_hesse_not_quoted` and the report prints *"not quoted: MINOS did not
converge on the Asimov data set"* instead of a number.

The grouped impacts share nuisance parameters between categories — `mu_ttbar`, the eμ trigger efficiency and
`mu_Z` are one degeneracy, so `NormFactors` and `Emu trigger` both contain it. Their quadrature sum is
therefore larger than the MINOS total. `step5_fit.py` records `grouped_impact_quadrature_sum` and
`grouped_impact_scale` (the factor that makes the sum equal the total), `RESULTS.md` prints both and says
that the uncertainty of the measurement is the MINOS one, and `results.json` offers the combination
`groups_rescaled` next to the raw `groups` (`docs/11-combination-inputs.md` §6).

## 3. Findings 3, 6 and 7: the priors that were estimates

* **eμ cross trigger.** `analysis_v4._insitu_sf` took `(SF ± stat) × (1 ± 2 %)` per leg with the same sign on
  both legs, so the 2 % of CMS arXiv:1801.03535 (a per-*channel* number) entered twice and linearly. The flat
  term is now applied to one leg only (`flat=0.0` for the muon leg); the per-leg statistical errors of the
  in-situ tables are unchanged. The eμ signal normalisation moves by 2.7–2.9 % per σ instead of 5.4 %.
  `ztautau_emutrig2x` repeats the fit with the variation doubled, i.e. with the old treatment, and
  `RESULTS.md` quotes it: that difference is the size of the choice.
* **eτh Ele27 turn-on.** `analysis_v4._ele27_sf` splits the in-situ electron-trigger uncertainty in two:
  `ElectronTrigger` above 35 GeV (bin statistics + the flat 2 %) and `ElectronTrigger_lowpt` below
  (bin statistics + `ele27_turnon_rel`, the step to the next pT bin of the measured table). The two are
  uncorrelated, so a mismodelled turn-on can no longer be traded against the plateau.
* **ℓτh fake composition.** The W+jets fraction of the application region comes from simulation and carried a
  ±20 % prior; the fit pulled it +2.6σ (eτh) and +1.5σ (μτh) and constrained it to ≈ 0.6 of the prior, i.e.
  the data said the prior was about twice too tight. It is now ±40 % (`config.LTAU_FF_FRAC_SYST`), which the
  fit still constrains to ≈ 25 %.

None of these three is a change of a central value: they change how far a parameter is allowed to move.

## 4. Finding 8: why `MET_Unclustered` was left alone

The unclustered-energy variation changes the ℓτh signal *yield* by −3.8 % / +2.9 % because it moves events
across two cuts that both use the MET: m_T(ℓ, MET) < 40 GeV in the region definition and the MET term of the
di-τ likelihood mass, which decides which m_tt bin an event lands in (and whether it stays inside 0–350 GeV).
That is a genuine acceptance effect of the same source, not an artefact of the template construction, so
splitting it into an uncorrelated "shape" and "normalisation" pair would give the fit a second free parameter
where the physics has one. It stays one nuisance parameter; its size, pull and post-fit constraint are in
`RESULTS.md`, and a combination sees it in the `MET` category (ρ = 0, our regions only).

The same reasoning does *not* apply to the other two constrained-below-0.5 parameters: `QCDScale` and
`Pileup` are constrained because the τhτh/ℓτh ratio is sensitive to them, which is the intended in-situ
behaviour of this fit.

## 5. What a reader of the result should carry away

The review's central point stands and is now in the report rather than in a review document: with the τh ID
scale factors free — which is the point of the four-channel fit — the eμ channel and the three τ channels
measure `mu_Z` in different ways and do not agree well. The combined fit is the correct likelihood answer,
and it is the number this channel quotes; it is also the number that a combination should treat as
correlated with the eμ trigger efficiency and with `mu_ttbar`. `docs/11-combination-inputs.md` §7 lists what
has to be known before the number is used, and `output/RESULTS.md` carries both sub-results and the
cross-check fits next to the headline.
