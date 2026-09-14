# Steps 4-5 -- signal extraction and the cross section

`scripts/step4_signal.py`, `scripts/step5_crosssection.py`

## The formula

```
                     N_obs - N_bkg
    sigma_fid  =  ---------------------
                     eff_total * L
```

Four ingredients, each from a different step:

| Symbol | Value from | Page |
|---|---|---|
| `N_obs` | step 1, opposite-sign yield in the window | [01](01-selection.md) |
| `N_bkg` | step 3, same-sign + e-mu | [05](05-backgrounds.md) |
| `eff_total` | step 2, tag-and-probe | [04](04-efficiency.md) |
| `L` | `brilcalc`, 16290.7 pb^-1 | below |

## Why fiducial, and what that means

The **total** `Z` production cross section is what most papers quote, and it is
*not* what this measures. Getting there requires extrapolating from the phase
space the detector can see to all pT and rapidity, through an acceptance factor
`A`:

```
    sigma_total = sigma_fid / A
```

`A` cannot be measured. It is a ratio of two generator-level numbers, and it
comes from a theory calculation or an MC sample. There is no MC here, so
inventing an acceptance would mean quoting a number whose dominant input was
guessed.

Instead the cross section is quoted **in the phase space actually measured**,
where `A = 1` by construction:

| Requirement | Value |
|---|---|
| muons | exactly 2, opposite sign |
| leading muon pT | > 26 GeV |
| subleading muon pT | > 20 GeV |
| muon abs(eta) | < 2.4 |
| muon ID | `Muon_mediumId` |
| muon isolation | `pfRelIso04_all < 0.15` |
| dimuon mass | 60-120 GeV |
| lepton level | **dressed** (muon + recovered FSR photons) |

That last row matters for any comparison with theory -- see
[02-fsr-photons.md](02-fsr-photons.md). Compare against a dressed-lepton
prediction, not a Born-level one.

Anyone with a DY sample can compute `A` for exactly this volume and divide;
nothing in the measurement has to be repeated.

## Signal extraction (step 4)

**Nominal: counting.** `N_sig = N_obs - N_bkg`. No shape assumption, which is
the right choice when the lineshape is distorted by FSR in a way that cannot be
modelled without a generator.

**Cross-check: a fit.** A Voigt profile (Breit-Wigner with the PDG Z width,
convolved with a Gaussian resolution) plus an exponential background, over the
full window.

The fit returns a peak at ~90.8 GeV against the PDG 91.19, and a resolution of
~1.4 GeV. **The chi2/ndf is poor (~6), and the fitted yield is ~5% below the
count. Both are expected.** A single Voigt profile has no radiative tail, so the
exponential background bends to absorb the FSR shoulder and steals events from
the signal integral. The fitted "background" is consequently ~30x the
data-driven estimate and is an upper bound, not a measurement.

This is why the fit is **not** used as the signal-extraction systematic --
|fit - count| compares two different quantities. It is used to validate the
peak position, the resolution, and the absence of a large unaccounted smooth
background. Replacing the Voigt with a Crystal Ball convolved with a
Breit-Wigner would fix the tail and make the fit a genuine alternative; that is
listed as a next step.

The ~0.4 GeV shift of the fitted peak below PDG is the muon momentum scale plus
the residual FSR tail pulling the fit left. A dedicated Rochester correction
would address the first part.

## The efficiency (step 5)

Built per muon, then per event:

```
    eff_total = eff_reco^2
              * <eff_ID>_1 <eff_ISO>_1
              * <eff_ID>_2 <eff_ISO>_2
              * eff_trigger(event)
```

The `<...>` are the 2D tag-and-probe maps **averaged over the observed (pT, eta)
distribution of the signal muons**, taken from the maps step 1 fills. A muon in
a poorly-performing corner of the detector is therefore weighted by how often it
actually occurs, rather than by the flat average over the map.

Leading and subleading muons are averaged separately, because their kinematics
differ -- the subleading isolation efficiency comes out visibly lower, since
softer muons are less isolated.

**Approximation:** correlations between the two muons' efficiencies are
neglected. They are small, because the two muons are well separated in the
detector.

## Luminosity

```
    L = 16290.713420 pb^-1   (Run2016G + Run2016H, runs 278820-284044)
```

From `brilcalc lumi -c web --begin 278820 --end 284044 -i datasets/GRL/GRL.txt`,
computed by the z-ee subgroup and reused here so the three channels share a
normalisation. Uncertainty 1.2% (CMS 2016 legacy calibration).

**This is eras G and H only** -- roughly 16 of the ~36 fb^-1 CMS recorded in
2016. The Open Data release used here contains no other eras. Do not compare the
event counts with a full-2016 analysis without accounting for that.

## Uncertainties

Propagated as independent relative contributions, summed in quadrature.

The largest are, in order: **luminosity (1.2%)**, **muon reconstruction
efficiency (0.8%, external input)**, and the **trigger method bias (0.5%,
assigned not measured)**. All three are systematic, and all three are limits of
the available inputs rather than of the data: 11 million signal events give a
statistical uncertainty of only ~0.03%.

That ordering is the main message of the measurement. More data would not
improve this result. Better inputs would -- specifically `TrigObj` branches in
the skim and an orthogonal dataset for the trigger efficiency.

## Outputs

- `output/plots/step4_fit.png` -- the fit with a pull panel
- `output/plots/step5_systematics.png` -- the uncertainty breakdown
- `output/plots/step6_summary.png` -- the headline plot
- `output/RESULTS.md`, `output/results.json`
