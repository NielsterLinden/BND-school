# Overview

A measurement of the fiducial `pp -> Z -> mu+ mu-` cross section at 13 TeV,
using CMS 2016 Open Data (SingleMuon, Run2016G+H), entirely from data -- no
simulation anywhere in the chain.

## The measurement in one line

```
                     N_obs - N_bkg
    sigma_fid  =  ---------------------
                     eff_total * L
```

## The chain

| Step | Script | Produces | Doc |
|---|---|---|---|
| 1 | `step1_selection.py` | three regions (OS / SS / e-mu), cutflow, mass spectra | [01](01-selection.md) |
| 2 | `step2_efficiency.py` | tag-and-probe ID / isolation / trigger efficiencies | [04](04-efficiency.md) |
| 3 | `step3_backgrounds.py` | data-driven backgrounds | [05](05-backgrounds.md) |
| 4 | `step4_signal.py` | background-subtracted yield + fit validation | [06](06-cross-section.md) |
| 5 | `step5_crosssection.py` | the cross section and its uncertainties | [06](06-cross-section.md) |
| 6 | `step6_report.py` | `results.json`, `RESULTS.md`, summary plot | -- |

Two topics cut across the steps and have their own pages:

- [02-fsr-photons.md](02-fsr-photons.md) -- final-state radiation
- [03-fake-leptons.md](03-fake-leptons.md) -- non-prompt and fake muons
- [07-systematics.md](07-systematics.md) -- every uncertainty, and how far to trust it

## The three design decisions

**Fiducial, not total.** The acceptance correction that converts a fiducial
cross section to a total one is a generator-level quantity. With no MC
available, quoting a total cross section would mean quoting a number whose
dominant input was guessed. The fiducial volume is defined so that `A = 1`
exactly; anyone with a DY sample can divide by their own acceptance later.

**Data-driven everything.** Efficiencies from tag-and-probe on the Z peak,
backgrounds from same-sign and e-mu control regions. This is forced by the
absence of MC, but it is also the more instructive way to do it.

**Counting, not fitting.** The nominal yield is a background-subtracted count.
The fit is a cross-check on the peak position and resolution, not an alternative
measurement -- see [06](06-cross-section.md) for why conflating the two would be
wrong.

## What limits the result

Not statistics (11 million signal events, ~0.03%). The total uncertainty is
~1.6%, and it is dominated by the **luminosity (1.2%)**, the **external muon
reconstruction efficiency (0.8%)**, and an **assigned trigger-method systematic
(0.5%)**. More data would not help. Better inputs would -- specifically
`TrigObj` branches in the skim, and an orthogonal dataset for an unbiased
trigger efficiency.
