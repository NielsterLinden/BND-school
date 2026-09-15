# 02 — The correlation model

This is the only judgement the combination makes. Everything else is arithmetic on published
numbers, so this page is where a reader should push back.

The model lives in `comb/model.py` as two dictionaries, `CORRELATION` (in-fit systematics, keyed
by the `Category` strings of `fitting/CONVENTIONS.md`) and `ACC_CORRELATION` (acceptance, which is
outside both fits). A category present in a channel's grouped-impact table but absent from the map
raises a `KeyError` — a new systematic must be given a correlation, never default to zero.

## ρ = 1: derived once, used by both channels

| category | μμ | ττ | why |
|---|---:|---:|---|
| `Luminosity` | 23.5 pb | 26.1 pb | the same 16 393.381 pb⁻¹ from the same normtag over the same certified runs, with the same ±1.2 %. There is no sense in which the two channels measured the luminosity independently. |
| `Pileup` | 5.1 pb | 24.0 pb | both build the data pileup profile from the same Open Data luminosity table with σ_mb = 69.2 mb and reweight with the same prescription. Both inherit the same ~4 % bias (`z-mumu/REVIEW.md` F6). |
| `L1 prefiring` | 10.6 pb | 5.7 pb | same maps, same per-event weight. |
| `Background normalisation` | 2.4 pb | 11.0 pb | literally the same nuisance parameters (`XS_TTbar` 6 %, `XS_SingleTop`, `XS_WW`, `XS_WZ`, `XS_ZZ` 10 %) applied to the same MC samples. |
| `Signal modelling` — PDF, α_s, QCD scale, PS ISR/FSR | 5.1 pb | 94.2 pb | the same NNPDF3.1 Hessian members, the same 7-point scale envelope and the same shower weights of the *same* aMC@NLO sample. The effect on each channel's C factor differs; the underlying parameter is one. |

## ρ = 0: different objects, different methods, different events

| category | why |
|---|---|
| `Muon efficiency`, `Muon momentum` vs `Tau ID`, `Tau trigger`, `Tau energy scale` | different objects. μμ derives its muon scale factors in-house from tag-and-probe on the Z peak; ττ takes TauPOG DeepTau scale factors. No common input. |
| `Gammas` | per-bin MC statistics. Both channels use the same DY sample, but the μμ signal region and the τhτh signal region select *disjoint* events (opposite-sign dimuons vs two hadronic taus, with an explicit e/μ veto in ττ). Disjoint events → independent Poisson fluctuations. |
| `Fakes` | μμ uses a same-sign muon fake factor with prompt subtraction; ττ uses a jet→τh fake factor binned in (era, DM, N_jets, p_T) with an OS/SS correction. Different regions, different objects, different statistics. |
| `MET` | ττ only. |
| `Electron efficiency` | μμ only, and zero (it applies to the eμ validation region, which is not in the fit). |
| `Lineshape model` | the ±0.7 % from `z-mumu/REVIEW.md` F3/F4. It is a property of the μμ mass template and of the powheg/aMC@NLO comparison in the μμ phase space; nothing analogous enters the ττ fit, which fits m_ττ from 0 to 350 GeV. |
| `Data statistics` | the two channels read **disjoint primary datasets** (`SingleMuon` and `Tau`) selected by orthogonal triggers. `rho_override` in `model.build` deliberately does not touch this source: its ρ is a fact, not a modelling choice. |

## The one genuine split: `SigModel`

Both channels call their generator systematic `SigModel`, so a TRExFitter MultiFit would
correlate them by name. They are not the same comparison:

| channel | comparison | effect on C |
|---|---|---:|
| μμ | powheg vs aMC@NLO (two NLO generators) | 0.2 % |
| ττ | madgraph **LO** vs aMC@NLO | 7.3 % |

An LO-vs-NLO difference and an NLO-vs-NLO difference are not one nuisance parameter, so the
category is split before the correlation is applied:

```
f      = impact(SigModel) / impact(Signal modelling category)
s_gen  = s_category · f            ρ = 0
s_th   = s_category · sqrt(1 − f²) ρ = 1
```

which preserves the published category total exactly. `impact(SigModel)` comes from each
channel's NP ranking, except for the μμ counting extraction where only the template's
normalisation difference survives (0.17 %, from `sigmodel_powheg_over_nlo` in the fit-result
JSON). Treating `SigModel` as correlated instead is the `sigmodel_correlated` variation: −1.4 pb.

The same decorrelation is needed if the MultiFit is ever run; `fit/comb.config` documents the
`DecorrSysts`/`DecorrSuff` lines that implement it.

## Acceptance

| component | ρ | why |
|---|---:|---|
| PDF | 1 | same NNPDF3.1 set, same Hessian eigenvectors |
| α_s | 1 | same members (0.116 / 0.120) |
| QCD scale | 1 | same 7-point envelope of the same sample |
| PS ISR, FSR | 1 | same shower weights (ττ only quotes them) |
| MC statistics | 0 | disjoint generated events pass the two fiducial volumes |

The QCD-scale choice is the conservative one. The two fiducial volumes respond very differently
to a scale variation — 0.32 % for μμ, 3.4 % for ττ, because the double 40 GeV visible-τ cut sits
on the Z p_T tail — and one could argue they are different observables of the same variation
rather than one shared parameter. Decorrelating it is the `acc_scale_decorrelated` variation:
+1.0 pb, a twentieth of the total uncertainty. The choice does not matter at this precision.

## What the resulting correlation is

Summing the ρ = 1 sources gives **ρ(μμ, ττ) = 0.15**. That number is what produces the negative
ττ weight in the BLUE ([03](03-method.md)); the bracketing variations `rho_zero` and `rho_one`
span 6.5 pb, a fifth of the total uncertainty, so the whole correlation model is a sub-dominant
assumption.
