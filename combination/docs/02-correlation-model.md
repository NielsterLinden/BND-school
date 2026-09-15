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
| `Luminosity` | 22.7 pb | 19.9 pb | the same 16 393.381 pb⁻¹ from the same normtag over the same certified runs, with the same ±1.2 %. There is no sense in which the two channels measured the luminosity independently. |
| `Pileup` | 2.6 pb | 19.8 pb | both build the data pileup profile from the same Open Data luminosity table with σ_mb = 69.2 mb and reweight with the same prescription. Both inherit the same open issue (z-mumu's profile is fitted to N_PV rather than taken from the official `puWeights` file). |
| `L1 prefiring` | 9.9 pb | 3.4 pb | same maps, same per-event weight. |
| `Background normalisation` | 1.6 pb | 34.0 pb | literally the same nuisance parameters (`XS_TTbar` 6 %, `XS_SingleTop`, `XS_WW`, `XS_WZ`, `XS_ZZ` 10 %) applied to the same MC samples. |
| `Signal modelling` — PDF, α_s, QCD scale, PS ISR/FSR | 5.7 pb | 34.1 pb | the same NNPDF3.1 Hessian members, the same 7-point scale envelope and the same shower weights of the *same* aMC@NLO sample. The effect on each channel's C factor differs; the underlying parameter is one. |

## ρ = 0: different objects, different methods, different events

| category | why |
|---|---|
| `Muon efficiency`, `Muon momentum` vs `Tau ID`, `Tau trigger`, `Tau energy scale` | different objects. μμ derives its muon scale factors in-house from tag-and-probe on the Z peak; ττ takes TauPOG DeepTau scale factors. No common input. |
| `Gammas` | per-bin MC statistics. Both channels use the same DY sample, but the μμ signal region and the τhτh signal region select *disjoint* events (opposite-sign dimuons vs two hadronic taus, with an explicit e/μ veto in ττ). Disjoint events → independent Poisson fluctuations. |
| `Fakes` | μμ uses a same-sign muon fake factor with prompt subtraction; ττ uses a jet→τh fake factor binned in (era, DM, N_jets, p_T) with an OS/SS correction. Different regions, different objects, different statistics. |
| `MET` | ττ only. |
| `Electron efficiency` | μμ only, and zero (it applies to the eμ validation region, which is not in the fit). |
| `Data statistics` | the two channels read **disjoint primary datasets** (`SingleMuon` and `Tau`) selected by orthogonal triggers. `rho_override` in `model.build` deliberately does not touch this source: its ρ is a fact, not a modelling choice. |

## The split of `Signal modelling`, and why it is now one-sided

Until z-tautau v2.0 both channels fitted a generator systematic called `SigModel`, so a
TRExFitter MultiFit would have correlated them by name even though they are not the same
comparison:

| channel | comparison | status |
|---|---|---|
| μμ | powheg vs aMC@NLO (two NLO generators), two-sided, built inside a common 50 < m_LHE < 120 GeV window | **fitted**; 0.40 % of μ_Z |
| ττ | madgraph **LO** vs aMC@NLO, C_LO/C_NLO = 0.867 | **reported, not fitted** since v2.1 |

z-tautau removed its version deliberately (`z-tautau/docs/07`, its `REVIEW.md` §3.4): the LO
sample is simply the worse model of the visible-τ p_T spectrum whose slope at the 40 GeV threshold
is what the number measures, an uncertainty built as "nominal minus a worse model" is a
placeholder, it double-counts `QCDScale`/`PS_ISR`, and under its old name it would have been
correlated in a MultiFit with the μμ NLO-vs-NLO comparison. Its `Signal modelling` category is
now pure PDF/α_s/scale/PS.

The split survives because the μμ side still needs it — and because a channel that reinstates its
own generator parameter would be handled without a code change:

```
f      = impact(SigModel) / impact(Signal modelling category)
s_gen  = s_category · f            ρ = 0
s_th   = s_category · sqrt(1 − f²) ρ = 1
```

which preserves the published category total exactly. `impact(SigModel)` comes from each channel's
NP ranking — zero for ττ, 0.40 % for μμ — except for the μμ counting extraction, where only the
template's normalisation difference survives (0.28 %, `meta.sigmodel.C_ratio_powheg_over_nlo` in
the fit-result JSON). The `sigmodel_correlated` variation is therefore **exactly null now**
(0.0 pb): there is nothing on the ττ side to correlate it with. It is kept in the table to make
that visible rather than silent.

The same decorrelation would still be needed if the MultiFit were run *and* z-tautau switched its
`SigModel_tautau` on (`config.SIGMODEL_IN_FIT = True`); `fit/comb.config` documents the
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
+1.6 pb, a twentieth of the total uncertainty. The choice does not matter at this precision.

## What the resulting correlation is

Summing the ρ = 1 sources gives **ρ(μμ, ττ) = 0.114**. That number is what produces the negative
ττ weight in the BLUE ([03](03-method.md)) — barely, since it only just exceeds the ratio of the
two uncertainties, 35.4/322.0 = 0.110. The bracketing variations `rho_zero` and `rho_one` span
8.7 pb, a quarter of the total uncertainty, so the whole correlation model is still a sub-dominant
assumption: `mumu counting`, a single choice inside one channel, moves the result by more
(+11.2 pb).
