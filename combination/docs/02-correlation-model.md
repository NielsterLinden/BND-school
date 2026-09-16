# 02 — The correlation model

This is the only judgement the combination makes. Everything else is arithmetic on published
numbers, so this page is where a reader should push back.

The model lives in `comb/model.py` as two dictionaries, `CORRELATION` (in-fit systematics, keyed
by the `Category` strings of `fitting/CONVENTIONS.md`) and `ACC_CORRELATION` (acceptance, which is
outside the μμ and ττ fits), plus `STRUCTURAL_RHO` for the two sources that are not categories at
all. A category present in a channel's grouped-impact table but absent from the map raises a
`KeyError` — a new systematic must be given a correlation, never default to zero.

Sizes below are in pb on σ(60 < m < 120 GeV), from `output/combination_result.json`.

## ρ = 1: derived once, used by every channel

| category | μμ | ττ | ee | why |
|---|---:|---:|---:|---|
| `Luminosity` | 22.7 | 13.7 | 18.8 | the same 16 393.381 pb⁻¹ from the same normtag over the same certified runs, with the same ±1.2 %. There is no sense in which the channels measured the luminosity independently. |
| `Pileup` | 2.6 | 10.7 | 2.6 | all three build the data pileup profile from the 2016 luminosity information with σ_mb = 69.2 mb and reweight with the same prescription. (ee uses the official `Collisions16_UltraLegacy_goldenJSON` `puWeights`; μμ fits the profile to N_PV because the official file is unreachable for its sample — the *input* is the same measurement either way, and decorrelating it would change the combination by well under a pb.) |
| `L1 prefiring` | 9.9 | 3.6 | 12.4 | the same per-event `L1PreFiringWeight` branch of the same NanoAOD. |
| `Background normalisation` | 1.6 | 49.9 | 7.0 | literally the same nuisance parameters (`XS_TTbar` 6 %, `XS_WW`/`XS_WZ`/`XS_ZZ`/`XS_WJets` 10 %, `XS_DYtautau` 5 %) applied to the same MC samples. |
| `Signal modelling` — PDF, α_s, QCD scale, PS ISR/FSR | 5.7 | 24.4 | 7.7 | the same NNPDF3.1 Hessian members, the same 7-point scale envelope and the same shower weights of the *same* aMC@NLO sample. The effect on each channel's C factor differs; the underlying parameter is one. **Caveat:** on the ee side this category is not a C-factor variation but the un-renormalised envelope ([01](01-inputs.md)), so ρ = 1 is the conservative reading of a quantity that should not have been in the fit in that form at all. |

## ρ = 0: different objects, different methods, different events

| category | why |
|---|---|
| `Muon efficiency`, `Muon momentum` (μμ) vs `Electron efficiency` (ee) vs `Tau ID`, `Tau trigger`, `Tau energy scale` (ττ) | different objects, and each channel is the only user of its own. μμ derives its muon scale factors in-house from tag-and-probe on the Z peak; ee takes EGM-POG reco and Medium-ID scale factors from `jsonpog-integration`; ττ takes TauPOG DeepTau scale factors. No common input. |
| `Gammas` | per-bin MC statistics. All three use the same DY sample, but the three signal regions select *disjoint* events (two muons / two hadronic taus with an e-μ veto / two electrons). Disjoint events → independent Poisson fluctuations. |
| `Fakes` | μμ uses a same-sign muon fake factor with prompt subtraction; ττ a jet→τh fake factor binned in (era, DM, N_jets, p_T) with an OS/SS correction. ee has no data-driven fake estimate. Different regions, different objects, different statistics. |
| `MET` | ττ only. |
| `Data statistics` | the three channels read **disjoint primary datasets** (`SingleMuon`, `Tau`, `Electron`). ττ vetoes electrons and muons explicitly; the only process that can enter both the μμ and the ee selection is ZZ → 4ℓ, 61 of 10.38 M μμ SR events (0.0006 %, checked by z-mumu, `fitting/CONVENTIONS.md` §6). `rho_override` in `model.build` deliberately does not touch this source: its ρ is a fact, not a modelling choice. |
| `Fit residual` | ee only. The part of the z-ee published MINOS error its own grouped impacts do not add up to (0.92 % of μ_Z, 18.0 pb) — a property of *that fit's* correlation matrix, so it cannot be shared with anything. `rho_override` does not touch it either. See [01](01-inputs.md). |

## The split of `Signal modelling`, and why it is now a no-op

Only z-mumu still fits a generator systematic:

| channel | comparison | status |
|---|---|---|
| μμ | powheg vs aMC@NLO (two NLO generators), two-sided, built inside a common 50 < m_LHE < 120 GeV window | **fitted**; 0.40 % of μ_Z |
| ττ | madgraph **LO** vs aMC@NLO, C_LO/C_NLO = 0.841 | **reported, not fitted** since v2.1 |
| ee | — | none |

z-tautau removed its version deliberately (`z-tautau/docs/07`, its `REVIEW.md` §3.4): an
uncertainty built as "nominal minus a worse model" is a placeholder, it double-counts
`QCDScale`/`PS_ISR`, and under its old name it would have been correlated in a MultiFit with the
μμ NLO-vs-NLO comparison. The split survives because the μμ side still needs it — and because a
channel that reinstates its own generator parameter would be handled without a code change:

```
f      = impact(SigModel) / impact(Signal modelling category)
s_gen  = s_category · f            ρ = 0
s_th   = s_category · sqrt(1 − f²) ρ = 1
```

which preserves the published category total exactly. `impact(SigModel)` comes from each channel's
NP ranking — zero for ττ and ee, 0.40 % for μμ — except for the μμ counting extraction, where only
the template's normalisation difference survives (0.28 %, `meta.sigmodel.C_ratio_powheg_over_nlo`).
The `sigmodel_correlated` variation is therefore **exactly null** (0.0 pb): there is nothing on the
other two sides to correlate it with. It is kept in the table to make that visible rather than
silent.

## Acceptance

Outside the μμ and ττ likelihoods; inside the ee one, so every ee entry here is zero.

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
+0.4 pb. The choice does not matter at this precision.

## What the resulting correlations are

| pair | ρ |
|---|---:|
| μμ, ee | **0.436** |
| μμ, ττ | 0.125 |
| ττ, ee | 0.087 |

ρ(μμ, ee) = 0.44 is by far the largest correlation in the combination and it is what makes the
two precise channels gain only 10.7 % rather than the 29 % two independent measurements of equal
precision would give. It is almost all luminosity: 22.7 × 18.8 pb of shared uncertainty out of a
35.4 × 39.6 pb product.

It also stays comfortably below the ratio of the two uncertainties (35.4 / 39.6 = 0.89), so both
channels get positive weights (μμ +0.596, ee +0.404). The ττ weight is −0.0003, negative for the
usual BLUE reason ([03](03-method.md)); its effect on the central value is 0.06 pb.

The bracketing variations `rho_zero` and `rho_one` span 9.9 pb, under a third of the total
uncertainty — so the correlation model is **not** the limiting assumption here. The z-ee signal
template normalisation is, at 99 pb.
