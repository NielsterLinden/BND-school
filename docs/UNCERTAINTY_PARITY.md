# Uncertainty parity with a published measurement

How a channel makes its uncertainty budget comparable, row by row, with a published measurement of the
same quantity, so that "we agree with CMS at n σ" means something. Worked out on Z → μμ against
CMS-SMP-20-004 (16 Sep 2026): `z-mumu/docs/16-uncertainties-vs-cms.md`,
`z-mumu/scripts/v2_7_reco_tnp.py`, `v2_8_theory_acceptance.py`, `v2_9_cms_parity.py`, `z-mumu/zmumu/acceptance.py`.
The agent that runs this procedure for another channel is `.claude/agents/uncertainty-parity-auditor.md`;
ready-to-paste prompts for z-ee and z-tautau are at the end of this file.

> **The channel results were frozen on 17 Sep 2026.** Z → μμ went through this procedure and its
> measurements were promoted into the frozen result (tag `zmumu-freeze-2026-09-17`). A run for z-ee or
> z-tautau after that date is a study next to the frozen result: it produces tagged variants and a
> comparison, and never overwrites the channel's result files or the inputs the combination reads.

## Why

A published measurement and ours can disagree because the physics differs, or because one budget has a
row the other lacks. The second kind is not a result. Before quoting a pull, every row of the published
table must be accounted for in ours: done the same way, done better, done partially, or not done. The
last two are closed by measuring what can be measured and by an explicit, bracketed estimate for the rest.

## The procedure

1. **Get the published table.** Download the paper (`curl -sSL -o paper.pdf https://arxiv.org/pdf/<id>`,
   `pdftotext -layout`), and check HEPData (`https://www.hepdata.net/record/<ins>?format=json`) and the
   collaboration's public-results page for per-channel or additional tables. Record the table number,
   the column, fiducial vs total, and whether the channel is fitted alone or combined with others
   (CMS-SMP-20-004 fits e and μ together and publishes no per-flavour breakdown). Read the systematics
   *text* too: sources that are "small" or folded into another row are listed only there.
2. **Read our budget from the source, not from memory.** The fit result JSON (grouped impacts and the
   per-NP ranking), the TRExFitter config (what each NP is), the acceptance code (which denominator!),
   the scale-factor and fake-factor docs. Note for every NP how it is derived and what it is correlated with.
3. **Build the cross-reference table**, one row per published source (split a row when ours splits it,
   keep the published number on the first of the split rows). Status of each row:

   | status | meaning |
   |---|---|
   | `similar` | same kind of evaluation, comparable scope |
   | `better` | ours is measured where theirs is assigned, or ours is more complete / more precise |
   | `partial` | we assign instead of measure, cover only part of the source, or estimate it |
   | `missing` | not in our budget at all |
   | `not-applicable` | the source cannot affect our observable (say why) |
   | `ours-only` | we carry it, they do not need it (say why) |
   | `measured-now` | was `partial`/`missing`, closed in this pass by a measurement |

4. **For every `partial` or `missing` row decide: measure, or estimate.** Measure whenever the data or
   simulation in reach can do it within the session (a tag-and-probe on the unskimmed NanoAOD, a
   generator-level pass over the parent MC, a data/MC shape comparison that constrains a model, a
   different denominator). "Not in the skim" is not a reason: the parents are on dCache/EOS.
   Look for data constraints on the model before calling something theory-only (Z → μμ: the measured
   pT(μμ) spectrum sized the resummation row; the recovered-FSR-photon rate sized the QED FSR estimate).
5. **Estimate only what cannot be measured**, and for each such source write down
   1. *why it cannot be measured* (the missing sample, collection, tool or computing, concretely);
   2. *why the published number cannot be cited* (merged with another source, other channel or phase
      space, other definition, internal note);
   3. three options, each carried to the end: **a** the published number, **b** zero, **c** our
      back-of-the-envelope estimate as (size of the effect we can compute) × (plausible model difference,
      ideally sized by data), with the derivation written out.
6. **Add the measurements to the result without touching the channel's baseline.** A fit variant with a
   `--tag` (Z → μμ did this first as the tagged fit `zmumu_recosf`), acceptance items outside the fit as
   the channel already does. The combination (a TRExFitter MultiFit, `combination/combLieke/`) reads each
   channel's baseline config, fit inputs and acceptance metadata (`combLieke/config/channels.json`);
   promoting a variant is the channel's decision, taken with the combination group. For Z → μμ it was
   promoted on 17 Sep 2026: the reconstruction SF became the default of `v2_5_fit.py`, the acceptance rows
   moved into `zmumu/acceptance.py` (the `acceptance` block of the fit result and of `zmumu.root.meta.json`),
   and each new row got a combination parameter (`Acc_PTZ_mumu`, `Acc_Generator_mumu`, `Acc_PS_FSR`,
   `Acc_QEDFSR_mumu`).
7. **Write the parity file and run the shared module:**
   `python fitting/uncertainty_parity.py <parity.json> --plot-dir <plots> --prefix <name>` (or call
   `fitting.uncertainty_parity` from the channel script). It gives the total for each option, the pull
   with and without the shared theory rows correlated (`rho_reference`), one budget-plus-comparison
   figure per option, an overview of the three options, and the markdown table.
8. **Document** in the channel's `docs/`: the explicit lists (similar / better / partial / missing /
   n.a. / ours-only), what was measured and how, each estimate with its three options, the comparison,
   and what would change the conclusion. Point the channel `handoff.md` and `CLAUDE.md` step table to it.

## The parity file

Schema in the docstring of `fitting/uncertainty_parity.py`. Rules that matter:

* percentages are of the respective central value; `reference_pct: null` when the published table has no
  row for it (its number sits inside another row: say which in `reference_how`);
* `group`: `luminosity` rows are kept apart from `syst` in the summary; everything else is systematic;
* `rho_reference`: 1 for theory rows computed with the same generator and PDF set as the reference
  (scales, PDF), 0 otherwise; the module quotes the pull both ways;
* an `estimate` row must carry `why_not_measured`, `why_not_cited`, `options` {a, b, c} and `derivation`.

## Pitfalls met on Z → μμ

* **The acceptance denominator.** The acceptance theory uncertainties were computed against
  m > 50 GeV while the cross section is quoted in 60–120 GeV. Recompute on the quoted denominator.
* **"Born level" from LHE particles is generator-dependent.** LHE leptons precede the shower recoil;
  powheg and aMC@NLO agree on the dressed acceptance to 0.04 % but differ by 1.4 % in an LHE-lepton
  "Born" volume. Use dressed leptons, or the last pre-FSR copies, never LHE kinematics for cuts.
* **"Not measurable in NanoAOD" was wrong** for the muon reconstruction efficiency: NanoAODv9 keeps
  stand-alone muons (`Muon_isStandalone`) and isolated tracks (`IsoTrack`, cleaned of loose muons).
  Check the branch list of a parent file before accepting such a statement.
* **Merged published rows.** "Resum. + FSR", "Efficiency (stat)" etc. merge sources and lepton flavours.
  Say so, and do not split them by guesswork: option a takes the whole row.
* **Correlation with the reference.** Theory rows computed from the same generator and PDF set are
  largely common to both measurements; a pull that ignores this overstates the tension.
* **Keep dCache readers ≤ 8 in total** (a second job streams from EOS instead).

## Prompts for the other channels

Spawn the agent (`.claude/agents/uncertainty-parity-auditor.md`) from the repository root with one of
these. They are complete: the agent reads everything else itself. (`combination/docs/` was removed from
the working tree on 16 Sep 2026 when the combination became a MultiFit; its comparison page with the
per-channel gap list is still in git, hence the `git show`.)

### Z → ee

> Run the uncertainty-parity procedure (`docs/UNCERTAINTY_PARITY.md`) for **z-ee**. Reference
> measurements: **CMS-SMP-20-004** (arXiv:2408.03744, Tables 4/7/10/13; e and μ fitted together, the
> electron-specific rows are prefiring (2017 ECAL), Efficiency (stat/syst) and the electron momentum
> scale in the text) and, because it publishes Z → ee alone, **ATLAS 1603.09222** (Table 2 and the
> per-channel systematic table; 66–116 GeV, scale by 1.0143 to 60–120 as in
> `git show 17d497c:combination/docs/05-vs-published.md`). Our inputs: `z-ee/fit.config`, `z-ee/Zee_fit.tar.gz`,
> `z-ee/handoff.md`, how the combination reads z-ee (`combination/combLieke/mf/ee_input.py`,
> `combination/combLieke/README.md` "The ee channel"), the gaps listed in
> `git show 17d497c:combination/docs/05-vs-published.md` §2 (electron trigger SF, electron energy
> scale/resolution, charge misidentification, FSR, multijet) and the open ECAL-gap problem in the
> repository guide `.claude/CLAUDE.md` (no barrel–endcap gap veto, placeholder ID SF there, `ElectronID` ±5.9 %). Measure what can be
> measured on the Open Data (electron trigger and reco/ID efficiencies by tag-and-probe on the `Electron`
> dataset parents, energy scale/resolution from the Z peak, charge flips from same-sign Z → ee,
> multijet from anti-isolated or same-sign data); estimate the rest with the three options. Write
> `z-ee/docs/uncertainties-vs-published.md`, a parity file and plots via `fitting/uncertainty_parity.py`;
> do not change the z-ee baseline — produce tagged variants. Commit locally, do not push.

### Z → τhτh

> Run the uncertainty-parity procedure (`docs/UNCERTAINTY_PARITY.md`) for **z-tautau**. Reference
> measurement: **CMS Z → ττ at 13 TeV, arXiv:1801.03535** (five ττ final states fitted together, τh ID
> efficiency constrained to 2.2 % in the fit; 60–120 GeV Born mass). Download it and its HEPData, take
> the systematic table of the inclusive cross section, and separate what is specific to τhτh from what
> the eτh/μτh/eμ/μμ final states buy them (the τh ID constraint: `git show 17d497c:combination/docs/05-vs-published.md`
> §2 "The one that should embarrass us"). Our inputs: `z-tautau/CLAUDE.md`, `z-tautau/docs/`,
> `z-tautau/output/results.json` (its `nominal_variant`), `z-tautau/fit/`, `z-tautau/review/REVIEW_v3.md`,
> `z-tautau/external/*.json` (TauPOG SFs). Measure what can be measured (e.g. the τh ID efficiency from a
> μτh tag-and-probe on `SingleMuon` parents, τh energy scale from the visible-mass peak, trigger
> efficiency with the orthogonal-trigger method, the fake-factor closure), estimate the rest with the
> three options. Write `z-tautau/docs/uncertainties-vs-published.md`, a parity file and plots via
> `fitting/uncertainty_parity.py`; do not change the z-tautau baseline — produce a tagged variant
> (`BND_TAUTAU_WP`-style redirection, never overwrite `output/results.json`). Mind the `/data/atlas`
> quota noted in `z-tautau/CLAUDE.md`. Commit locally, do not push.
