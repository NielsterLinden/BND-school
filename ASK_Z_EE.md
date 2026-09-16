# What z-ee still has to change

For the Z → e⁺e⁻ group, so the channel lines up with Z → μμ and Z → τhτh and can go into the
combination without a caveat. Plain list, biggest first. Details and numbers:
`combination/docs/01-inputs.md` and `combination/docs/05-vs-published.md`.

**Already done (16 Sep 2026):** reference cross section 1954.1 pb ✓, negative bins ✓.

---

## Must fix — this one changes the answer by 12 %

1. **Don't let the QCD-scale uncertainty change the total yield.**
   `h_mass_ScaleUp/Down` and `h_mass_PDFUp/Down` are currently ±5.9 % and ±0.3 % on the
   *normalisation*. That is the uncertainty on the predicted cross section, and it is the same
   thing `mu_signal` measures — so the fit cannot tell them apart. It pulls `QCDScale` to −1.9σ,
   shrinks the prediction by 11 %, and `mu_signal` comes out 8 % too high.
   **Fix:** rescale each LHE replica histogram to the nominal total before writing it, so only the
   *shape* varies. One line in the notebook.

## Should fix — missing systematics the other channels and CMS/ATLAS all have

2. **Electron trigger efficiency.** No scale factor and no uncertainty for
   `HLT_Ele27_WPTight_Gsf`. Measure it by tag-and-probe like z-mumu does, or assign a number.
3. **Electron energy scale and resolution.** No correction, no nuisance parameter. This is the one
   systematic that moves the shape of the peak you are fitting.
4. **FSR.** No photon recovery and no FSR uncertainty. Electrons radiate more than muons, so this
   matters more here than in z-mumu, which does both.
5. **Charge misidentification.** You require opposite sign; electron charge flips are a few per
   mille in the barrel and up to ~1 % in the endcap. ATLAS quotes 0.1 % for this.
6. **Fake / multijet background.** No data-driven estimate and no uncertainty. CMS quotes 0.14 %.

## Should fix — fit quality

7. **Use coarser bins.** 60 × 1 GeV over-constrains the shape parameters: `Pileup` is pulled
   −2.9σ, `L1Prefiring` +3.6σ and `ElectronID` is squeezed to 8 % of its input width. z-mumu had
   the same problem and fixed it by fitting in 5 GeV bins.
8. **Run the goodness of fit.** Use `trex-fitter ...s` (saturated model). Yours is the only
   channel with no GoF, so nobody can say whether the fit describes the data.
9. **Fix the `Wjets` template.** 24 % MC statistics and ~30 bins that were negative before being
   floored to 10⁻⁶. Flooring stops TRExFitter crashing; it does not make the template usable.
   Smooth it, or merge W+jets into the other backgrounds.

## Should fix — conventions (`fitting/CONVENTIONS.md` §1–§4)

10. **Rename things** so the MultiFit can correlate by name:
    - POI `mu_signal` → **`mu_Z`**
    - `DY_ee` → **`DYee`**, and give it `Type: SIGNAL` (it is currently BACKGROUND)
    - region `SR` → **`ee_SR`** (names must be unique across channels)
    - `LUMI` → **`Lumi`**; `ElectronRECO` → **`ElectronReco`**
    - `XS_VV` → split into **`XS_WW`, `XS_WZ`, `XS_ZZ`** (that is how μμ and ττ name them)
11. **Write one input file** `z-ee/fit/fitinputs/zee.root` with histograms named
    `ee_SR__<sample>` and `ee_SR__<sample>__<syst>Up|Down`.
12. **Publish a results file** `z-ee/fit/results/zee_fit_result.json` from `fitting/run_trex.py`
    (μ_Z, MINOS errors, grouped impacts, σ_fid, σ(60–120), n_obs, n_bkg, A, C, L). Today the
    combination reverse-engineers all of that out of `Zee_fit.tar.gz` with
    `combination/tools/extract_zee.py` — delete that script once the JSON exists.
13. **Quote σ(60 < m_LHE < 120 GeV).** `handoff.md` still says 2124.6 pb, which is μ × 6077.22/3,
    i.e. the m > 50 GeV normalisation. On the combination's footing it is μ × 1954.1 pb.
14. **Make `HistoPath` repo-relative.** It points at `/project/atlas/users/lvdurenw/...`, so
    nobody else can rerun the fit.

---

### What it is worth

| | σ(Z → ee) | combined | χ²/ndf |
|---|---|---|---|
| as published now | 2054 ± 40 pb | 1980 ± 32 pb | 9.67/2 (p = 0.008) |
| item 1 alone, undone by hand | 1833 pb | 1882 ± 30 pb | 7.39/2 (p = 0.025) |
| without z-ee | — | 1931 ± 35 pb | 0.37/1 (p = 0.54) |

Item 1 is the whole story. Items 2–14 are what a reviewer will ask for next.
