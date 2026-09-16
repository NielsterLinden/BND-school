# Chapter brief — Z→μμ (4_zmumu)

Status: approved (2026-09-15, deck owner) — **worked example**, written after the fact from the
decisions recorded in `../03-sections-storyboard.md` and the docstring of
`scenes/s4_zmumu_story.py`. Other chapters should not copy its story, only see how a brief reads.

## 1. Scope
- Section folder / number: `4_zmumu` (4); chained after the delivered 4-01 `mumu_process`, 4-02 `mumu_detector`
- Time budget: ≈ 4 min in the talk → 8 clips of 8–13 s (the speaker talks over held frames)
- Audience: PhD students, mostly outside CMS
- **The one sentence**: the raw simulation is ~6 % off; the fakes are tiny (0.037 %) and cannot
  explain that, the measured corrections do, and the fit gives σ(60–120) = 1931 ± 33 pb.

## 2. Story
- Spine nodes: comparison, backgrounds (fakes), corrections, fit. Selection and files are left to section 2.
- Method in depth: the fake factor as "1 in 10 passes" (f = N_tight/N_anti = 1/9), `z-mumu/docs/13-fake-factor.md`
- Left out: tag-and-probe mechanics (section 2 shows them), momentum calibration, the counting cross-check
- Tone: real data from the first clip of the chain on

## 3. Numbers on screen
Result version: z-mumu v2 shape fit, 12 × 5 GeV (`z-mumu/handoff.md`, `combination/result.md`, 15 Sep 2026).

| symbol | printed | source | frozen in |
|--------|---------|--------|-----------|
| event | Run 278969, Event 1437904187 | skim, SR | `zmumu_events.json` |
| N | 10,378,567 | SR data, 60 bins | `zmumu_sr_stack.json` `data.total` |
| data/pred (raw) | 0.944 | raw stack | `zmumu_sr_stack.json` |
| N_fake | 3870 ± 80 | `z-mumu/RESULTS_v2.md` | `zmumu_fakes.json` `yields.sr_fakes` |
| corrections | ⟨w_PU⟩ 0.994, ⟨w_L1⟩ 0.980, SF_ID 0.980, SF_iso 1.006, ε_trig 0.907/0.923 | docs/11, 12 | `zmumu_corrections.json` |
| μ_Z | 0.988 ± 0.016 | fit result | `zmumu_fit.json` |
| σ_fid | 790.1 ± 0.2_stat ± 8.5_syst ± 9.6_lumi pb | handoff | `zmumu_fit.json` `sigma_fid` |
| σ(60–120) | 1931 ± 33 pb vs 1953.9 pb | handoff | `zmumu_fit.json` |

## 4. Entry and exit
- Opens on: the final frame of 4-02 `mumu_detector`; its schematic tracks fade out
- Event display: four real SR events, ten real same-sign pairs
- Ends on: σ(60–120) point beside the purple prediction line on a 1850–2050 pb axis

## 5. Clips
| letter | clip | the one idea |
|--------|------|--------------|
| a | `mumu_a_event` | a real event becomes the first entry of a log m_μμ plot |
| b | `mumu_b_rain` | the plot fills to 10.4 M events |
| c | `mumu_c_stack` | the uncorrected simulation is ~6 % high |
| d | `mumu_d_control` | the same-sign control region |
| e | `mumu_e_tenpairs` | 1 in 10 passes → f = 1/9, the real FF map |
| f | `mumu_f_transfer` | the fakes slide into the stack; the ratio does not move |
| g | `mumu_g_corrections` | the corrections walk the prediction down; data never moves |
| h | `mumu_h_fit` | rebin, pulls, μ_Z, σ beside the prediction |

## 6. Look
- Layout: slice parked left (0.42 at (−4.6, −0.7)), plot right, one-row key — the defaults of 06 §B3 come from here
- Motifs: event → entry → rain, stack + ratio, control-region box, tally, value grid, pull plot
- Colours: gold fills, `shade(gold, .30)` lines; fakes pale slate; theory purple

## 7. Deviations
- Log axis from 10^1.5 so the fakes are a visible wedge; h's y axis per 5 GeV

## 8. Process
- Critique by the deck owner on 480p previews; reviews at delivery only
