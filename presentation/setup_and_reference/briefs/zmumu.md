# Chapter brief — Z→μμ (4_zmumu)

Status: approved (2026-09-17, deck owner) — **final pass** on the frozen result. Replaces the brief of
15 Sep (fake-factor story) and the re-cut of 16 Sep (five tag-and-probe clips). Other chapters should
not copy its story, only see how a brief reads.

## 1. Scope
- Section folder / number: `4_zmumu` (4); 21 clips, **clip number = play order** (renumbered on 17 Sep, see §7)
- Time budget: ≈ 5 min in the talk; one idea per clip, the speaker talks over the held end frames
- Audience: PhD students, mostly outside CMS
- **The one sentence**: the uncorrected prediction is ~6 % off; the Z peak itself tells us which muons are
  real, so we measure how often a real muon survives a cut, in data and in simulation, and correct the
  simulation by the ratio; with the other corrections the prediction matches, and the fit gives
  σ(60–120) = 1931 ± 30 pb.

## 2. Story
- Spine nodes: event → histogram → comparison → corrections (tag and probe first) → fit → cross sections.
- Method in depth: **tag and probe**, `z-mumu/docs/17-how-tag-and-probe-works.md` (and docs/12): tag,
  probe, the Z mass as the label; why a probe passes or fails tight ID (muon stations reached); pass / fail
  spectra of one cell; the identical procedure on simulation (2 × 2: pass | fail × data / simulation);
  a scan over the ten barrel p_T cells fills ε_ID(p_T); SF = ε_data/ε_sim below it; the SF points become a
  table, a scan over |η| (side view of the detector) completes it; the table is applied event by event.
- Said, not shown as its own clip: the background under the failing probes (drawn as a grey band under
  the fitted peak in every panel, so ε is the fitted number); systematic variants; the trigger as a count.
- Left out: the fake factor (Z→ττ tells it; here the fakes are a wedge in the prediction from its first
  appearance), momentum calibration, fit stability.
- Tone: real data from the first clip of the chain on. The data never move; only the prediction does.

## 3. Numbers on screen
Result version: **z-mumu frozen result, 17 Sep 2026** (`FREEZE.md`, tag `zmumu-freeze-2026-09-17`,
`z-mumu/handoff.md` "Frozen result"): v2 shape fit, 12 × 5 GeV, measured reconstruction SF, 60–120 acceptance.

| symbol | printed | source | frozen in |
|--------|---------|--------|-----------|
| events | Run 278969 Event 1437904187 (first entry); four more real SR events, masses 86.73, 88.86, 89.70, 94.61 GeV | skim, SR | `zmumu_events.json` `sr.events[0, 14, 4, 8, 12]` |
| N | 10,378,567 | SR data, 60 bins | `zmumu_sr_stack.json` `data.total` |
| data/pred ladder | 0.944 → 0.971 (ID table) → 0.968 (isolation, trigger, reconstruction) → 0.974 (pile-up) → 0.994 (prefiring) | the frozen weights in the order told here; `final` == the fit input | `zmumu_sr_stack.json` `recut.totals` |
| row factors | muon efficiency × 0.976, extra collisions × 0.994, trigger fired early × 0.980, prediction × 0.950; identification 0.972, isolation 1.007, trigger 0.996, reconstruction 1.0002 | mean event weights | `recut.step_factors`, `recut.raw_weighted_means`, `zmumu_tnp.json` `reco` |
| tag-and-probe pairs | run 278969: event 351671459 (m = 90.45 GeV, probe 3 muon stations, passes), event 351264945 (tracker-only probe, 1 station, fails) | first data skim file | `zmumu_tnp.json` `events` |
| cell | 40 < p_T < 45 GeV, \|η\| < 0.9: N_pass 1,999,122, N_fail 86,537, ε_data 0.9585, ε_sim 0.9711, SF 0.987 | docs/12, docs/17 | `zmumu_tnp.json` `cell`, `cells_barrel` |
| ε_ID(p_T), SF(p_T), 10 × 4 SF map (mean 0.980) | per cell, all four \|η\| columns | `tnp_result.json` | `zmumu_tnp.json` `id` |
| simulated events | three real DY events, w = SF·SF = 0.974, 0.961, 0.971 | first DY skim file | `zmumu_tnp.json` `apply_examples` |
| μ_Z | 0.988 ± 0.014 | fit result | `zmumu_fit.json` `poi` |
| σ_fid, counting | (10,378,567 − 68,799) / (0.792 × 16.39 fb⁻¹) = 794.5 pb | fit meta | `zmumu_fit.json` `counting`, `C`, `lumi_pb` |
| σ_fid, fit | 790.2 ± 0.2_stat ± 6.1_syst ± 9.6_lumi pb | handoff | `zmumu_fit.json` `sigma_fid` |
| A | 0.409 (1/A = 2.44) | `zmumu/acceptance.py` | `zmumu_fit.json` `sigma_60_120.A` |
| σ(60–120) | 1931 ± 30 pb vs 1953.9 pb | handoff | `zmumu_fit.json` `sigma_60_120` |

## 4. Entry and exit
- Opens on: 4-01 `mumu_process` (unchanged); 4-02 `mumu_into_detector` shrinks the diagram into the
  interaction point of the slice, then the two muons appear; the story chain opens on its last frame.
- Ends on: σ(60–120) point beside the purple prediction line on a 1850–2050 pb axis (unchanged; the object
  section 6 collects).

## 5. Clips (number = play order)
| # | clip | the one idea |
|---|------|--------------|
| 4-01 | `mumu_process` | the Drell-Yan diagram, μμ flavour |
| 4-02 | `mumu_into_detector` | the collision happens *there*: the diagram shrinks into the detector, two muons come out |
| 4-03 | `mumu_event_to_mass` | a real event; its two p_T values become one number, m_μμ |
| 4-04 | `mumu_first_entry` | that number is the first entry: a block on a count axis |
| 4-05 | `mumu_fill` | four more events (two share a bin), then 10.4 M: the axis goes logarithmic, blocks rain over the whole plot |
| 4-06 | `mumu_prediction` | the uncorrected prediction slides in under the data |
| 4-07 | `mumu_ratio` | data/pred = 0.944 |
| 4-08 | `mumu_tnp_tag_probe` | a tag, a probe, and the Z mass as the label: the probe is a real muon |
| 4-09 | `mumu_tnp_pass_fail` | the question (tight ID) and why one probe passes and another fails; pass / fail spectra; ε_data |
| 4-10 | `mumu_tnp_grid` | the same on simulation: pass \| fail × data / simulation |
| 4-11 | `mumu_tnp_scan` | scan p_T: the spectra change, ε_ID(p_T) fills point by point |
| 4-12 | `mumu_tnp_scale_factor` | SF = ε_data / ε_sim |
| 4-13 | `mumu_tnp_table` | the SF points become a table; the \|η\| scan completes it |
| 4-14 | `mumu_tnp_apply` | the table weights every simulated event: 0.944 → 0.971 → 0.968 |
| 4-15 | `mumu_corr_pileup` | extra collisions: → 0.974 |
| 4-16 | `mumu_corr_prefiring` | the trigger that fired one crossing early: → 0.994 |
| 4-17 | `mumu_corrected` | prediction × 0.950: the final data vs prediction, the fit input |
| 4-18 | `mumu_rebin` | 60 × 1 GeV → 12 × 5 GeV |
| 4-19 | `mumu_fit` | the fit: pulls, μ_Z = 0.988 ± 0.014, ratio flat |
| 4-20 | `mumu_sigma_fid` | σ_fid = (N − B)/(C·L): counting 794.5 pb, the fit 790.2 ± 0.2 ± 6.1 ± 9.6 pb |
| 4-21 | `mumu_sigma_total` | the acceptance A = 0.409: σ_fid/A = 1931 ± 30 pb beside the prediction |

## 6. Look
- Layout: slice parked left (0.42 at (−4.6, −0.7)), plot right, one-row key — the defaults of 06 §B3 come from here.
  Tag and probe: tag solid gold with a cyan trigger ring, probe dashed grey until the mass confirms it;
  pass panels twice as high as fail panels, each on its own linear scale (the fail peak is ~20× smaller).
- Motifs: event → entry block → rain of blocks, stack + ratio, 2 × 2 panel grid, scan band, value grid,
  side view of the detector (`eta_fan`) for |η| and for the acceptance, pull plot.
- Colours: gold fills, `shade(gold, .30)` lines; simulation outlines gold; fit model and SF points method cyan;
  fitted background grey; fakes pale slate; theory purple.

## 7. Deviations (authorised by the deck owner, 17 Sep 2026)
- **Renumbered.** The rule "clip numbers are never renumbered" (02-deliverables-and-naming.md) is set aside for
  this section once: 4-01…4-21 in play order. The superseded MP4s are in
  `clips/4_zmumu/00_archive/pre_recut_2026-09-17/`, their registry rows in `clips/CLIPLIST_4_zmumu_pre_recut.tsv`.
  4-01 `mumu_process` kept its file (`_v2`) and number.
- **Words on screen.** "pass" / "fail" under the panels, "data" / "simulation" beside the rows, and the
  correction rows in plain words ("muon efficiency", "extra collisions", "trigger fired early", with the
  technical name small below) are an exception to "no narrative text in a clip".
- The first entries use a linear count axis (a count of 1 is below the log floor); it becomes the log axis
  (from 10^1.5, so the fakes are a visible wedge) when the rain starts; 4-18 on: y axis per 5 GeV.
- The ten pairs of the acceptance picture and the side view are schematic; A is the frozen number.

## 8. Process
- Designed up front, one 480p pass per scene to catch crashes, one delivery render; critique by the deck owner
  on the delivered clips. Reviews (physics, style, continuity) only on request.
