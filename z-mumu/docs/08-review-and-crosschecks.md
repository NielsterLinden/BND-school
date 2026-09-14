# Review and cross-checks

`scripts/xcheck_efficiency.py`, `scripts/mc_acceptance.py`

A review of the steps 1-6 measurement by a second person (Niels, 14 Sep 2026,
on `stbc-i2`), with every concern turned into a test that can be rerun. Nothing
in steps 1-6 was changed: the committed result is still what `run_all.py`
produces, and this page says what the tests found and what the numbers become
if the corrections are adopted.

## Summary

| # | Check | Verdict | Effect on sigma_fid |
|---|---|---|---:|
| 1 | Rerun steps 1-6 from scratch | reproduced exactly, every intermediate number | 0 |
| 2 | Input integrity (skim vs parent vs CERN catalogue) | complete and uncorrupted | 0 |
| 3 | Trigger efficiency (reference-trigger method) | **biased high by ~1%**, confirmed in simulation | +1.5% |
| 4 | ID / isolation tag-and-probe | **background in the failing probes** biases both low | -1.9% |
| 5 | Folding the efficiency over the signal | per-event 1/eps weighting needed | +0.5% |
| 6 | Luminosity | **brilcalc value lacks the normtag**, 0.63% low (all channels) | -0.6% |
| 7 | "A = 1 by construction" | not exact: 1.1% kinematic migration | +1.1% |
| 8 | Acceptance for the combination | computed from DY NLO and LO, with PDF/scale | -- |

The corrections partly cancel. Adopting all of them:

```
committed:  sigma_fid = 773.2 +/- 0.2 (stat) +/- 11.9 (syst) pb
revised:    sigma_fid = 776.9 +/- 0.2 (stat) +/- 14.8 (syst) pb     (dressed muons)
prediction: 799.6 pb  (DY amcatnloFXFX, normalised to 6077.22 pb)  -> measured / predicted = 0.972
```

The central value barely moves because corrections of 1-2% in both directions
happen to cancel; each one matters on its own. The uncertainty grows because
the tag-and-probe background (0.9%) and the method closure in simulation (0.8%)
are now included; the assigned 0.5% trigger systematic is replaced by measured
ones that are smaller (0.2% + 0.2%).

## 1. Reproduction

`run_all.py` in a clean copy, with a `.venv` built from `requirements.txt`,
gives 773.2 +/- 0.2 +/- 11.9 pb, the same cutflow, the same efficiencies to four
digits and the same 21,510,663 tag-and-probe pairs. About 60 s on `stbc-i2`.

## 2. Inputs

- The parent NanoAOD (`config.PARENT_DIR`) holds **323,952,013** events,
  exactly the catalogue number, of which **80,191,719** have >= 2 muons --
  exactly the skim. The skim is complete.
- All 152 parent files on dCache match the CERN catalogue in size and adler32
  (read from dCache's stored checksum).
- Step 2's own `_trigger_efficiency` and `_tag_and_probe`, run on the parent,
  return step 2's numbers exactly, and the step 1 selection on the parent gives
  exactly 10,777,373 events. The cross-checks below therefore look at the same
  events as the measurement.

**Operational note.** dCache NFS reads of the parent stalled repeatedly under
10-24 parallel readers (uninterruptible I/O; `dd` on an uncached file hanging),
although the files are intact. The cross-check scripts therefore stream from
CERN EOS over xrootd by default and run each file in its own subprocess
(`zmumu/batch.py`), resumable. opendata.cern.ch rate-limits HTTP range requests
(status 429), so use xrootd, not HTTP.

## 3. Trigger efficiency

**The concern.** The reference paths (`HLT_IsoMu27`, `HLT_Mu27`, `HLT_Mu50`,
`HLT_Mu45_eta2p1`) share the L1 seed and the L3 muon reconstruction with
IsoMu24. If IsoMu27 fired, IsoMu24 essentially always did, so the ratio is close
to one by construction and mostly measures the online isolation. The committed
plot shows exactly that: 1.000 at 26-50 GeV (IsoMu27 dominates the reference),
falling to 0.987 above 120 GeV (Mu50 dominates, which has no isolation). And
`docs/04` justifies ~0.998 with `1 - (1 - eps)^2`, which only holds if *both*
muons can fire -- the subleading muon is often below 24 GeV.

**The test in simulation**, where the truth is the fraction of offline-selected
events that fired:

| DY sample | truth | step-2 reference method | TrigObj tag-and-probe |
|---|---:|---:|---:|
| NLO | 0.9863 | 0.9967 (+1.05%) | 0.9878 (+0.15%) |
| LO  | 0.9868 | 0.9968 (+1.01%) | 0.9887 (+0.19%) |

**The measurement in data.** The parent NanoAOD has `TrigObj_*`. A muon is
matched if an HLT muon object with `filterBits & (2|8)` (Iso | IsoTkMu) and
`pT > 24` GeV lies within dR < 0.1. Validation on 1.6 M events: such an object
exists in 99.99% of events that fired the OR; in non-fired events it exists in
5%, from L1_SingleMu20-seeded paths sharing the Iso filter, so a match only
counts in fired events.

- Tag: selected muon, pT > 26, matched. Probe: selected muon, pT > 20, opposite
  sign, 70 < m < 110 GeV. Per-muon efficiency in (pT, eta), `config.TRIG_PT_BINS`.
- **Plateau per-muon efficiency 0.900** (G 0.8999, H 0.9002), with the known
  2016 eta structure (0.81 in the forward endcap, 0.89 at |eta| < 0.3, 0.93 in
  between): `output/plots/xcheck_trigger_tnp_vs_{pt,eta}.png`.
- Folded over the signal muons with 1 - (1 - eps1)(1 - eps2) and a 1/eps event
  weight: **event efficiency 0.9835**, against 0.9977 committed.
- **Model check:** with this map, the predicted fraction of triggered signal
  events in which *both* muons are matched is 0.7706; observed 0.7707. No
  muon is matched in only 0.004% of triggered signal events.
- Systematics: vetoing objects whose L1 seed is in [20, 22) GeV gives 0.9814
  (0.21%); the simulation closure above (0.16%). The assigned 0.5% is dropped.

## 4. Background in the ID / isolation tag-and-probe

**The concern.** Probes are loose (`isTracker | isGlobal`, pT > 20) and pass and
fail are simply counted in 70-110 GeV. Non-Z pairs land mostly in the fail
sample. `step2_tagprobe_mass.png` shows a flat continuum under the failing-ID
peak.

**What the spectra show** (`output/plots/xcheck_{id,iso}_fail_mass_lowpt.png`).
At 20-25 GeV the opposite-sign continuum under the failing-isolation peak is
roughly ten times the same-sign one -- heavy flavour (bb -> mu mu) is mostly
opposite-sign -- so same-sign subtraction under-subtracts. A template fit
(fail = s * pass shape + exponential, 60-120 GeV) removes more, but also absorbs
any extra width of the failing-probe peak. The two bracket the truth.

| probe pT [GeV] | 20-25 | 25-30 | 30-35 | 35-40 | 40-50 | 50-60 | 60-80 | 80-120 | 120-200 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| ID raw | 0.942 | 0.966 | 0.976 | 0.980 | 0.983 | 0.979 | 0.971 | 0.956 | 0.927 |
| ID same-sign subtracted | 0.974 | 0.978 | 0.982 | 0.984 | 0.985 | 0.984 | 0.984 | 0.985 | 0.986 |
| ID template fit | 0.984 | 0.984 | 0.984 | 0.984 | 0.985 | 0.985 | 0.986 | 0.987 | 0.989 |
| ID: background share of fail probes (fit) | 74% | 49% | 30% | 18% | 13% | 27% | 50% | 70% | 84% |
| iso raw | 0.820 | 0.877 | 0.919 | 0.951 | 0.975 | 0.982 | 0.984 | 0.984 | 0.986 |
| iso same-sign subtracted | 0.828 | 0.879 | 0.920 | 0.951 | 0.975 | 0.983 | 0.985 | 0.985 | 0.986 |
| iso template fit | 0.856 | 0.899 | 0.929 | 0.955 | 0.977 | 0.985 | 0.988 | 0.990 | 0.991 |

(Tags trigger-matched. The raw ID efficiency *falling* above 50 GeV is not
physics -- medium ID is flat in pT -- it is background, and it disappears
after subtraction.)

Also measured: requiring the tag to be trigger-matched lowers the raw
efficiencies by 0.1% (ID) and 0.2% (iso) -- without it, some events fired only
on the probe, which biases the probe towards passing.

**In simulation** (no background), step 2's tag-and-probe reproduces the
truth-matched efficiencies to 0.15% (ID) and 0.37% (iso) per muon in the LO
sample, 0.12% / 0.10% in NLO. So the method is sound; the background is the
problem.

**What is used.** Step 2's (pT, eta) maps are rescaled per pT bin by
corrected/raw, for both corrections, and the midpoint of the two folded
efficiencies is taken, with half the difference (**0.88%**) as a systematic,
plus the simulation closure (**0.81%**, twice the per-muon LO deviations in
quadrature). A proper fix is a simultaneous pass/fail fit per bin with a signal
shape that allows for the failing-probe resolution.

## 5. Folding the efficiency

Step 5 averages each map over the *observed* muons and multiplies the averages.
The selected sample is the *efficiency-weighted* one, so the unbiased estimator
is N_true = sum over events of 1/eps_event. With per-event
eps = id(mu1) iso(mu1) id(mu2) iso(mu2), the ID x iso efficiency goes from 0.8634
to 0.8589: **+0.52%** on sigma. (The trigger and ID x iso factors are folded
separately, which neglects their correlation -- second order.)

## 6. Luminosity

`config.LUMI_PB = 16290.713420` comes from `brilcalc lumi -c web --begin 278820
--end 284044 -i GRL.txt`, **without `--normtag`**, i.e. the online luminosity,
not the calibrated one the 1.2% uncertainty refers to. The CMS Open Data
luminosity record for 2016 ([recid 1059](https://opendata.cern.ch/record/1059))
gives the normtag-corrected values for exactly these lumisections:
Run2016G 7.653261 + Run2016H 8.740119 = **16.393381 fb^-1**. Summing its
per-lumisection table (`pp_2016lumibyls.csv`) over `datasets/GRL/GRL.txt`
reproduces this to the last digit. The value in use is **0.63% low**, which
inflates every channel's cross section by 0.63% -- z-ee and the combination
share the number. `config.LUMI_PB_NORMTAG` records the corrected value.

## 7. The fiducial volume is not "A = 1 by construction"

The efficiency corrections undo losses from ID, isolation, trigger and
reconstruction. They do not undo **migration across the volume boundary**:
the analysis cuts on reco muon pT and FSR-recovered mass, the volume is defined
on dressed generator muons (photons within dR < 0.1 added to the muon pT, which
a reco PF muon does not have), and resolution moves events across the pT and
mass thresholds.

In simulation: events with exactly two *prompt* reco muons passing only the
kinematic cuts (no ID, isolation or trigger), divided by the generator-level
fiducial events and by reco_eff^2, is **0.9891** (NLO) / 0.9901 (LO). The
efficiency-corrected yield is therefore 1.1% below the dressed fiducial
yield. Of the reco-selected events, 99.55% are inside the generator volume.

As a check of the whole chain in simulation, N_reco / (eps_method * N_fid) with
per-muon inclusive efficiencies is 0.982 with the step 2 methods, 0.991 with the
TrigObj trigger, and **1.002** after the migration correction (NLO).

Also: the fiducial table lists `mediumId` and `pfRelIso04_all < 0.15`, but
those are corrected for by the efficiencies, so they are not part of the volume
that is measured. The volume is kinematic: two muons, pT, eta, charge, mass.

## 8. Acceptance and the fiducial prediction

From the full DY samples (`scripts/mc_acceptance.py`), streamed from EOS:
71,839,442 NLO and 82,448,537 LO events. In both, the summed genWeight equals
the `Runs` genEventSumw, and the LHE mu mu fraction is 0.3338 (NLO) / 0.3337 (LO).

Generator-level volume: exactly two `GenDressedLepton` muons without tau
ancestor, pT > 20, |eta| < 2.4, leading pT > 26, opposite charge,
60 < m < 120 GeV, in an LHE Z -> mu mu event.

| | NLO (amcatnloFXFX) | LO (madgraphMLM) |
|---|---:|---:|
| predicted sigma_fid = 6077.22 pb x fiducial fraction | 799.6 pb | 825.4 pb |
| **A = sigma_fid / (6077.22/3)** | **0.3947** | 0.4074 |
| MC stat | 0.05% | 0.03% |
| PDF (NNPDF3.1 Hessian, LHA 325300) | 0.52% | 0.44% |
| alpha_s (+/-0.0015) | 0.03% | 0.06% |
| scale (7-point) | 0.32% | 0.82% |
| A w.r.t. LHE mu mu (m > 50) | 0.3942 | 0.4070 |
| A w.r.t. LHE mu mu, 60 < m < 120 GeV | 0.4092 | 0.4219 |

**LO and NLO differ by 3.2%** in A, far more than either sample's scale and
PDF uncertainties. NLO is the natural nominal choice; whether to assign the
difference as a modelling uncertainty is a decision for all three channels
together, since it is correlated between them.

**Convention.** The combination notebook expects sigma near 6077.22/3 pb,
i.e. sigma(Z/gamma* -> mu mu, m_LHE > 50 GeV). The first row of A matches that.
`handoff.md` had proposed a 60-120 GeV denominator, which gives a number 3.7%
larger and a cross section for a different mass range. The three channels must
use the same one.

## 9. Revised result

```
sigma_fid = N_sig / (eps_total * L * M)
          = 10,763,262 / (0.8544 * 16393.4 pb^-1 * 0.9891)
          = 776.9 +/- 0.2 (stat) +/- 14.8 (syst) pb
```

| Source | Relative |
|---|---:|
| luminosity | 1.200% |
| eff: ID/iso background method | 0.883% |
| eff: ID/iso tag-and-probe closure (MC) | 0.813% |
| eff: reconstruction (external) | 0.803% |
| eff: trigger (L1-band veto variation) | 0.208% |
| muon momentum scale | 0.200% |
| eff: trigger tag-and-probe closure (MC) | 0.155% |
| migration (NLO vs LO) | 0.101% |
| background estimate | 0.056% |
| statistical | 0.031% |
| others (T&P stat, FSR) | < 0.01% each |
| **total** | **1.91%** |

Against the NLO prediction (799.6 pb, uncertainty 0.6% from A plus the
unquoted normalisation uncertainty of 6077.22 pb): measured / predicted = 0.972.

## 10. For the combination

| Variable | Committed | Revised | Note |
|---|---:|---:|---|
| `n_obs` | 10,777,373 | 10,777,373 | |
| `n_bkg` | 14,110 | 14,110 | |
| `acc_eff` | 0.8546 (eps only) | **0.3336** | A (NLO) x eps x migration |
| `lumi_pb` | 16290.713 | **16393.381** | normtag; applies to every channel |
| sigma(Z -> mu mu, m > 50) | -- | **1968 pb** | target 2025.7 pb |

With `lumi_pb` left at 16290.713 and `acc_eff = 0.3336`, the notebook would
return 1981 pb instead.

Add `acceptance` (0.61% for NLO: PDF, alpha_s, scale, MC stat) to the
correlated sources, and agree on the NLO-vs-LO question (item 8).

## What was not checked, or remains open

- **Reconstruction efficiency** stays the external 0.996 +/- 0.004 per muon.
  Simulation truth is 0.9986; the data value still cannot be measured in NanoAOD.
- **ID / isolation background**: a proper simultaneous pass/fail fit per
  (pT, eta) bin would replace the bracket in item 4 and shrink its 0.9%.
- **Joint folding** of the trigger and ID x iso efficiencies per event.
- **L1 prefiring**, **Rochester corrections**, the Voigt fit model, R_OS/SS and
  the e-mu `k` factor: unchanged from `handoff.md`.
- **The 6077.22 pb normalisation** and its uncertainty were not checked against
  the CMS cross-section TWiki (the z-ee handoff already flags this).
- Smaller documentation inconsistencies: `handoff.md` quotes chi2/ndf as
  "908/115" (908 is chi2/ndf; chi2 is 104,469), `docs/06` says "~6"; the step 4
  docstring says the exponential bends *down*, `handoff.md` says *up*.

## How to reproduce

```bash
cd z-mumu
python3 -m venv .venv && .venv/bin/pip install -r requirements.txt   # includes xrootd
.venv/bin/python run_all.py                                          # steps 1-6, ~1 min

# needs the cernopendata-client listings in filelists/:
#   cernopendata-client get-file-locations --recid <id> --protocol xrootd --verbose \
#       > filelists/<name>_<recid>.sizes.tsv
.venv/bin/python scripts/xcheck_efficiency.py        # 152 files from EOS, ~10 min on 8 workers
.venv/bin/python scripts/mc_acceptance.py            # NLO + LO from EOS, ~15 min on 8 workers
```

Both scripts resume after an interruption and accept `--summarise-only` to
redo the numbers and plots from the stored results.
