| source | CMS-SMP-20-004 [%] | ours [%] | status | how we do it | how they do it |
|---|---:|---:|---|---|---|
| Luminosity | 2.30 | 1.20 | ours better | 1.2 % for the 2016 legacy calibration (CMS-LUM-17-003), normtag luminosity 16393 pb^-1 | 2.3 % for the 2017 low-pileup runs (LUM-17-004) |
| μR and μF scales | 0.66 | 0.31 | similar | 7-point envelope on A(60-120) 0.29 % ⊕ in the fit on C 0.10 %; a CMS-like 25/25 GeV volume gives 0.38 % | largest acceptance shift of the 6 non-extreme variations |
| PDF + αs | 0.43 | 0.51 | similar | NNPDF3.1 Hessian eigenvectors and αs ± 0.0015 on A(60-120) and in the fit on C | NNPDF3.1 prescription on the acceptance |
| Boson pT (resummation) and generator on A | 0.12 | 0.38 | measured now | new: A with the pT(Z) spectrum reweighted to our measured pT(μμ) (+0.38 %; PS ISR 0.07 %), powheg vs aMC@NLO (0.04 %), PS FSR (0.03 %) | DYTURBO NNLO+NNLL vs aMC@NLO on A, together with PYTHIA vs PHOTOS |
| Generator and parton shower on C | — | 0.42 | ours only | in the fit: powheg vs aMC@NLO lineshape (SigModel), PS ISR/FSR weights, all on C | no generator or shower term on the efficiency in the table |
| QED FSR model (PHOTOS vs PYTHIA) | — | a 0.12 / b 0.00 / **c 0.10** | partial / eyeballed | not measurable here: estimated (three options) | inside Resum. + FSR (0.12 %) |
| L1 trigger prefiring | 0.34 | 0.52 | similar | NanoAOD L1PreFiringWeight Up/Dn (muon + ECAL maps, stat ⊕ 20 %) | prefiring maps, stat ⊕ ±20 % of the probability |
| Lepton ID, isolation, trigger efficiency | 0.32 | 0.39 | similar | own T&P fits, stat ⊕ syst (bkg shape, fit range, tag, powheg template, closure) shifted coherently over cells (conservative for the stat part); no charge binning (second order for Z: SF(μ+)SF(μ-)) | T&P in charge × pT × η bins, stat uncorrelated per bin; syst: fit model, PHOTOS FSR in templates, tag |
| Muon reconstruction (tracking, stand-alone) | — | 0.26 | measured now | new: T&P on the unskimmed NanoAOD (stand-alone and isolated-track probes), SF per muon 1.0001 ± 0.0013; was 1 ± 0.004 assigned (0.78 % on μ) | T&P, inside Efficiency (stat)/(syst) |
| QCD multijet / non-prompt | 0.16 | 0.04 | similar | fake factor (anti-isolated muons), stat + 30 % method | misidentification factor (iso/anti-iso), stat + ±10/20 % prompt contamination |
| MC statistics | 0.16 | 0.36 | similar | Barlow-Beeston lite gammas per bin ⊕ MC stat on A | Barlow-Beeston |
| EW + tt̄ cross sections | 0.04 | 0.07 | similar | 6 % tt̄, 10 % tW/WW/WZ/ZZ, 5 % Z→ττ | 5 % top and diboson; τ backgrounds via theory ratios |
| Hadronic recoil calibration | 0.05 | 0.00 | n/a | no p_T^miss in the Z selection | W recoil calibration, enters Z via the joint W+Z fit |
| Muon momentum scale and resolution | — | 0.35 | similar | Z-peak calibration per |η| (Rochester corrections not obtainable), shape NPs | Rochester corrections; shape NP, 'small', not tabulated |
| Pileup | — | 0.15 | ours only | N_PV-matched profile, σ_minbias ± 4.6 % | low-pileup runs (<μ> ≈ 3): not needed |

| option | stat ⊕ syst [%] | total [%] | total [pb] | ours − CMS-SMP-20-004 [pb] | pull (uncorrelated) | pull (theory correlated) |
|---|---:|---:|---:|---:|---:|---:|
| 3a: published value | 1.047 | 1.59 | 30.7 | -21.0 | -0.36σ (±57.5) | -0.38σ (±54.7) |
| 3b: set to zero | 1.040 | 1.59 | 30.7 | -21.0 | -0.36σ (±57.5) | -0.38σ (±54.7) |
| 3c: our estimate | 1.045 | 1.59 | 30.7 | -21.0 | -0.36σ (±57.5) | -0.38σ (±54.7) |
