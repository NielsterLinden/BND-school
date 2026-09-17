# Chapter brief — theory (1_theory)

Status: approved (2026-09-17, deck owner: "whatever you think is best, make the animations" after the storyline page
https://claude.ai/artifact/HUzg3pjsKLBETaEADph4dY). The owner's steer: the section is about the **Z cross section**,
which the final plot shows, not the Z mass; lepton universality in; open on the SM but get deep early.

## 1. Scope
- Section folder / number: `1_theory` (1); scene `scenes/s1_theory_story.py`
- Time budget: ≈ 3 min in the talk, ≈ 70 s of clips (8 clips of 4–12 s), the rest spoken over held last frames
- Audience: PhD particle physicists
- **The one sentence**: the Z cross section brings the whole SM together in one number — the proton and QCD make it,
  the electroweak sector fixes how it couples, three lepton flavours must see the same value — and it is precise enough
  to separate PDF sets once the luminosity is under control.

## 2. Story (beats → clips)
| t | beat | clips |
|---|------|-------|
| 0:00–0:25 | one number, the whole SM | 1-02 `sm_to_process` |
| 0:25–1:05 | what the prediction is made of: factorisation, the α_s series, our band, the PDF-set spread | 1-03 `factorisation`, 1-04 `prediction` |
| 1:05–1:25 | a σ is only as good as its luminosity (Z counting; CMS 2024 is lumi-limited) | 1-05 `lumi_limited` |
| 1:25–1:45 | the Z under discoveries: H→ττ, Z→νν + jet | 1-06 `z_everywhere` |
| 1:45–2:25 | three flavours, one Z | 1-07 `lepton_universality` |
| 2:25–3:00 | what we measure: the 60–120 GeV window, the empty final plot | 1-08 `mass_window`, 1-09 `final_frame` |

Left out on purpose: the LEP lineshape / m_Z / N_ν, m_W and sin²θ_eff tensions (mass and couplings, not σ); σ vs √s and
the strange sea (possible later modules); our measured values (section 6).

## 3. Numbers on screen
All from `data/theory_reference.json` (`data/extract_theory_reference.py`, LCG env, 15 checks), asserted at scene import.

| value | clip | source |
|-------|------|--------|
| 1953.9 +56.1 −81.7 pb; scale +2.5/−3.9 %, α_s 1.3 %, PDF 0.7 % | 1-04, 1-08, 1-09 | `combination/combLieke/output/result.json` `prediction` |
| CT18 1921 +30 −33, MSHT20 1935 +23 −27, NNPDF3.1 1940 +15 −21, NNPDF4.0 1970 +11 −14 pb | 1-04 | arXiv:2408.03744 Table 5 (re-read 17 Sep 2026) |
| CMS 1952 ± 4 ± 18 ± 45 pb | 1-05 | `combination/combLieke/config/references.json` |
| x ∈ [6·10⁻⁴; 0.08]; LO 1970, NLO 1979, NNLO 1991, N³LO 2022 | 1-03 | computed; papers in the extractor |
| σ𝓑(H→ττ) = σ𝓑(Z→ττ)/560 | 1-06 | LHC HXSWG YR4 × BR(ττ) vs 1953.9 pb |
| ×207, ×17; 𝓑(Z→ℓℓ) 3.366 %; g_V^ℓ −0.037; Γ_μμ/Γ_ee 1.0009 ± 0.0028, Γ_ττ/Γ_ee 1.0019 ± 0.0032 | 1-07 | PDG 2024; Phys. Rept. 427 (2006) 257 |
| m_ℓℓ spectrum; 0.965 σ(m > 50) | 1-08 | `z-mumu/output/v2/gensums.json` `h_lhe_mll` |

## 4. Entry and exit
- Opens on: a white frame (the SM table); 1-07 opens on 1-01's final frame.
- Ends on: 1-09 = the final plot's σ axis (1650–2400 pb) with the aMC@NLO line and band and four empty rows
  (Z→ee green, Z→μμ gold, Z→ττ red | Z→ℓℓ slate). **Section 6 must open on this frame** (6-01 `three_to_one` is still
  the schematic skeleton; rebuild it on `F_*` / `fx()` of the scene).

## 5. Look
- Theory purple for predictions and σ̂, cyan for the PDFs ("what is tested"), slate for the CMS measurement, red for the
  one highlight (lumi); e / μ / τ in the channel colours where the three flavours appear together (1-07, 1-09) and for
  the μ / τ objects of 1-06.
- Words on screen limited to names on points (aMC@NLO, CMS 2024, PDF sets, LEP + SLD) and physics labels (stat, syst,
  lumi, scale, PDF, LO … N³LO, jet).
- Schematic: the diagram, the protons, the m_ττ shapes of 1-06 (Gaussians, only the area ratio is a number), the small
  event display.

## 6. Process
- Critique by the deck owner on the delivered v1; physics / style / continuity reviews only when asked.
- Git: leave for the owner (other sessions share the tree); push only when asked.
