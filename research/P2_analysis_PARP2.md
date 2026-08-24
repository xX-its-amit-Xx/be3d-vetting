# BE3D Run & Benchmark: PARP2 (olaparib/niraparib resistance, 3D) + PARP1↔PARP2 conservation

> Runner: Queen "PARP1/PARP2" (round 2). Companion to `P2_analysis_PARP1.md`. Same Coelho/Dincer 2024
> base-editing resistance screen; same pipeline; **A0 discrimination test leading and honest.**

## 0. Target, data, structure
- **Gene / protein:** PARP2 (UniProt **Q9UGN5**, 583 aa canonical). Domains (canonical numbering):
  WGR 104–201, **HD 231–348**, **CAT 356–583** (NAD+/inhibitor pocket; H-Y-E triad **H428 / Y462 / E558**).
- **Structure:** AlphaFold **AF-Q9UGN5-F1-model_v6** (canonical WT, 583 aa). Monomer mode.
- **Screen:** Coelho/Dincer 2024 `ST2 BE z-scores`, **1,618 PARP2 guide rows**; PARPi arm in **MHH-ES-1,
  ABE only (809 guides)**. Drugs present: **olaparib** (`L2FC_MHHES1_Olap_…zscore`), **niraparib**
  (`…Nirap…`). No talazoparib/rucaparib. RESISTANCE = positive z.
- **⚠ Numbering offset (important):** the screen's `Amino_Acid_Position` uses a **PARP2 transcript numbered
  canonical − 13** (isoform-2 / short-isoform convention). Raw match to canonical Q9UGN5 was only 14 %; after
  applying **+13 (screen → canonical/AF)** the refAA match rose to **86 %** — clean. All results below are in
  canonical Q9UGN5 numbering (screen pos + 13). Converter `scratchpad/convert_parp.py` (`offset=13`).
  Per drug: 798 scored rows → Missense 502 (334 distinct positions, span 14–580), Silent 68,
  Splice-donor 81, Intron 15, No-Mutation 132.
- **QA (BE-QA):** cases vs controls KS is **significant (D=0.28, p=1.3×10⁻⁴)** — a genuine QA *pass*
  (knockouts shift vs neutral).

## 1. Sweep
AF primary; `structure_radius` ∈ {4,6,8} + r6 `sum`; `nRandom=1000`; clustering 6 Å. **4 configs × 2 drugs
= 8 runs, all EXITCODE 0.** `real_output/P2_PARP2/runs/…`.

---

## 2. A0 DISCRIMINATION BENCHMARK — HONEST (weaker than PARP1)

Ground truth (canonical): **CATALYTIC** (NAD+/inhibitor pocket + triad) =
{428,429,430,455,462,469,470,473,558} (M=9). **No residue-specific PARPi-resistance set or curated
tolerant set is available for PARP2** (Pettitt-style saturation resistance was done on PARP1), so the
resistance-direction and discrimination-gap tests cannot be run for PARP2 — reported as N/A, not as a pass.
N = 583 structured residues; **all 9 catalytic residues are scored** (in an edited-residue neighborhood).

### 2a. Base rate — HIGH
| drug | #scored | pos p<0.05 | pos p<0.01 | pos p<0.001 |
|---|---|---|---|---|
| olaparib | 201 | 83 (41 % of scored) | 66 (33 %) | 44 (22 %) |
| niraparib | 207 | 58 (28 %) | 49 (24 %) | 33 (16 %) |

### 2b. Enrichment vs chance (r6 baseline)
- **NEGATIVE (depletion) vs CATALYTIC — NOT significant.** olaparib R/E = 1.4–1.9 (p 0.23–0.45);
  niraparib R/E = 0.9–1.7 (p 0.26–0.70). Only 2–3 of the 9 catalytic residues reach neg-significance.
- **UNION vs FUNCTIONAL (=CATALYTIC here):** R/E ≈ 1.0 (p 0.5–0.8) — **at chance.**
- **Precision@10:** 0.0 for both catalytic and functional, both directions.

### 2c. Beats baseline?
AlphaMissense (top-tertile) vs CATALYTIC: **R/E 2.94, p 5×10⁻⁵** — i.e. a no-screen predictor *does*
concentrate on the pocket, while **BE3D does not**. **BE3D loses to AlphaMissense on PARP2.**

### 2d. HONEST VERDICT (PARP2)
**BE3D's discrimination essentially fails on PARP2 at the residue level.** Neither direction significantly
concentrates on the catalytic pocket, precision@10 is zero, and a sequence-only baseline beats it. The most
likely cause is **coverage/power**: the PARPi arm is **ABE-only with just 809 guides → ~200 scored residues**
over a 583-aa protein, so the catalytic domain is thinly and non-uniformly tiled. The signal that exists is
diffuse (top hits scatter across WGR/CAT/other; e.g. olaparib positive top-10 = 36,41,86,87,93,194,313,481,
581,582). **Do not present PARP2 hotspots as validated;** this run is under-powered.

---

## 3. Per-inhibitor overlap (olaparib vs niraparib)
Positive hits, r6: p<0.01 → #olap 66, #nirap 49, shared 30, **Jaccard 0.35**; p<0.001 → **Jaccard 0.33**.
Higher overlap than PARP1 (0.22) — for PARP2 the two inhibitors' (weak) maps are **more similar**, both
concentrating in **WGR + CAT**. Given the failed discrimination, this concordance is of limited significance
(likely reflects shared coverage rather than shared biology).

---

## 4. PARP1 ↔ PARP2 cross-paralog conservation (the paralog question)
Are the two paralogs' resistance/functional surfaces spatially conserved? Structural equivalents of the
PARP1 catalytic pocket (P1→P2): 862→428, 863→429, 864→430, 889→455, 896→462, 903→469, 904→470, 907→473,
988→558 (agent-curated from the conserved ART fold).

| observation | PARP1 | PARP2 |
|---|---|---|
| catalytic pocket direction | **depletion (negative)** — all 9 equiv. residues neg-significant | depletion, but only **3/9** (430,455,462) reach significance |
| resistance (positive) at the pocket | **none** of the pocket is positive | **none** positive |

**Answer:** the catalytic pocket is **structurally conserved** and, in **both paralogs, falls in the
depletion (essential) direction — never the resistance direction** — a consistent cross-paralog signature
(mutating the NAD+ pocket is deleterious, not protective, in both). BE3D **recovers this conserved pocket
strongly in PARP1 (p≈10⁻¹⁰) but only weakly in PARP2 (3/9, NS)** — the conservation is real in the biology
and structure; BE3D's *power* to detect it is paralog-dependent and limited by PARP2's sparse ABE-only
coverage. So: spatially conserved surface, **yes**; equally recoverable by BE3D, **no**.

## 5. Ground-truth citations
- PARP2 catalytic domain & H-Y-E triad (H428/Y462/E558), NAD+/inhibitor pocket: Oliver et al. 2004
  *Nucleic Acids Res* (PARP2 CAT); Obaji et al. 2018/2021 *Nat Commun* (PARP2 DNA-dependent activation);
  PARP2–inhibitor co-structures. Cross-paralog ART fold conservation: Langelier et al. 2018.
- No PARP2-specific saturation resistance map exists (Pettitt 2018 profiled PARP1); hence no PARP2
  RESISTANCE/TOLERANT ground-truth sets.

## 6. Competitor comparison
Same landscape as PARP1 (§6 there). For PARP2 specifically, **AlphaMissense is the better tool** here
(recovers the pocket, R/E 2.9, p 5×10⁻⁵, with no screen), because BE3D's screen-based signal is
power-limited. BE3D's niche (drug/editor-specific escape, trapping allostery) needs a deeper PARP2 screen
(add a CBE arm; more guides) to be realized.

## 7. Deliverables
- `outreach/PARP2/PARP2_hotspots.tsv` (+ `_full`) — flagged, **with the explicit under-power caveat**.
- `outreach/PARP2/PARP2_Olap_G2P.tsv`, `PARP2_Nirap_G2P.tsv`.
- `outreach/PARP2/PARP2_brief.md` + `PARP2_email.md`.
- `real_output/P2_PARP2/analysis/{discrim,cross}.json`, `real_output/P2_PARP1/analysis/conservation.json`.
