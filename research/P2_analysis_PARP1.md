# BE3D Run — Decision Log & Benchmark: PARP1 (olaparib vs niraparib resistance, 3D)

> Runner: Queen "PARP1/PARP2" (round 2). Real public base-editing resistance screen (Coelho/Dincer
> 2024) → BE3D LFC3D + 3D clustering → per-inhibitor 3D resistance/essentiality atlas, benchmarked
> **honestly** against a curated, BE3D-independent ground truth (A0 discrimination test leading).

## 0. Target, data, structure
- **Gene / protein:** PARP1 (UniProt **P09874**, 1014 aa). Domains used (canonical numbering):
  ZnF1/2/3 225–359, BRCT 385–476, WGR 542–638, **HD (helical/autoinhibitory) 662–779**,
  **CAT (ART/catalytic) 788–1014** (NAD+/inhibitor pocket; H-Y-E triad H862 / Y896 / E988).
- **Structure:** AlphaFold **AF-P09874-F1-model_v6** (canonical WT, full length). Primary & only
  analysis structure (monomer mode; DSSP placeholder feeds characterization only).
- **Screen:** Coelho, Dincer et al. 2024 *Nat Genet* (DOI 10.1038/s41588-024-01948-8), Suppl. Table S4,
  sheet `ST2 BE z-scores`. **4,208 PARP1 guide rows.** The PARP-inhibitor arm was run in **MHH-ES-1**
  (Ewing sarcoma), **ABE only** (2,104 guides). Score columns present for PARP:
  `L2FC_MHHES1_Olap_plasmid_average_zscore` (**olaparib**) and `…_Nirap_…` (**niraparib**).
  **There is NO talazoparib or rucaparib arm in this dataset** — only olaparib & niraparib exist, so
  exactly two PARPi screens were run (not the four named in the brief). RESISTANCE = **positive** z
  (enrichment under drug).
- **Converter:** `scratchpad/convert_parp.py` (same schema as KRAS: category from
  `most_severe_consequence`; edits reconstructed from `Amino_Acid_Position`+`Edited_AA`+`New_AA`).
  PARP1 numbering matches canonical (**98.4 %** of missense guides' refAA match P09874 — clean, offset 0).
  Per drug: 1,981 scored rows → Missense 1,062 (693 distinct positions, span 1–1014), Silent 112,
  Splice-donor 134, Intron 39, No-Mutation 634.
- **QA (BE-QA):** cases=[Nonsense, Splice-donor] vs controls=[No Mutation, Silent].
  KS is **significant** (knockouts shift vs neutral) — a genuine QA *pass*, unlike KRAS. Interpretation
  still uses both directions (positive = resistance; negative = depletion/essentiality).

## 1. Parameter sweep (focused)
AF primary; `structure_radius` ∈ {4, 6, 8} + an r6 `sum` aggregation; `nRandom=1000`; clustering 6 Å;
p 0.05/0.01/0.001. **4 configs × 2 drugs = 8 runs, all EXITCODE 0.** Driver `scratchpad/run_grid_parp.sh`;
outputs `real_output/P2_PARP1/runs/<id>_<drug>/PARP1/`. Hits are stable across radius (counts grow
smoothly with radius; the domain-level conclusions below hold at every radius and for sum vs mean).

---

## 2. A0 DISCRIMINATION BENCHMARK (the primary test) — HONEST

Ground-truth sets (curated independently of BE3D, canonical numbering, cited in §5):
**CATALYTIC** (NAD+/inhibitor pocket + H-Y-E triad) = {862,863,864,871,878,888,889,896,903,904,907,988} (M=12);
**RESISTANCE** (trapping-loss: HD allostery + ZnF/WGR DNA-binding + Pettitt-2018 + Coelho-2024 resistance
residues) M=33; **FUNCTIONAL** = CATALYTIC ∪ RESISTANCE (M=44); **TOLERANT** (BRCT body + benign V762A) M=86.
N = 1014 structured residues; chance expectation E = K·M/N.

### 2a. Base rate — HIGH (the headline caveat)
At the **positive/resistance** direction, BE3D calls a large fraction of the protein significant:

| drug | #scored | pos p<0.05 | pos p<0.01 | pos p<0.001 |
|---|---|---|---|---|
| olaparib | 479 | 168 (35 % of scored, 17 % of protein) | 133 (28 %/13 %) | 100 (21 %/10 %) |
| niraparib | 628 | 284 (45 %/28 %) | 239 (38 %/24 %) | 200 (32 %/20 %) |

**With 10–28 % of the whole protein flagged as "resistance," any single overlap with a known site is weak
evidence.** This must frame everything below.

### 2b. Enrichment vs chance (Fisher / hypergeometric), r6 baseline
- **POSITIVE (resistance) vs RESISTANCE set — at / near chance.**
  - olaparib: p<0.01 → R=3, E=4.3, **R/E=0.69, OR=0.65, p=0.83 (NS)**; p<0.001 → R/E=0.92 (NS).
  - niraparib: p<0.01 → R=14, E=7.8, **R/E=1.8, OR=2.5, p=0.012**; p<0.001 → R/E=2.0, p=0.006.
  → **olaparib's resistance signal is indistinguishable from random; niraparib's is only weakly enriched.**
- **NEGATIVE (depletion) vs CATALYTIC pocket — STRONG, robust.**
  - olaparib: p<0.01 → R=12/12, E=2.0, **R/E=6.1, OR=∞, p=2.8×10⁻¹⁰**; p<0.001 → R/E=6.4, p=3.9×10⁻⁹.
  - niraparib: p<0.01 → **R/E=4.1, OR=8.8, p=4.3×10⁻⁴.**
  → **BE3D's depletion signal recovers the entire NAD+/inhibitor catalytic pocket** (essentiality).
- **UNION (either direction) vs FUNCTIONAL:** R/E ≈ 1.4–1.7, p = 0.005–0.035 — modest, and driven almost
  entirely by the negative/catalytic component.

### 2c. Discrimination gap (union hit-rate: FUNCTIONAL vs TOLERANT)
Functional residues are hit **2.6–5×** more often than tolerant residues (olaparib gap 4.1–5.0×,
niraparib 2.6–3.5×). Real but modest — inflated by the high base rate and by functional residues being
buried/core.

### 2d. Precision@10 (top residues by |LFC3D|)
- **Positive** top-10 (resistance rank): precision@10 = **0.2** vs RESISTANCE and vs FUNCTIONAL (both drugs).
  Top olaparib positives = {**1, 41, 45, 46**, 274, 275, 427, 428, 504, 834}; niraparib = {**1,13,19,28,39,41,45,46**,479,590}.
  The very top hits **do** include real resistance residues — **M1** (start-loss LOF → resistance, a
  Coelho hit) and the **ZnF1 DNA-binding cluster 43–46** (Pettitt) — but 8/10 are not in the curated set.
- **Negative** top-10 (depletion rank): {680–685 (HD helix), 827, 884, 890, 984 (CAT)} — the
  autoinhibitory HD helix + catalytic core, exactly the essential module.

### 2e. Beats a trivial baseline? — **MIXED; loses to AlphaMissense for resistance**
Two no-screen predictors on the SAME protein (top-tertile call, matched to BE3D hit counts):

| predictor | vs CATALYTIC | vs RESISTANCE | vs FUNCTIONAL |
|---|---|---|---|
| burial (top-tertile CA coordination) | R/E 1.97, p 0.01 | R/E 1.11, p 0.36 (NS) | R/E 1.37, p 0.035 |
| **AlphaMissense** (top-tertile mean pathogenicity) | R/E **2.94, p 2×10⁻⁶** | R/E **2.14, p 5×10⁻⁶** | R/E **2.34, p 3×10⁻¹⁰** |
| **BE3D negative (depletion)** | **R/E 6.1, p 3×10⁻¹⁰** ✅ | — | — |
| **BE3D positive (resistance)** | — | R/E 0.7–2.0 (olap NS; nirap p 0.01) ❌ | 1.4–1.7 |

**Verdict on baselines:** BE3D's **depletion** direction clearly **beats both burial and AlphaMissense at
finding the catalytic pocket** (R/E 6 vs 2.9 vs 2.0) — a real structure-function win. BE3D's
**resistance** direction **does not beat AlphaMissense**: a sequence-only predictor with *no screen at all*
recovers the curated PARPi-resistance residues better (R/E 2.1, p 5×10⁻⁶) than BE3D's positive signal
(olaparib at chance; niraparib R/E 1.8–2.0).

### 2f. HONEST VERDICT (PARP1)
**BE3D is genuinely intelligent in the depletion/essentiality direction and only marginal in the
resistance direction.** It recovers the PARP1 catalytic pocket and the autoinhibitory HD helix as a tight,
robust depletion hotspot (p ≈ 10⁻¹⁰, beats every baseline) — a strong internal positive control that the
LFC3D spatial aggregation works. But the **resistance** readout — the biologically interesting question —
has a **high base rate (10–28 % of the protein)** and **does not beat chance/AlphaMissense**: it *does*
surface the real resistance residues (M1, ZnF 43–46/119/120, WGR 589/591, BRCT 390/391) at the top of the
rank and at 4/4 sweep persistence, but with **poor specificity** (precision@10 ≈ 0.2). Trust the *ranked*
positive hits as hypotheses, not the significance labels as a discovery set.

---

## 3. Per-inhibitor 3D overlap (olaparib vs niraparib) — the "distinct maps?" question
Positive/resistance hit sets, baseline r6:

| threshold | #olaparib | #niraparib | shared | **Jaccard** |
|---|---|---|---|---|
| p<0.01 | 133 | 239 | 71 | **0.24** |
| p<0.001 | 100 | 200 | 54 | **0.22** |

Domain distribution of positive hits (p<0.01): both inhibitors concentrate in the **DNA-binding module —
ZnF (olap 24 / nirap 62) + WGR (8 / 36)** — i.e. loss of DNA engagement → loss of trapping → resistance,
with smaller HD (6/13) and CAT (13/10) components. **Niraparib's resistance map is roughly 2× broader than
olaparib's** and extends much further into the ZnF/WGR DNA-binding surface.

**Answer:** the two PARP inhibitors have **only moderately overlapping** 3D resistance maps (Jaccard ≈ 0.22–0.24,
comparable to KRAS sotorasib-vs-adagrasib ≈ 0.12), but this is driven mainly by **niraparib producing a
larger, more DNA-binding-centric hotspot** rather than by a clean *orthosteric-vs-distal* split. It is **not**
the crisp drug-specific escape geometry KRAS showed; the high base rate blurs the contrast. Honest read:
**a difference in breadth/DNA-binding involvement, not a cleanly separable per-drug pocket signature.**

**HD allosteric domain:** yes — a minority of *resistance* (positive) hits fall in the HD autoinhibitory
domain (olap 6, nirap 13 at p<0.01; 3 persist ≥3/4 runs), consistent with HD mutations that release
autoinhibition/reduce trapping. But the **dominant** HD signal is in the **depletion** direction (32 HD
residues neg-significant), i.e. HD perturbation is mostly deleterious here.

---

## 4. Ground-truth classification of robust hotspots (TP / novel / FP)
Robust set (p<0.01 baseline + persistence). Full table: `outreach/PARP1/PARP1_hotspots.tsv` (+ `_full`).
- **TP — catalytic pocket (depletion):** all 12 CAT residues **862,863,864,871,878,888,889,896,903,904,907,988**
  (H-Y-E triad + NAD pocket), strong negative LFC3D, 3–4/4 persistence. *De novo* rediscovery of the active site.
- **TP — HD autoinhibitory helix (depletion):** the contiguous 678–701 / 767–779 HD stretch (deep negative,
  4/4). Matches the Langelier/Pascal autoinhibition module.
- **TP — resistance (positive, validated):** **M1** (start-loss LOF), **F44/D45, K119/S120, G161/V164**
  (ZnF DNA-binding; Pettitt 2018), **W246**, **W589/R591** (WGR), **L390/S391** (BRCT), **Y848** (CAT) —
  14 residues, mostly 4/4 persistent.
- **NOVEL-CANDIDATES:** ~229 positive residues persistent across the sweep with no prior annotation —
  treat as hypotheses only (the high base rate means many are passengers).
- **LIKELY FALSE POSITIVES:** ~70 residues with pLDDT < 70 (disordered N-terminal/linker regions) — flagged.
- **FALSE NEGATIVES:** because the PARPi arm is **ABE-only**, C→T (CBE) escape alleles are unreachable; the
  catalytic pocket appears only as depletion (not resistance) because catalytic-dead PARP1 is still trapped
  (trapping toxicity), so pocket point-mutants deplete rather than enrich — a direction artifact, not a miss.

---

## 5. Ground-truth citations
- **Catalytic H-Y-E triad & NAD+/inhibitor pocket** (H862, Y896, E988, G863, S904, Y907, K903, G888, Y889,
  D766): Langelier et al. 2012 *Science* 336:728 (PARP1 CAT structure); Ruf et al. 1996/1998 (chicken PARP
  CAT, catalytic Glu); olaparib co-crystal 5DS3 (nicotinamide pocket Y907/S904/G863). 
- **HD autoinhibitory / allosteric domain & PARP-trapping** (662–779; L698/L701; reverse/allostery):
  Langelier et al. 2018 *Nat Commun* 9:844; Dawicki-McKenna et al. 2015 *Mol Cell*; Zandarashvili et al.
  2020 *Science* 368 (allosteric retrapping).
- **PARPi-resistance residues** (ZnF 43–45/119–120, WGR 591, HD 688/742/743, CAT 848/925): Pettitt et al.
  2018 *Nat Commun* 9:1849 (genome-wide PARP1 mutational resistance screen). **Coelho/Dincer 2024** report
  PARP1 start-loss (M1) and additional base-edit resistance alleles.
- **DNA-binding ZnF residues**: Ali et al. 2012; Eustermann et al. 2015 *Mol Cell*.
(Literature curation performed by a dedicated sub-agent via PubMed/paperclip; residue identities verified
against the P09874 sequence.)

## 6. Competitor comparison
| tool | input signal | on PARP1 | advantage | disadvantage |
|---|---|---|---|---|
| **BE3D** (this run) | BE resistance screen LFC → LFC3D | catalytic+HD depletion cluster (p≈10⁻¹⁰); broad, low-specificity resistance map | sees screen-based essentiality + editor-reachable resistance; drug-specific arms | resistance dir. high base rate, doesn't beat AlphaMissense; ABE-only coverage |
| **AlphaMissense** | seq/structure predictor (no screen) | R/E 2.1–2.9 for catalytic & resistance sets | zero-cost, no screen; beats BE3D on resistance | can't see drug/editor-specific escape or trapping mechanism |
| **3dhotspots.org / cBioPortal** | somatic mutation recurrence | PARP1 not a recurrent cancer hotspot gene → little signal | decade of benchmarking | blind to drug-resistance & functional-screen signal |
| **ProTiler-Mut** (3D-RRA) | tiling screen → 3D | conceptually similar; would also flag the catalytic cluster | direct 3D tiling rival | recurrence/RRA null differs from BE3D's randomization |

## 7. Deliverables
- `outreach/PARP1/PARP1_hotspots.tsv` (+ `_full`) — flagged robust hotspots.
- `outreach/PARP1/PARP1_Olap_G2P.tsv`, `PARP1_Nirap_G2P.tsv` — interactive 3D view input.
- `outreach/PARP1/PARP1_brief.md` + `PARP1_email.md`.
- Machine-readable: `real_output/P2_PARP1/analysis/{discrim,cross,conservation}.json`, `am_baseline.json`.
