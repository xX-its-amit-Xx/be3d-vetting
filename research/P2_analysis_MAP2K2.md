# BE3D Run — Analysis: MAP2K2 / MEK2 (trametinib resistance, 3D) + MEK1↔MEK2 spatial conservation

> Runner: Queen "MEK1/MEK2". Same drug class (trametinib), paralog of MEK1.
> NEW question: are MEK2 resistance hotspots **spatially conserved** with MEK1's — i.e. is the
> drug-escape surface structurally conserved across paralogs? (a cross-paralog comparison).

## 0. Target & data
- **Gene / protein:** MAP2K2 / MEK2, UniProt **P36507**, 400 aa. **Primary structure: AlphaFold `AF-P36507-F1-model_v6`** (WT, canonical numbering == screen).
- **Screen:** Coelho 2024 `ST2 BE z-scores`, **trametinib** arm in HT29 (`L2FC_HT29_Tram_plasmid_average_zscore`); RESISTANCE = positive z; pooled ABE+CBE. Converter `real_output/P2_MAP2K2/input/convert_mek.py`: 2,146 rows → Missense 837 (369 distinct positions, span 1–400), Silent 186, Nonsense 37, Splice-donor 170, Intron 54, No Mutation 862. refAA-validated (0.9% dropped).
- **Pocket definition:** MEK2 has no trametinib crystal here; the pocket + ATP site were transferred from MEK1 7JUR by **structural superposition** (§4), giving a MEK2 `pocket_dist.json` in the common frame.

## 1. QA
- cases=[Nonsense,Splice-donor], controls=[No Mutation,Silent]: **KS D=0.190, p=6.2e-6 (207 vs 1,048) — SIGNIFICANT.** Screen ACCEPTED; positive/resistance direction used.

## 2. Sweep (AF, n1000; r4/6/8 mean + r6 sum) — 4 runs, all EXITCODE 0
| config | base%@.05 | base%@.001 | FUNC OR@.01 | LIT OR@.01 | gap@.01 | prec@10(LIT) |
|---|---|---|---|---|---|---|
| r4 mean | 28.8 | 23.2 | 3.48 | 2.56 | 4.5× | 0.3 |
| **r6 mean (baseline)** | **32.0** | **24.5** | **3.13** | **3.66** | **8.9×** | **0.2** |
| r8 mean | 30.0 | 23.8 | 1.91 | 3.86 | 6.9× | 0.2 |
| r6 sum | 27.8 | 20.5 | 2.14 | 3.81 | 3.2× | 0.3 |

## 3. A0 DISCRIMINATION BENCHMARK — baseline r6, positive direction
Sets (independent; MEK2 pocket from superposition, LIT = MEK1 validated set +4 offset):
- **ALLO6** (M=34): 82–84,101,103,119,122,130–133,145,147,192–194,199,210–217,219–221,223,224,227,234,238,244.
- **LIT** (M=19, +4): 57,60,61,102,103,107,119,122,125,128,132,133,145,147,207,211,213,215,219.
- **FUNC** = ALLO6∪LIT (M=42). **TOLERANT** (M=17): 35,36,74,76,273,277,281,335,336,339,340,365,371,374,376,388,390.

| metric | p<0.05 | p<0.01 | p<0.001 |
|---|---|---|---|
| #sig / N=400 (**base rate**) | 128 (**32.0%**) | 115 (28.8%) | 98 (**24.5%**) |
| Enrichment vs **FUNC** (R;R/E;OR;p) | 23; 1.7×; 2.9; 1.1e-3 | 22; 1.8×; **3.1**; 5.6e-4 | 21; 2.0×; **3.65**; **1.3e-4** |
| Enrichment vs LIT | 11; 1.8×; 3.1; 0.015 | 11; 2.0×; 3.7; 0.006 | 10; 2.2×; 3.7; 0.006 |
| Enrichment vs pure ALLO6 pocket | 18; 1.65×; 2.6; 0.007 | 17; 1.7×; 2.7; 0.005 | 16; 1.9×; 3.1; 0.002 |
| **Discrimination gap** (FUNC%/TOL%) | 54.8/11.8 = 4.6× | 52.4/5.9 = **8.9×** | 50.0/5.9 = 8.5× |
| **Precision@10 / @20** (LIT; base 0.10) | — | — | 0.20 / 0.20 (FUNC@20 0.30) |
| Baseline "buried" (35%) vs FUNC | OR 2.85 | | |
| Baseline "AlphaMissense≥0.9" (38%) vs FUNC | **OR 44.5** (recovers 40/42) | | |

**HONEST VERDICT (MEK2): signal is REAL but NOISIER than MEK1.** Enrichment is significant on all three sets — including the **pure structural allosteric pocket** (OR 2.6–3.1, p≤0.007), which MEK1 did *not* reach — and functional ≫ tolerant (gap 5–9×). **But the base rate is very high (24–32%): BE3D flags a quarter-to-a-third of the protein**, so any single overlap is weak evidence, and **precision@10 is only 2× base** (vs 4.4× for MEK1). The N-terminal 1–32 disordered tail is heavily over-called (25 robust low-pLDDT hits → likely false positives). **Does NOT beat AlphaMissense** on raw enrichment (AM OR 44.5, recovers 40/42). **Net: positive but use p<0.001 + a pLDDT>60 filter; treat as hypothesis-generating, not a precise hotspot list.**

## 4. MEK1 ↔ MEK2 SPATIAL-CONSERVATION COMPARISON (the NEW question)
Method: pairwise-aligned the two AF sequences (identity 318/392 = 81% over the aligned core), superposed MEK2-AF onto MEK1-AF (**RMSD 2.14 Å over 392 CA**), and transferred the 7JUR QOM/ANP ligands into the common frame (7JUR→MEK1 RMSD 1.17 Å). Data: `targets/mek_crossparalog.json`.

**(a) The pocket itself is essentially identical.** All **27/27** MEK1 allosteric-pocket residues map to MEK2 with a **clean +4 offset**, every one is the **SAME amino acid**, and CA–CA distances after superposition are **≤0.6 Å**:
`K97→K101, I99→I103, L115→L119, L118→L122, I126→I130, V127→V131, G128→G132, F129→F133, I141→I145, M143→M147, D208→D212, F209→F213, V211→V215, S212→S216, L215→L219, M219→M223, R234→R238 …` (C121→**C125**, the paralog cysteine of the C121S resistance mutation). The trametinib escape *architecture* is structurally conserved.

**(b) The RESISTANCE HOTSPOTS are spatially conserved.** Mapping MEK1's robust hits (pos p<0.01) into MEK2 numbering: **44/69 MEK1 hits (64%) fall on the SAME aligned position that is ALSO a MEK2 hit** (Jaccard 0.31 in the aligned frame). For the pocket specifically, **6/7 MEK1 pocket hits are shared**: MEK1 {99,115,118,128,129,143,219} ↔ MEK2 {119,122,132,133,147,223}, all called in both screens. Both paralogs also light up the **N-lobe helix-A activating cluster** (MEK1 44–61 ↔ MEK2 48–65) and the catalytic-loop gatekeeper (MEK1 M143 ↔ MEK2 M147).

**VERDICT: YES — the MEK1 and MEK2 trametinib-escape surfaces are spatially conserved.** The same allosteric wall (β3–αC + catalytic loop + activation segment) and the same N-lobe activator light up in both paralogs, at structurally superposable positions with identical residues. This is a genuine cross-paralog structural conservation of the drug-escape surface, recovered independently from two separate screens.

## 5. Robust hotspots — 95 residues; `outreach/MAP2K2/MAP2K2_hotspots.tsv`
- **TP pocket-resistance / structural (13):** C125 (=MEK1 C121), L119, L122, F133, M147 (gatekeeper), V215-region, E207 — the conserved allosteric wall.
- **TP activating / 3dhotspots (17):** helix-A 49–65 (incl. F57, Q60, K61), Y134, E207.
- **Novel pocket-shoulder (16)** and **novel distal (24 structured):** incl. C-lobe 268,269,278,317–319,374,375 (paralog-consistent with MEK1's 357/371/372 distal patch).
- **LIKELY FALSE POSITIVES (25):** N-terminal **1–32** (disordered, low pLDDT) + 225,226,293,294,302. The N-terminal over-calling is the dominant failure mode here — explicitly flagged.
- **FALSE NEGATIVES:** base-editing-inaccessible pocket codons; some ALLO6 residues below threshold.

## 6. Competitor note
Same as MEK1: **3dhotspots** would report only the MEK2 N-lobe activating codons (recurrence-based; misses the allosteric pocket and the paralog-conservation story); **AlphaMissense** flags the conserved pocket wall drug-agnostically (OR 44.5) but gives no resistance direction and over-calls. BE3D's unique contribution is the **directional, drug-specific, cross-paralog-conserved escape surface**.

## 7. Caveats
High base rate (24–32%); N-terminal FP over-calling; MEK2 pocket is transferred by superposition (RMSD 2.14 Å) not a native trametinib crystal; per-guide z assigned to all edited residues. Grounding refs as in P2_analysis_MAP2K1.md (Khan 2020; Emery 2009; Zhang 2025; Chang 2017; Coelho 2024).
