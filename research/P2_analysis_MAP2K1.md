# BE3D Run — Analysis: MAP2K1 / MEK1 (trametinib resistance, 3D)

> Runner: Queen "MEK1/MEK2". Real public base-editing screen (Coelho/Dincer 2024) → BE3D LFC3D +
> 3D clustering → a 3D map of **trametinib** (allosteric/type-III MEKi) resistance on MEK1.
> NEW question: does MEKi resistance form a **tight 3D shell around the allosteric pocket**
> (distinct from the ATP site), and are any resistance residues **distal** (allosteric-of-allosteric)?

## 0. Target & data
- **Gene / protein:** MAP2K1 / MEK1, UniProt **Q02750**, 393 aa. **Primary structure: AlphaFold `AF-Q02750-F1-model_v6`** (WT, canonical numbering == screen).
- **Pocket-defining structure:** PDB **7JUR** chain C — MEK1 with **ANP** (AMP-PNP, ATP-site) + **QOM** (a trametinib-class allosteric inhibitor) + Mg. 7JUR numbering == AF (303/304 CA match). This crystal lets us split the **allosteric pocket** (<6 Å of QOM) from the **ATP site** (<6 Å of ANP), independently of BE3D.
- **Screen:** Coelho, Dincer et al. 2024 *Nat Genet* (DOI 10.1038/s41588-024-01948-8), sheet `ST2 BE z-scores`. The only MEK-inhibitor arm present is **trametinib** in HT29 (BRAF-V600E CRC): `L2FC_HT29_Tram_plasmid_average_zscore`. RESISTANCE = **positive** z (enrichment under drug). Pooled ABE+CBE.
- **Converter:** `real_output/P2_MAP2K1/input/convert_mek.py`. 2,271 guide rows → Missense 722 (337 distinct positions, span 1–393), Silent 173, Nonsense 31, Splice-donor 110, Intron 31, No Mutation 1,204. refAA validated against the AF sequence (0.9% mismatch dropped).

## 1. QA (BE-QA)
- cases = [Nonsense, Splice-donor], controls = [No Mutation, Silent]. **KS D=0.211, p=1.9e-5 (141 cases vs 1,377 controls) — SIGNIFICANT.** Unlike the KRAS resistance screens (null QA), MEK knockouts *do* separate from neutrals in the trametinib arm (MEK is essential downstream of BRAF-V600E in HT29). **Screen ACCEPTED.** Downstream interpretation uses the **positive/resistance** direction.

## 2. Sweep (AF, nRandom 1000; structure_radius 4/6/8 mean + r6 sum; clustering 6 Å) — 4 runs, all EXITCODE 0
`real_output/P2_MAP2K1/runs/AF_r{4,6,8}_mean_n1000/`, `AF_r6_sum_n1000/`.

| config | base%@.05 | base%@.001 | FUNC OR@.01 | LIT OR@.01 | gap@.01 | prec@10(LIT) |
|---|---|---|---|---|---|---|
| r4 mean | 20.6 | 13.0 | 5.44 | 8.33 | 6.9× | 0.4 |
| **r6 mean (baseline)** | **20.1** | **13.2** | **2.76** | **7.49** | **5.1×** | **0.4** |
| r8 mean | 27.0 | 20.1 | 5.12 | 14.0 | 4.3× | 0.3 |
| r6 sum | 20.9 | 15.0 | 5.32 | 15.57 | ∞ | 0.5 |

Signal is **stable across the sweep**; larger radius inflates the base rate (r8 → 20% at p<0.001) without improving discrimination. Baseline = r6 mean.

## 3. A0 DISCRIMINATION BENCHMARK (the primary test) — baseline r6, positive direction
Independent residue sets (curated WITHOUT BE3D, citations in §6):
- **ALLO6** (structural allosteric pocket, <6 Å QOM, M=27): 97,99,115,118,126–129,141,143,188–190,195,206–213,215,216,219,234,240.
- **LIT** (literature-validated trametinib/allosteric-MEKi resistance, M=19): 53,56,57,98,99,103,115,118,121,124,128,129,141,143,203,207,209,211,215.
- **FUNC** = ALLO6 ∪ LIT (the escape surface, M=35).
- **TOLERANT** (solvent-exposed, non-pocket, non-ATP, AM-benign, ordered, M=15): 28,273,277,280,289,292,302,304,328,354,357,363,366,382,385.

| metric | p<0.05 | p<0.01 | p<0.001 |
|---|---|---|---|
| #sig / N=393 (**base rate**) | 79 (**20.1%**) | 69 (17.6%) | 52 (**13.2%**) |
| Enrichment vs **LIT** (R; R/E; OR; p) | 11; 2.9×; **6.2**; **2.2e-4** | 11; 3.3×; **7.5**; 5.5e-5 | 10; 4.0×; **8.8**; **2.6e-5** |
| Enrichment vs FUNC | 13; 1.85×; 2.6; 0.011 | 12; 1.95×; 2.8; 0.010 | 10; 2.2×; 3.0; 0.009 |
| Enrichment vs pure ALLO6 pocket | 8; 1.47×; 1.75; **0.15 (ns)** | 7; 1.5×; 1.7; 0.18 ns | 6; 1.7×; 2.0; 0.13 ns |
| **Discrimination gap** (hit FUNC% / TOL%) | 37.1/6.7 = **5.6×** | 34.3/6.7 = 5.1× | 28.6/**0.0** = ∞ |
| **Precision@10 / @20** (LIT; base 0.09) | — | — | **0.40 / 0.25** |
| Baseline "call all buried" (36% of protein) vs FUNC | OR 3.4 | | |
| Baseline "AlphaMissense≥0.9" (45% of protein) vs FUNC | **OR 24.5** (recovers 33/35 FUNC) | | |

**HONEST VERDICT (MEK1): BE3D beats chance — clearly on the literature-validated resistance set** (OR 6–9, p<1e-4; precision@10 = 4.4× base; functional hit-rate ≫ tolerant, gap 5×→∞). But two honest caveats: (1) the **base rate is high (13–20%)**, so overlap at p<0.05 is weak evidence — credibility rests on **p<0.001**. (2) BE3D does **not** enrich for the *pure structural* allosteric pocket beyond chance (ALLO6 p≈0.15); its signal is carried by the **validated-resistance positions** (which include the N-lobe activating cluster). **BE3D does NOT beat the drug-agnostic AlphaMissense baseline on raw enrichment** — AM (≥0.9) recovers 33/35 functional residues wholesale (it flags the entire buried conserved kinase core). BE3D's ADDED value over AM is (a) a **resistance DIRECTION** AM cannot provide, (b) a **3.4× tighter call set** (52 vs 177 residues), and (c) recovery of the drug-relevant subset. **Net: genuine but modest discrimination; strongest as a resistance-direction filter, not a standalone pathogenicity predictor.**

## 4. Robust hotspots (pos p<0.01 in ≥3/4 runs) — 54 residues; see `outreach/MAP2K1/MAP2K1_hotspots.tsv`
- **TP — pocket-disrupting trametinib resistance (9):** C121, L115, L118, F129, M143 (gatekeeper), Q116-region, I99, I216, M219 — the αC/β3 + catalytic-loop + activation-segment **allosteric wall**. C121/L115/L118/M143 are AM-LPath (0.99+) and directly validated (Emery 2009, Wagle 2011, Zhang 2025).
- **TP — activating / 3dhotspots N-lobe cluster (17):** the **helix-A negative-regulatory region 44–61** (incl. F53, Q56, K57) + G128, Y130, E203. These are the recurrence-based 3dhotspots set (Chang 2017).
- **Novel pocket-shoulder (14):** αC/β3 residues 112–123 (6–12 Å from QOM) — the pocket rim; structurally plausible resistance shoulder.
- **Novel DISTAL, structured (the answer to the NEW question):** **357, 371, 372** — a **C-lobe αH/αI patch, pLDDT 92–98, 18–28 Å from the allosteric pocket.** Candidate *allosteric-of-the-allosteric*. **Caveat:** this stretch (F371/A372/G373 region) is a documented **non-trametinib** MEKi (binimetinib/selumetinib/cobimetinib) resistance region (Zhang 2025) — so under trametinib it may reflect general MEKi cross-resistance rather than a trametinib-specific escape. Treat as hypothesis. Also 199–202 = activation-loop N-term, adjacent to validated E203.
- **LIKELY FALSE POSITIVES (4):** **284–287** (αF–αG loop, **pLDDT 40–46**, disordered) — low-confidence geometry, flag as spurious. (Also low-pLDDT C-terminal 387–389.)

**FALSE NEGATIVES:** P124 (validated resistance) — the required codon change is base-editing-inaccessible / no guide; V211, F209, I141, I103 (pocket) are in ALLO6 but not called (single-guide coverage / signal below threshold). BE3D's misses are dominated by **editor reachability**, not geometry.

## 5. NEW question — answered
**YES.** Trametinib resistance forms a **tight 3D shell on the allosteric-pocket wall** (β3–αC 112–131, catalytic loop 142–144, activation segment 199–220), **distinct from the ATP catalytic core** — the ATP-only residues (74, 81, 82, 95, 144–153) are essentially *not* called; the called pocket residues (99,115,118,128,129) are **ALLO-only, not ATP**. On top of that, BE3D recovers the **N-lobe activating cluster (44–61)** that 3dhotspots finds. And it nominates a **structured distal C-lobe patch (357/371/372)** as a candidate allosteric-of-the-allosteric site (caveated above). So the map cleanly separates: allosteric shell (drug-specific) vs ATP site (not resistance) vs N-lobe activator vs a distal candidate.

## 6. Competitor comparison
| tool | input signal | on MEK1 | advantage | disadvantage |
|---|---|---|---|---|
| **BE3D** (this run) | trametinib screen LFC, 3D-aggregated + null | allosteric shell + N-lobe activator + distal candidate; **directional (resistance)** | drug-specific escape geometry; controls; sees allosteric+distal sites | high base rate; needs p<0.001; doesn't beat AM on raw enrichment |
| **3dhotspots.org** (Chang 2017, PMC5260099) | patient mutation **recurrence** | R49/A52/F53/Q56/K57/G128/Y130 (N-lobe activators only) | decade of benchmarking; no screen needed | **misses the allosteric drug-resistance pocket** (C121,L115,I141,V211,F209 rare in cohorts); flags trametinib-*sensitive* variants (Q56P, Y130C) |
| **AlphaMissense** | seq/structure pathogenicity, **drug-agnostic** | flags the buried conserved pocket (M143,I141,L115,C121,V211,F209) for the *wrong* reason | no screen needed; strong raw enrichment (OR 24) | non-directional; flags 45% of protein; no resistance concept; mishandles gain-of-function activators |

## 7. Confidence caveats
- Missense z assigned per-guide to all edited residues (within-guide causality unresolved); robustness over 337 positions + 4 configs mitigates. Direction: positive/resistance only. Distal 357/371/372 rests on AF (high pLDDT) but overlaps a literature *non-trametinib* MEKi region — hypothesis, not claim. 284–287 low-pLDDT → flagged FP. Literature grounding: Khan 2020 *Nature* PMID 32927473 (7JUR/trametinib), Emery 2009 *PNAS* PMID 19915144, Wagle 2011 *JCO* PMID 21383288, Zhang 2025 (PMC12022096), Chang 2017 *Genome Med* PMC5260099, Coelho 2024 *Nat Genet*.
