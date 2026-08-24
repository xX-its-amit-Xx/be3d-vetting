# BRCA1 base-editing screen → 3D functional hotspots (BE3D) — 1-page brief

**For:** clinical-genetics / VUS interpretation labs. **Tool:** BE3D (Iqbal Lab, Broad; Liau Lab, Harvard) — aggregates base-editor tiling-screen log-fold-change over 3D structural neighborhoods (LFC3D) with a permutation null, then spatially clusters significant residues into hotspots. **This is an independent analysis, not a BE3D-author endorsement.**

## Input
- **Screen:** Cuella-Martin et al. 2021, *Cell* (PMID 33592168), DDR base-editing screens (Table S5). BRCA1, CBE (BE3-FNLS), competitive **dropout ± olaparib and ± cisplatin**, MCF10A. 519 guides → two BE3D screens (olaparib, cisplatin) meta-aggregated. **Assay direction:** dropout — loss-of-function = negative LFC (verified: nonsense mean −0.33 vs synonymous +0.09; drug arms roughly double the LOF dropout via synthetic lethality). QA passed (nonsense vs no-mutation KS p = 3e-3 / 9e-4).
- **Structure:** AlphaFold P38398 (full length, primary); PDB 1JM7 (BRCA1–BARD1 RING) and 1T29 (BRCT) for domain detail.
- **Gold standards (independent of BE3D):** Findlay et al. 2018 SGE (PMID 30209399; per-residue function scores, RING 1–101 + BRCT 1631–1855) and ClinVar; AlphaMissense as the predictor competitor.

## Headline discrimination results (whole-protein, primary run r6/mean/n1000)
| Metric | Value | Read |
|---|---|---|
| Genome-wide significant-residue base rate | 7.4% (p<.05), 4.9% (p<.001) | Not "everything is a hotspot" |
| Domain hit-rate: RING / BRCT / central-IDR | 12.8% / 20.1% / 5.3% | Concentrates on folded functional domains (~2.4–3.8×) |
| Whole-protein enrichment for validated functional residues | OR 4.3, p 6×10⁻⁸ | Real, but domain-driven |
| **Trivial pLDDT>70 ("is it folded") baseline** | **enr 5.1 (≥ BE3D's 3.1)** | **BE3D does *not* beat this at domain level** |
| Findlay-SGE residue-level enrichment (within RING+BRCT) | enr 1.4, OR 1.8, p≈0.06 | Weak fine-grained resolution |
| Per-residue concordance with Findlay function score (Spearman) | ≈ 0.05 | No quantitative rank agreement |

**Verdict:** BE3D discriminates at the **domain** level (functional domains ≫ disordered middle) but that is largely recapitulated by a folded-vs-disordered / AlphaMissense baseline; it adds only **modest residue-level** signal. It is complementary evidence, not a standalone classifier.

## Hotspot calls
- **TP (validated):** **RING Zn ligand C64** (robust 4/5 runs) — Findlay-LOF and ClinVar-pathogenic; **an AlphaMissense blind spot** (called likely-benign) → the clearest place a functional screen adds value. Plus a coherent **BRCT Findlay-LOF core** (H1686/V1687/D1692, R1751/A1752/R1753, F1734) and the BARD1-interface helix (89–95). Adjacent Zn ligand C61 is borderline (no direct base-editable guide → partial false negative).
- **Likely false positive:** **BRCT-linker 1801–1804** (T1802/G1803/V1804) — BE3D's single strongest cluster, but **Findlay scores these tolerant**. Either a base-editor-specific effect or spurious; do not report as functional without orthogonal support.
- **Novel candidate (hypothesis):** central cluster **214–222** — reproducible multi-guide dropout, but in low-pLDDT disordered structure; unvalidated.
- **False negatives:** LOF residues with no base-editable missense guide in the screen (CBE cannot install the required codon change) are systematically missed — a base-editing coverage limit, not a BE3D failure.

## Competitor context
| Tool | Signal | On BRCA1 | Edge / limit |
|---|---|---|---|
| **BE3D** | BE tiling screen → 3D | RING Zn + BRCT hotspots from function; catches C64 | Functional readout, germline-capable; weak residue resolution, BE-coverage-limited |
| 3dhotspots / cBioPortal | Somatic recurrence | ~nothing (BRCA1 is germline, not recurrently somatic) | Blind to germline tumor suppressors → BE3D's niche |
| AlphaMissense | Seq/structure predictor | Flags RING+BRCT, benign middle; **misses C64** | No assay needed; localized blind spots |
| Findlay SGE | Saturation genome editing | Gold-standard per-residue function, RING+BRCT | Orthogonal truth; only assayed exons |
| ProTiler-Mut | CRISPR tiling → 3D/PPI | (conceptual rival) | Adds PPI/substructure inference |

## Bottom line for VUS work
Use BE3D as **supporting spatial functional evidence (toward ACMG PS3), alongside** Findlay SGE and AlphaMissense — most useful where predictors disagree or are silent (the C64 case). Weight down solitary hotspots in disordered/low-pLDDT regions and any single-guide-driven cluster.
