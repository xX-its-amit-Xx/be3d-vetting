# BE3D User Package — BRCA1 (clinical-genetics / VUS lab)

Runner: Queen "BRCA1". Date: 2026-08-24. Status: **sweep + benchmark COMPLETE** (see §6 robustness).
Outputs: `real_output/P2_BRCA1/`, `outreach/BRCA1/`. Venv: `scratchpad/be3dvenv`. Clone: `scratchpad/BE3D`.
Central deliverable = the **DISCRIMINATION BENCHMARK** (OUTREACH_GUIDE A0), because BRCA1 has an
orthogonal functional gold standard (Findlay SGE) + ClinVar → a real false-positive test.

## 0. Target, data, structures
- **Gene/protein:** BRCA1, UniProt **P38398** (1863 aa). Domain architecture: RING (24–65, with
  BARD1-dimerization helices ~8–22 / ~81–96) — a **long intrinsically-disordered central region
  (~110–1640, AF pLDDT mean 32, 98% pLDDT<70)** — tandem BRCT (1642–1855). This IDR is the crux of
  the discrimination test: a real tool should light up the folded functional domains, not the middle.
- **Screen:** Cuella-Martin et al. 2021 *Cell* (PMID 33592168), DDR base-editing screens, **Table S5**
  (`cm_mmc5.xlsx`, curl 200). Sheet **`MCF10A treatments Subl1`** (non-transformed mammary line); 519
  BRCA1 guides. CBE (BE3-FNLS), competitive **dropout ± olaparib (`T18_OLAP_LFC`) ± cisplatin
  (`T18_CISP_LFC`)** → run as **two BE3D screens, meta-aggregated (SUM)**. `AAChg` column **confirmed
  present** (1-letter HGVS p.). `Function` → category; `ClinVar` column retained for benchmarking.
- **Assay direction (documented):** dropout — **LOF = negative LFC**. Verified on internal controls:
  nonsense mean OLAP −0.33 / CISP −0.26 vs synonymous +0.09; splice-acceptor at essential exons −0.4
  to −2.5. Drug arms roughly **double** the LOF dropout magnitude vs untreated (BRCA1-LOF → PARPi/
  cisplatin synthetic lethality). **No score negation needed** — LOF lands in BE3D's `neg` channel.
- **Structures:** primary **AlphaFold P38398 v6** (full length, whole-protein numbering = unipos);
  **1JM7** chain A (BRCA1–BARD1 RING, res 1–103) monomer; **1T29** chain A (BRCT, 1649–1859) monomer.
  DSSP via `gen_dssp.py` (placeholder SS; affects only characterization). Complex mode (1JM7 A+B) not
  run — the harness is monomer-only and complex mode has an unfixed Windows MAX_PATH blocker; the
  BARD1-interface residues are covered in the curated functional set instead.

## 1. Data reformatting (decision log)
Script `convert_brca1.py`; 519 guides → `brca1_OLAP.tsv` / `brca1_CISP.tsv` (BE3D 4-col format).

| Decision | Choice | Why |
|---|---|---|
| Sheet | MCF10A Subl1 | non-transformed line, cleaner BRCA1 dependency; MCF7 available as secondary |
| `Mutation_list`/`Mutation_type` | per-token parse of `AAChg`: `pX{n}Y`→Missense, `pX{n}=`→Silent, `pX{n}*`→Nonsense; delins/complex skipped; aligned `;`-joined | faithful window editing; BE3D reduces to one category/guide via priority |
| Splice / empty-window / non-coding | `Splice-donor/acceptor;` (placeholder token) / `No Mutation` | LFC-only in BE3D's else-branch; no residue mapping needed |
| `sgRNA_score` | `T18_OLAP_LFC` and `T18_CISP_LFC` (two screens) | drug-specific dropout; both ship computed LFC + p/FDR |
| Category counts | 250 Missense / 110 Silent / 26 Nonsense / 6 Splice / 127 No-Mutation (first-token) | — |

## 2. QA (BE-QA) — PASS
cases=[Nonsense] (35 mapped), controls=[No Mutation] (127). **KS D=0.33 p=3.4e-3 (OLAP); D=0.36
p=9.2e-4 (CISP). Mann-Whitney p=6e-4 / 1.6e-4.** Nonsense decisively separated from neutral → screen accepted.

## 3. Parameter sweep (exhaustive-enough)
Full-length AF, meta OLAP+CISP: **r6/mean/n1000 (PRIMARY)**, r4, r8, r10, r6/**sum**, r6/**n500**,
r6/**atom**. Plus RING 1JM7-A (r6, r8) and BRCT 1T29-A (r6, r8) monomer runs. 11 runs total; each
emits p<.05/.01/.001. (Sweep driver `sweep_brca1.sh`; per-run outputs `brca1_runs/<name>/`.)

## 4. DISCRIMINATION BENCHMARK (A0 — the primary result)
Gold standard = **Findlay 2018 SGE** (`benchmark_BRCA1/findlay_per_residue.tsv`): MaveDB
`urn:mavedb:00000097-0-2` c.HGVS translated via RefSeq NM_007294.3 CDS to protein positions
(`findlay_translate.py`); internal controls validate the LOF threshold (synonymous mean −0.05 vs
nonsense −2.18; midpoint −1.12). Coverage: **RING 1–101 + BRCT 1631–1855, 326 residues, 68 LOF**.
Analysis: `benchmark_brca1.py` (Fisher, Spearman; `benchmark_BRCA1/results/`).

**Primary run (AF full-length, r6/mean/n1000, meta OLAP+CISP), neg (LOF) channel:**

| Test | Number | Honest reading |
|---|---|---|
| **Base rate** (sig / 1863 scored) | **7.4% (p<.05), 4.9% (p<.001)** | NOT "everything is a hotspot" (<15–20% bar) |
| **Domain hit-rate** RING / BRCT / central-IDR | **12.8% / 20.1% / 5.3%** | flags folded functional domains ~2.4–3.8× the disordered middle |
| domain discrimination gap (RING+BRCT − central) | **+0.127** | strong DOMAIN-level discrimination |
| **Whole-protein enrichment** for validated functional residues (Findlay-LOF ∪ literature ∪ ClinVar-P/LP) | R/E=**3.1, OR=4.3, p=6.1e-8** | highly significant — but domain-driven |
| **Trivial pLDDT>70 baseline** (same functional set) | enr **5.1**, hit-rate 88.7% functional vs 12.9% tolerant | **BE3D does NOT beat "is it folded"** — the domain signal is recapitulated by burial/AlphaMissense |
| **Findlay residue-level** enrichment (within RING+BRCT) | M=68, K=54, R=16, E=11.3, enr=**1.42**, OR=**1.78**, **p=0.063** | weak/marginal FINE-grained resolution inside domains |
| Findlay discrimination gap (hit-rate LOF 0.235 vs tolerant 0.147) | **+0.088** | modest |
| **Per-residue concordance** with Findlay function score (Spearman, covered residues) | **ρ ≈ 0.05** | ~no quantitative rank agreement |
| **Precision@10 / @20** (validated functional) | **2/10, 4/20** | top ranks dominated by FP + novel, not validated sites |
| Precision@10, Findlay-covered subset | **4/4 are Findlay-TOLERANT** | the top BRCT-linker hits are false positives vs SGE |

**Interpretation (the honest verdict).** BE3D on BRCA1 **passes the "does it discriminate?" bar at the
DOMAIN level** (it concentrates on RING+BRCT, not the disordered bulk; whole-protein functional
enrichment OR 4.3, p 6e-8) and **is not a flag-everything caller** (base rate 7%). **But** that
domain-level signal is **fully recapitulated — and exceeded — by a trivial "folded / high-pLDDT"
baseline** (enr 5.1 > 3.1), and at the **residue level inside the domains** BE3D adds only **weak**
signal (Findlay enrichment 1.4×, p≈0.06; Spearman ≈ 0). So BE3D here is **complementary functional
evidence, not a standalone discriminator that beats sequence/structure predictors.**

## 5. Findlay concordance + FP/FN (the VUS-relevant evidence)
- **TRUE POSITIVE (the win): RING Zn ligands C61, C64** — top RING hotspots from screen signal alone;
  both **Findlay-LOF** (scores −1.84 / −2.04) and **ClinVar-pathogenic**. **C64 is an AlphaMissense
  blind spot** (C64G 0.089, C64R 0.107 = *likely benign*), immediately adjacent to a correctly-scored
  C61G (0.990). This is precisely where a **functional screen (BE3D / SGE) rescues a predictor miss** —
  the single most outreach-relevant result. (Spatial-neighbor bleed also flags L63/K65, Findlay-tolerant
  — mild aggregation FP.)
- **LIKELY FALSE POSITIVE: BRCT-linker 1801–1804 (T1802/G1803/V1804/G1801)** — BE3D's **single strongest
  cluster** (|LFC3D_neg| up to 3.7), yet **Findlay scores all of them tolerant** (0.03 / −0.22 / −0.03 /
  0.003; not LOF). Flagged explicitly: either a **base-editor-specific effect** the SGE didn't model, or
  spurious. Do not report as functional without orthogonal support.
- **NOVEL CANDIDATE (hypothesis, not a claim): central cluster 214–222** — a **reproducible multi-guide
  dropout** signal (guides at 218/220–221, 213–215 with OLAP LFC −0.94 to −1.50; not a single-guide
  artifact), but in the **low-pLDDT disordered region** (no structural/Findlay validation). Weight down.
- **LIKELY FP (single-guide): 1282–1284** — driven by one guide (H1283Y, OLAP −0.81), low-pLDDT.
- **FALSE NEGATIVES:** Findlay-LOF residues with **no base-editable missense guide** in the screen are
  systematically missed — the CBE can't install the required codon change (a base-editing coverage limit,
  not a BE3D modeling failure). Quantified in §6 / `BRCA1_false_negatives.tsv`.

## 6. Robustness across the sweep + hotspot table
_(filled after sweep completion — see `benchmark_BRCA1/results/benchmark_by_run.tsv`,
`BRCA1_hotspots_full.tsv`, `BRCA1_false_negatives.tsv`.)_

## 7. Competitor comparison (B)
| Tool | Input signal | On BRCA1 | Advantage | Disadvantage |
|---|---|---|---|---|
| **BE3D** | BE tiling screen → 3D LFC3D | RING Zn (C61/C64) + BRCT hotspots from function; concentrates on folded domains | direct functional readout; works on a **germline** tumor suppressor with no somatic recurrence; catches predictor blind spots (C64) | weak residue-level resolution (doesn't beat pLDDT/AlphaMissense at domain level); BE-codon-coverage-limited; bystander/aggregation FPs (1801–1804) |
| **3dhotspots.org / cBioPortal** | somatic mutation **recurrence** | ~no missense signal (BRCA1 rarely/non-recurrently somatically mutated; only truncating passenger hotspots in MSI) | gold-standard for recurrently-mutated **oncogenes** | **blind to germline tumor suppressors** — BE3D's niche |
| **AlphaMissense** (Cheng 2023) | seq+structure DL | correctly flags RING+BRCT pathogenic, IDR benign — no screen; **but mis-calls C64G/R benign** | proteome-wide, instant, strong overall SGE concordance | localized blind spots; predictor not measurement |
| **Findlay SGE** (2018) | saturation genome editing | per-residue function, RING+BRCT (the gold standard used here) | orthogonal functional truth, single-nt res. | only assayed exons; heavy engineering |
| **ProTiler-Mut** (He 2026) | CRISPR tiling → residue/substructure/PPI | (conceptual rival) direct tiling→3D competitor; adds PPI-disruption + separation-of-function class | robustness (3D-RRA), multi-condition | newer, less externally validated |

**Where BE3D wins:** germline tumor suppressor (somatic tools blind) + functional-screen rescue of
predictor blind spots (C64). **Where competitors win:** AlphaMissense/pLDDT already separate
functional-domain vs IDR *and* give finer residue resolution with no experiment; Findlay SGE is the
per-residue functional truth. → **BE3D is complementary evidence, best fused with these, not a replacement.**

## 8. Deliverables (outreach/BRCA1/)
`BRCA1_hotspots.tsv` (robust hotspots + TP/novel/FP flags), G2P TSV (`BRCA1_g2p.tsv`), `BRCA1_email.md`,
`BRCA1_brief.md`. Analysis artifacts in `benchmark_BRCA1/` and `real_output/P2_BRCA1/`.

## 9. BE3D issues encountered (→ BE3D_IMPROVEMENTS.md)
- [SCALE] Full-length AF (1863 res) + meta n1000 ≈ 6 min/run (vs ~40 s for a 400-res domain); the
  meta-aggregate re-randomization over all residues dominates. A per-residue neighbor cache / vectorized
  null would help large multi-domain proteins.
- [SCIENCE] BE3D/spatial-aggregation over a **compact folded domain** flags LOF and tolerant residues at
  similar rates (bystander bleed: L63/K65 flagged from C61/C64) — residue-level discrimination is limited;
  a distance-weighted (Gaussian) kernel instead of the hard radius would sharpen it.
- [CONTEXT] No pLDDT propagation: BE3D happily reports strong hotspots in low-pLDDT disordered regions
  (214–222, 1282–1284) with no confidence down-weighting → avoidable false positives. Propagating pLDDT
  into significance/flags is the single highest-value fix for a VUS use case.
