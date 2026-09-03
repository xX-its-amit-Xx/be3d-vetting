# BE3D — Base-Editing Field Trajectory & Readiness PRs (Phase 7)

**Question:** where is base editing going over the next ~3–5 years, and what would BE3D need to
ingest/represent so a scientist can *plug-and-play* the emerging data? Guiding ethos: **prioritize
what gets treatments to patients faster.** All primary sources cited with authors/year/DOI/PMID.

BE3D today assumes: **CBE/ABE-style** categories (Missense/Silent/Nonsense/Splice-donor/-acceptor/
No-Mutation), a **single per-guide score column**, a **nuclear protein structure**, and a reachability
report (PR #19) hard-wired to **CBE C→T / ABE A→G codon logic**. The field is moving away from every
one of those assumptions at once. Below: the cited survey, then a ranked set of readiness-PR proposals
with the flagship **editor-profile abstraction**.

---

## 1. Editor CHEMISTRIES beyond CBE/ABE — the transversion + multi-substrate wave

The canonical transitions (CBE C·G→T·A; ABE A·T→G·C) are now flanked by editors covering **most
transversions** and **simultaneous multi-base** outcomes. Feasibility is largely solved on paper;
**efficiency, product purity, and AP-site genotoxicity** are the limiters.

- **C→G / C→A glycosylase editors (CGBE / GBE).** Kurt et al. 2021, *Nat Biotechnol* 39:41–46,
  DOI 10.1038/s41587-020-0609-x (CGBE1 = APOBEC(R33A)+eUNG+nCas9; deaminate→excise uracil→AP-site→C→G).
  Zhao et al. 2021, *Nat Biotechnol* 39:35–40, DOI 10.1038/s41587-020-0592-2, PMID 32690970 (GBE:
  C→A in *E. coli*, C→G in mammalian cells). Purity/genotoxicity caveat: AP-site processing can create
  DSBs/deletions (*Nat Cell Biol* 2023, DOI 10.1038/s41556-023-01342-2).
- **Adenine transversions (A→C / A→T = A→Y).** Tong et al. 2023, *Nat Biotechnol* 41:1080–1084,
  DOI 10.1038/s41587-022-01595-6 (AYBE = ABE+MPG glycosylase excising inosine; up to ~72%). Independent
  A→C: Chen et al. 2023 ACBE / high-accuracy **ACBE-Q**, DOI 10.1038/s41587-023-01821-9.
- **Dual / "all-transition" (simultaneous C→T + A→G).** Grünewald et al. 2020 SPACE,
  DOI 10.1038/s41587-020-0535-y; Sakata et al. 2020 Target-ACEmax, DOI 10.1038/s41587-020-0509-0;
  Zhang et al. 2020 A&C-BEmax, DOI 10.1038/s41587-020-0527-y, PMID 32483363.
- **Expanded outcomes (T→G, G→Y), deaminase-free.** Ye et al. 2024 DAF-CBE/DAF-TBE (C→G, **T→G**),
  *Nat Biotechnol* 42:1538–1547, DOI 10.1038/s41587-023-02050-w, PMID 38168994; G→Y editor,
  *Natl Sci Rev* 2023 10:nwad143, DOI 10.1093/nsr/nwad143.
- **TadA becoming a universal, tunable scaffold.** ABE8e/TadA-8e (Richter et al. 2020,
  DOI 10.1038/s41587-020-0453-z); TadA→CBE (Neugebauer et al. 2023, DOI 10.1038/s41587-022-01533-6;
  Lam et al. 2023 **TadCBEs / dual CABE-Ts**, DOI 10.1038/s41587-022-01611-9); TadA→C→G
  (Chen et al. 2023 Td-CGBE, DOI 10.1038/s41587-022-01532-7).
- **PAM-flexibility as the reachability multiplier.** SpG/near-PAMless **SpRY** (Walton et al. 2020,
  *Science* 368:290, DOI 10.1126/science.aba8853); SpCas9-NG (Nishimasu et al. 2018,
  *Science* 361:1259, DOI 10.1126/science.aas9129). PAM flexibility decides *which* pathogenic base is
  addressable — the single biggest lever on BE3D's "reachability" gap.

**Consequence for BE3D:** the reachability report's hard-coded "can CBE/ABE install this?" is now one
row of a much larger matrix. Each editor installs a *different set of nucleotide substitutions* in a
*different window* under a *different PAM/targeting* constraint.

## 2. PRIME EDITING & its convergence with base editing — the "write-anything" layer

Prime editing (PE) installs **all 12 substitutions + small insertions + small deletions**, no DSB, no
donor (Anzalone et al. 2019, *Nature* 576:149, DOI 10.1038/s41586-019-1711-4, PMID 31634902). Engineering
has since made it efficient enough for screens and clinic: PEmax/PE4-5 (Chen et al. 2021, *Cell*
184:5635, DOI 10.1016/j.cell.2021.09.018, PMID 34653350); epegRNAs/tevopreQ1 (Nelson et al. 2022,
DOI 10.1038/s41587-021-01039-7); compact PE6 (Doman et al. 2023, DOI 10.1016/j.cell.2023.07.039); PE7
(Yan et al. 2024, *Nature* 628:639, DOI 10.1038/s41586-024-07259-6). "Write-large" wing: twinPE
(Anzalone et al. 2022, DOI 10.1038/s41587-021-01133-w), PASTE (Yarnall et al. 2023,
DOI 10.1038/s41587-022-01527-4), GRAND (Wang et al. 2022, DOI 10.1038/s41592-022-01399-1).

**Saturation PE screens — a fundamentally different data unit.** The SGE precedent is Findlay et al.
2018 (*Nature* 562:217, DOI 10.1038/s41586-018-0461-z, PMID 30209399). PE now generalizes it donor-free:
Erwood et al. 2022 SVI (*Nat Biotechnol* 40:885, DOI 10.1038/s41587-021-01201-1, PMID 35190686); Gould
et al. 2024 PE "sensor" libraries (DOI 10.1038/s41587-024-02172-9); Cirincione et al. 2024 multiplexed PE
dropout, ~240k epegRNAs, 7,996 nonsense across 1,149 essential genes (DOI 10.1038/s41592-024-02502-4);
Kim et al. 2025 PEER-seq EGFR drug-resistance saturation (DOI 10.1038/s41587-024-02465-z). Head-to-head
BE-vs-PE on one gene: Belli et al. 2024 multimodal EGFR scanning (DOI 10.1038/s41587-024-02439-1,
PMID 39533106).

**The load-bearing data-model difference (this is what BE3D must ingest):**
- **BE tiling = one score per guide *window*.** A guide edits any eligible base in a ~4–8-nt window →
  a *distribution* of outcomes (intended + bystander + mixed missense/silent/nonsense), restricted to
  the editor's chemistry. BE3D's `Mutation_type`/`Mutation_list` reflect exactly this many-outcomes-per-
  guide model.
- **PE/SGE saturation = one score per *precise variant/allele*.** A pegRNA templates one defined edit
  (specific substitution, insertion, or deletion at a specific codon), so the score maps to a single
  allele — no bystander deconvolution, and edit types (transversion, indel, frameshift, stop-gain
  anywhere) that BE tiling cannot reach.

## 3. ORGANELLE & non-nuclear editing — no PAM, no sgRNA, protein-guided dsDNA

- **DdCBE** (Mok et al. 2020, *Nature* 583:631, DOI 10.1038/s41586-020-2477-4, PMID 32641830): split
  DddA toxin + paired **TALEs** + UGI → C→T in mtDNA. **No guide RNA, no PAM**; targeting is paired
  protein–DNA half-sites with a spacer; strong **5′-TC** context. DddA6/DddA11 broaden context
  (Mok et al. 2022, DOI 10.1038/s41587-022-01256-8, PMID 35379962).
- **TALED** (Cho et al. 2022, *Cell* 185:1764, DOI 10.1016/j.cell.2022.03.039, PMID 35472301): mito
  **A→G**. **mitoBEs** strand-selective (Yi et al. 2024, DOI 10.1038/s41587-023-01791-y). Nuclear
  off-target caveat: Lei et al. 2022, *Nature* 606:804, DOI 10.1038/s41586-022-04836-5. Plastid editing:
  Li et al. 2021, *Mol Plant*, DOI 10.1016/j.molp.2021.07.007.

**Consequence for BE3D:** targets are the **13 mtDNA-encoded OXPHOS subunits** (+ rRNA/tRNA) on a
**heteroplasmy**-aware coordinate system; "guides" are TALE protein designs, not an sgRNA library, and
reachability is gated by **paired TALE half-sites + spacer + TC context**, not a PAM offset.

## 4. RNA base editing — transient, transcript-level, reversible

REPAIR (A→I read as G; Cox et al. 2017, *Science* 358:1019, DOI 10.1126/science.aaq0180, PMID 29070703);
RESCUE (C→U; Abudayyeh et al. 2019, *Science* 365:382, DOI 10.1126/science.aax7063, PMID 31296651);
endogenous-ADAR recruitment LEAPER (Qu et al. 2019, DOI 10.1038/s41587-019-0178-z), RESTORE (Merkle et
al. 2019, DOI 10.1038/s41587-019-0013-6), CLUSTER (Reautschnig et al. 2022,
DOI 10.1038/s41587-021-01105-0). **First-in-human RNA editing:** Wave WVE-006 (AIMer recruiting ADAR to
edit SERPINA1 Z-mRNA for AATD; RestorAATion-2, 2024). **Consequence for BE3D:** edits the *transcript
pool*, not the genome — dose/time-dependent, reversible, ADAR-expression-gated; permanent-genotype
LFC3D semantics do not directly apply and must be flagged.

## 5. READOUTS — from one LFC to a per-variant phenotype *vector*

- **Single-cell BE screens:** BE-Perturb-seq (Martin-Rufino et al. 2023, *Cell* 186:2456,
  DOI 10.1016/j.cell.2023.03.035, PMID 37044066); STING-seq/beeSTING-seq (Morris et al. 2023, *Science*
  380:eadh7699, DOI 10.1126/science.adh7699).
- **Multi-marker FACS:** Schmidt et al. 2023, *Nature* 624:390, DOI 10.1038/s41586-023-06835-6
  (proliferation + cytokine + cytotoxicity per allele); condition panels: Cuella-Martin et al. 2021,
  *Cell* 184:1081, DOI 10.1016/j.cell.2021.01.041, PMID 33606978 (four DNA-damaging agents → one profile
  per variant).
- **Protein abundance/stability:** VAMP-seq (Matreyek et al. 2018, *Nat Genet* 50:874,
  DOI 10.1038/s41588-018-0122-z, PMID 29785012).

**Consequence for BE3D:** a variant/guide increasingly carries **N score columns** (expression programs,
cell-state fractions, multiple markers, abundance, drug conditions), each with its own metadata — not a
scalar fitness LFC. BE3D's single `sgRNA_score`/`val_col` is the bottleneck.

## 6. CLINICAL TRAJECTORY — what's actually reaching patients, and the edit-type logic

In vivo/ex vivo editing crossed into pivotal-stage in 2025–2026. The edit types reaching patients map to
three classes (this is the "patient-impact" lens for every PR below):

| Target class | Edit type | Lead program | Status / source |
|---|---|---|---|
| Gene **knockdown** (install stop / disrupt splice) | ABE | VERVE-102 (PCSK9), VERVE-201 (ANGPTL3) | Human POC, LDL −62%, durable 18 mo — *NEJM* 2026 DOI 10.1056/NEJMoa2601283 |
| **Beneficial-gene de-repression** (promoter edit) | ABE | BEAM-101 (HBG1/2→HbF, SCD) | NEJM-published, RMAT — *NEJM* 2026;394:1824 DOI 10.1056/NEJMoa2504835 |
| **Precise point correction** | in vivo ABE | BEAM-302 (SERPINA1 E342K, AATD) | First in-human BE *correction*; AAT >protective threshold (Beam 2025–26; ERS 2026) |
| **Precise indel/arbitrary correction** | ex vivo **prime** | PM359 (NCF1 ΔGT, CGD) | First-in-human prime editing — *NEJM* DOI 10.1056/NEJMoa2509807 |
| **Multiplex T-cell KO** | ex vivo CBE ×4 | BEAM-201 (CD7/TRAC/PDCD1/CD52) | Early Phase 1/2 |

Delivery: liver-tropic **LNP** (Musunuru et al. 2021, *Nature* 593:429, DOI 10.1038/s41586-021-03534-y,
PMID 34012082; Rothgangl et al. 2021, DOI 10.1038/s41587-021-00933-4; GalNAc-LNP for LDLR-independent
uptake); dual-AAV split-intein for size; **eVLPs** (Banskota et al. 2022, *Cell* 185:250,
DOI 10.1016/j.cell.2021.12.021). Field-wide gating factor: **transient LNP hepatotoxicity** (Intellia
nex-z grade-4 liver AE / FDA hold Oct 2025, lifted 2026) — the main constraint on "faster to patients."

**Modality split now visible in clinic:** transition point-mutation in an editable window → ABE/CBE;
transversion → CGBE/AYBE; indel or otherwise-arbitrary fix → prime. **This is precisely the decision an
editor-profile abstraction would let BE3D make** ("your validated hotspot needs A→C — unreachable by
ABE, reachable by AYBE/ACBE-Q or prime editing"), turning BE3D from a CBE/ABE-only tool into a
modality-selection aid that stops teams chasing unreachable residues.

---

## 7. Readiness-PR proposals (ranked)

Ranked by patient-impact × evidence strength × additivity/testability. Effort S/M/L. All are additive
and unit-testable against BE3D's existing pipeline. The **editor-profile abstraction is the flagship**:
it is the smallest change that unlocks the whole chemistry + prime-editing + clinical-modality story,
and it directly generalizes the reachability report already shipped in PR #19.

| Rank | Idea | Motivating citation(s) | BE3D gap | Patient-impact rationale | Effort | Additive & testable? |
|---|---|---|---|---|---|---|
| 1 (flagship) | **EDITOR-PROFILE abstraction** — a declarative registry of installable nucleotide substitutions per editor (CBE C→T, ABE A→G, CGBE C→G, GBE C→A, AYBE A→C/A→T, ACBE-Q A→C, dual A&C, TadCBE C→T, DAF-TBE T→G, G→Y, PE=arbitrary+indel), each with window + targeting(PAM/TALE) + purity note; the reachability report (#19) becomes `is_reachable(from_aa,to_aa,editor)` parameterized over the registry | Kurt 2021 (10.1038/s41587-020-0609-x); Zhao 2021 (10.1038/s41587-020-0592-2); Tong 2023 (10.1038/s41587-022-01595-6); Chen 2023 ACBE (10.1038/s41587-023-01821-9); Grünewald 2020 (10.1038/s41587-020-0535-y); Lam 2023 (10.1038/s41587-022-01611-9); Walton 2020 (10.1126/science.aba8853); Anzalone 2019 (10.1038/s41586-019-1711-4) | Reachability + all editor logic hard-wired to CBE/ABE codon rules; no way to ask "which editor reaches this residue" | Directly the clinical modality split (BEAM-302 ABE-transition vs PM359 prime-indel): tells a team which editor class can install a needed fix, so they stop chasing unreachable residues and pick the right modality faster — the #1 time-to-patient lever BE3D already identified | M | Yes — pure lookup + registry; extends #19's 2 tests |
| 2 | **Configurable mutation-category vocabulary + arbitrary edit types** — data-driven category map (replace hard-coded Missense/Silent/Nonsense/Splice-donor/-acceptor) + first-class **indel / frameshift / stop-gain-anywhere** categories; warn on unmatched values | Anzalone 2019 (10.1038/s41586-019-1711-4); Erwood 2022 (10.1038/s41587-021-01201-1); Cirincione 2024 (10.1038/s41592-024-02502-4); Gould 2024 (10.1038/s41587-024-02172-9); already flagged in BE3D_IMPROVEMENTS (hard-coded Splice-donor/-acceptor silently empties a bare "Splice" bucket) | Categories/priorities hard-coded; prime-editing indels/frameshifts have no category and are silently dropped | Lets BE3D ingest prime-editing screens (the clinically-validated arbitrary-correction modality, PM359) without misclassifying edits — extends interpretation to the fastest-growing screen type | S–M | Yes — config + validation; small unit tests |
| 3 | **Per-variant / allele-resolved ingestion mode (PE/SGE)** — an input path where the unit is a *defined protein variant* (position + from→to + score), bypassing guide-window/bystander deconvolution; maps one score per allele onto structure | Findlay 2018 (10.1038/s41586-018-0461-z); Erwood 2022 (10.1038/s41587-021-01201-1); Kim 2025 (10.1038/s41587-024-02465-z); Belli 2024 (10.1038/s41587-024-02439-1) | Parser assumes guide→window→multiple categories; can't represent one score per precise allele | Saturation PE/SGE maps are the emerging gold-standard variant-effect data (BRCA1, EGFR-TKI resistance); ingesting them lets BE3D 3D-interpret clinical VUS maps directly | M | Yes — alternate parser feeding existing LFC3D; testable on a Findlay-style TSV |
| 4 | **Multi-phenotype / multi-score-column input** — accept N score columns per variant/guide; run LFC3D per column; emit per-phenotype hotspots + a cross-phenotype comparison | Martin-Rufino 2023 (10.1016/j.cell.2023.03.035); Schmidt 2023 (10.1038/s41586-023-06835-6); Cuella-Martin 2021 (10.1016/j.cell.2021.01.041); Matreyek 2018 VAMP-seq (10.1038/s41588-018-0122-z); Morris 2023 (10.1126/science.adh7699) | Single `sgRNA_score`/`val_col`; single-cell/FACS/VAMP-seq/condition-panel screens emit a phenotype vector | Multi-condition screens (e.g. drug panels) reveal *resistance vs essentiality vs stability* geometry per residue — separates on-pathway from off-pathway hotspots, cutting wasted wet-lab validation | M–L | Yes — loop existing pipeline over columns; testable |
| 5 | **Organelle / RNA editor caveat & targeting-model layer** — editor-class metadata + guardrails: mito editors = TALE-paired, no PAM/sgRNA, TC-context, heteroplasmy, 13 OXPHOS targets; RNA editing = transient/transcript-level/reversible/ADAR-gated → warn that genomic LFC3D semantics differ; reachability uses TALE half-site+spacer model, not PAM | Mok 2020 (10.1038/s41586-020-2477-4); Cho 2022 (10.1016/j.cell.2022.03.039); Yi 2024 (10.1038/s41587-023-01791-y); Cox 2017 (10.1126/science.aaq0180); Abudayyeh 2019 (10.1126/science.aax7063) | No concept of non-PAM/TALE targeting or transient RNA edits; would silently mis-apply PAM reachability & permanent-genotype assumptions | Prevents silent misinterpretation for mito-disease (OXPHOS) and RNA-editing (AATD, WVE-006) programs — correctness guardrail that protects downstream clinical decisions | S–M | Yes — metadata + warnings; cheap tests |

```json
[
 {"rank":1,"idea":"EDITOR-PROFILE abstraction: declarative registry of installable substitutions per editor (CBE C>T, ABE A>G, CGBE C>G, GBE C>A, AYBE A>C/A>T, ACBE-Q A>C, dual A&C, TadCBE C>T, DAF-TBE T>G, G>Y, prime=arbitrary+indel) with window/PAM-or-TALE/purity metadata; reachability report parameterized over the registry (generalizes PR #19)","citations":["Kurt 2021 10.1038/s41587-020-0609-x","Zhao 2021 10.1038/s41587-020-0592-2 PMID 32690970","Tong 2023 10.1038/s41587-022-01595-6","Chen 2023 ACBE 10.1038/s41587-023-01821-9","Grunewald 2020 10.1038/s41587-020-0535-y","Lam 2023 TadCBE 10.1038/s41587-022-01611-9","Walton 2020 SpRY 10.1126/science.aba8853","Anzalone 2019 10.1038/s41586-019-1711-4"],"gap":"Reachability and editor logic hard-wired to CBE/ABE codon rules; no way to query which editor can reach a residue","patient_impact":"Mirrors the clinical modality split (ABE-transition BEAM-302 vs prime-indel PM359); tells teams which editor class installs the needed fix so they stop chasing unreachable residues and pick the right modality faster","effort":"M","pr_candidate":true},
 {"rank":2,"idea":"Configurable mutation-category vocabulary + first-class arbitrary edit types (indel/frameshift/stop-gain-anywhere); warn on unmatched categories","citations":["Anzalone 2019 10.1038/s41586-019-1711-4","Erwood 2022 10.1038/s41587-021-01201-1 PMID 35190686","Cirincione 2024 10.1038/s41592-024-02502-4","Gould 2024 10.1038/s41587-024-02172-9"],"gap":"Mutation categories/priorities hard-coded (Missense/Silent/Nonsense/Splice-donor/-acceptor); prime-editing indels/frameshifts have no category and are silently dropped","patient_impact":"Enables ingestion of prime-editing screens (the clinically validated arbitrary-correction modality) without misclassifying edits","effort":"S-M","pr_candidate":true},
 {"rank":3,"idea":"Per-variant / allele-resolved ingestion mode for PE/SGE saturation screens: unit = defined protein variant (pos+from>to+score), bypassing guide-window bystander deconvolution","citations":["Findlay 2018 10.1038/s41586-018-0461-z PMID 30209399","Erwood 2022 10.1038/s41587-021-01201-1","Kim 2025 PEER-seq 10.1038/s41587-024-02465-z","Belli 2024 10.1038/s41587-024-02439-1 PMID 39533106"],"gap":"Parser assumes guide->window->multiple categories; cannot represent one score per precise allele","patient_impact":"Saturation PE/SGE variant-effect maps (BRCA1, EGFR-TKI) are emerging clinical VUS gold standard; ingesting them lets BE3D 3D-interpret them directly","effort":"M","pr_candidate":true},
 {"rank":4,"idea":"Multi-phenotype / multiple-score-column input: N score columns per variant/guide, LFC3D per column, per-phenotype hotspots + cross-phenotype comparison","citations":["Martin-Rufino 2023 BE-Perturb-seq 10.1016/j.cell.2023.03.035 PMID 37044066","Schmidt 2023 10.1038/s41586-023-06835-6","Cuella-Martin 2021 10.1016/j.cell.2021.01.041 PMID 33606978","Matreyek 2018 VAMP-seq 10.1038/s41588-018-0122-z PMID 29785012","Morris 2023 STING-seq 10.1126/science.adh7699"],"gap":"Single sgRNA_score/val_col; single-cell/FACS/VAMP-seq/condition-panel screens emit a phenotype vector","patient_impact":"Multi-condition screens separate resistance vs essentiality vs stability geometry per residue, cutting wasted wet-lab validation","effort":"M-L","pr_candidate":true},
 {"rank":5,"idea":"Organelle/RNA editor caveat & targeting-model layer: editor-class metadata + guardrails for TALE-paired no-PAM mito editors (TC context, heteroplasmy, 13 OXPHOS targets) and transient/reversible ADAR-gated RNA editing","citations":["Mok 2020 DdCBE 10.1038/s41586-020-2477-4 PMID 32641830","Cho 2022 TALED 10.1016/j.cell.2022.03.039 PMID 35472301","Yi 2024 mitoBE 10.1038/s41587-023-01791-y","Cox 2017 REPAIR 10.1126/science.aaq0180 PMID 29070703","Abudayyeh 2019 RESCUE 10.1126/science.aax7063 PMID 31296651"],"gap":"No concept of non-PAM/TALE targeting or transient RNA edits; PAM reachability and permanent-genotype assumptions applied silently","patient_impact":"Prevents silent misinterpretation for mito-disease (OXPHOS) and RNA-editing (AATD/WVE-006) programs; correctness guardrail protecting clinical decisions","effort":"S-M","pr_candidate":true}
]
```

## 8. Bottom line

Base editing is fanning out along four axes simultaneously — **chemistry** (transversions + multi-
substrate TadA + PAM-flex), **prime-editing convergence** (arbitrary edits, indels, per-variant
saturation screens), **compartment** (mito TALE editors, transient RNA editing), and **readout**
(single-cell/FACS/VAMP-seq phenotype vectors) — while the **clinic** validates a clean modality split
(ABE/CBE for editable transitions, CGBE/AYBE for transversions, prime for indels/arbitrary). BE3D's
CBE/ABE-only, single-score, PAM-gated, permanent-genotype assumptions are each contradicted by one of
these trends. The **editor-profile abstraction** is the flagship because it is the smallest, most
additive change that reframes BE3D's already-shipped reachability report (PR #19) into a general
modality-selection engine — the most direct "gets treatments to patients faster" lever in the set —
and every other PR (arbitrary categories, per-variant ingestion, multi-score, organelle/RNA guardrails)
composes cleanly on top of it.
