# BE3D — Personas & Use Cases, v2 (Emerging-Tech Gap-Fill)

**Research Queen "PERSONAS-V2".** This EXTENDS `Q2_personas.md` (personas A1–A11, I1–I9 + matrix). It does
**not** repeat Q2 — it adds **15 NEW personas** for emerging base-/prime-/RNA-/organelle-editing and the
therapeutic areas Q2 under-covered, plus **7 NEW use cases**, an updated fit matrix, and a cited
quantification of the incremental *patient* population these personas serve.

**Guiding ethos (user):** judge every fit by *"gets treatments to patients faster."* The high-value
additions below are the **therapeutic-interpretation** personas (in-vivo BE target-residue selection,
saturation-SNV → 3D for VUS, multiplex cell-therapy edit-site choice).

**Positioning reminder (from BE3D_BRIEF).** BE3D turns a *completed BE tiling screen* into **3D
structural hotspots** (LFC3D + spatial clustering + multi-screen meta-aggregation → G2P). It is a
**downstream interpretation layer** over *coding-sequence* screens mapped to protein structure. Two
structural facts gate several new personas:
1. **It maps to PROTEIN structure** → non-coding/regulatory edits (promoter/enhancer) are out of native scope.
2. **It assumes CBE/ABE-style installable substitutions** → new editor chemistries (DdCBE, TALED, prime,
   ADAR-RNA) need the **editor-profile abstraction** (the flagship Phase-7 readiness PR: installable
   substitutions per editor — CBE C>T, ABE A>G, CGBE C>G, AYBE transversion, prime arbitrary, mito
   dual-strand, ADAR A>I). Where a persona needs that (or non-coding / multi-phenotype ingestion), it is
   flagged **needs_readiness_feature**; otherwise it **works today**.

---

## 1. NEW PERSONA CARDS

Card format matches Q2: role · goals · input data (and BE3D-ingestibility today vs readiness) · why BE3D
helps · alternatives · decision criteria · sophistication · adoption likelihood.

### NEW ACADEMIC PERSONAS

---
**A12 — Mitochondrial-disease / organelle-editing researcher (DdCBE / TALED)**
- **Role:** Studies mtDNA-encoded OXPHOS disease; uses **DdCBE** (DddA-derived dsDNA cytosine editors, Liu/Broad + Mougous, 2020) and **TALED** (TALE-linked adenine deaminase, Kim/IBS, 2022) to install/correct mtDNA point mutations.
- **Goals:** Model pathogenic mtDNA variants (MELAS m.3243A>G, LHON, Leigh) and read out which residues of the 13 mtDNA-encoded OXPHOS subunits (ND1–6, COX1–3, ATP6/8, CYB) are functionally critical.
- **Input:** DdCBE/TALED tiling of mt-encoded proteins → per-residue functional/heteroplasmy-shift LFC. **BE3D today:** the 13 mt-proteins DO have structures (AlphaFold + cryo-EM complexes I/III/IV/V), so LFC3D mapping is conceptually native — **BUT** needs the **editor-profile** for DdCBE's *dual-strand* C•G→T•A window and the **mitochondrial genetic code** + heteroplasmy handling → **needs_readiness_feature**.
- **Why BE3D:** Cryo-EM DdCBE structures already guide a 2–3-nt editing window; LFC3D over the OXPHOS complex interfaces would surface which edited residues break assembly vs catalysis.
- **Alternatives:** Manual PyMOL on complex structures; bespoke heteroplasmy analysis.
- **Decision criteria:** Non-nuclear codon table; dual-strand reachability; organelle-complex (multi-chain) support.
- **Sophistication:** High. **Adoption: MED** — distinctive, small but fast-growing field; conversion gated on editor-profile readiness. [DdCBE: broadinstitute.org/news/new-molecular-tool; TALED: Cell 2022 S0092867422003890, correcting ~43% of 90 pathogenic mtDNA mutations sciencedaily 2022/04; DdCBE cryo-EM: Nature Commun s41467-023-42359-3]

---
**A13 — Prime-editing / saturation-SNV variant-mapping lab**
- **Role:** Runs **saturation prime editing** (PE tiling; Anzalone/Liu lineage) to install *arbitrary* SNVs across a clinically actionable gene and score each for function — the highest-resolution VUS assay.
- **Goals:** Resolve VUS (>1M in ClinVar, 2024) with dense per-variant maps; add spatial context.
- **Input:** Saturation-PE per-variant functional scores (e.g. PEER-seq scored 2,476 EGFR SNVs, ~99% of the kinase domain). **BE3D today:** the data are per-residue coding effects on a structured protein → LFC3D/clustering applies with only a **format adapter + prime editor-profile** (arbitrary substitutions) → **needs_readiness_feature (minor/format-level)**.
- **Why BE3D:** Adds 3D spatial corroboration + hotspot membership to per-variant PE maps (ACMG PS3 spatial support); rescues weak singletons via neighborhood aggregation.
- **Alternatives:** Per-variant heatmaps only; AlphaMissense/ESM overlay; direct SGE maps (A5's turf).
- **Decision criteria:** Ingest dense saturation tables; distinguish +/- effect; clinical defensibility.
- **Sophistication:** Very high. **Adoption: MED-HIGH** — PE is the fastest-rising VUS modality; strong pull, low switching cost once the adapter exists. [Prime editing: Nature Biotech s41587-020-0455-x; saturation-PE VUS: Cell Genomics S2666-979X(25)00070-9; PEER-seq EGFR: PMC12008803; drug-resistance saturation-PE: Nat Biotech s41587-024-02465-z]

---
**A14 — Single-cell / Perturb-BE multi-phenotype screener**
- **Role:** Couples base editing with single-cell readouts (Perturb-seq-style; scBE; PTM-centric BE screens, e.g. phosphosite editing during T-cell activation, Nat Methods 2024).
- **Goals:** Get *many* phenotypes per edited residue (transcriptome/epigenome/signaling), not one dropout LFC → per-phenotype functional hotspots.
- **Input:** Per-cell, per-phenotype effect matrices. **BE3D today:** BE3D ingests **one** score column per residue → multi-phenotype requires running it N times or a **multi-readout ingestion** feature → **needs_readiness_feature**.
- **Why BE3D:** Would let each phenotype axis produce its own 3D hotspot map, then MetaClust3D across phenotypes → residues that are functional-surface-specific vs pleiotropic.
- **Alternatives:** Per-phenotype 1D analysis; custom scripts.
- **Decision criteria:** Multi-readout ingestion; per-phenotype null models.
- **Sophistication:** Very high. **Adoption: MED** — cutting-edge but small cohort; gated on multi-phenotype readiness. [Perturb-seq; PTM-BE screens Nat Methods s41592-024-02256-z; Multiome Perturb-seq Cell Systems S2405-4712(24)00366-1]

---
**A15 — Neuroscience / neurodegeneration in-vivo BE lab**
- **Role:** Uses AAV/LNP base editing for neurological disease — **PRNP** (prion), **SOD1** (ALS), **APOE**, **MAPT**.
- **Goals:** Choose and validate the *right* coding edit (e.g. PRNP R37X knockdown) and understand which residues of a neuro target are functionally load-bearing.
- **Input:** BE tiling of a neuro gene → per-residue LFC. **BE3D today:** SOD1, PRNP, APOE are small, well-structured coding proteins → **works today** (CBE/ABE, AlphaFold/PDB native).
- **Why BE3D:** Maps functional hotspots onto the target to distinguish knockdown-tolerant vs structurally essential residues; supports the "which edit is safest/most effective" decision.
- **Alternatives:** Manual structure inspection; single-variant validation.
- **Decision criteria:** Small-protein handling; interpretable figures for grant/paper.
- **Sophistication:** High. **Adoption: MED-HIGH** — high-profile Broad prion result (AAV base editing, ~50% brain PrP reduction, ~52% lifespan extension in prion mice, Nat Med Jan 2025) makes this a visible, patient-facing wedge. [Nat Med s41591-024-03466-w; broadinstitute.org/news/gene-editing-extends-lifespan-prion]

---
**A16 — Immuno-oncology / T-cell engineering academic lab (multiplex BE in immune cells)**
- **Role:** Academic groups (e.g. UCL/GOSH lineage) making **multiplex base-edited CAR-T** — the Alyssa world-first used CBE stop-codons to disrupt CD7, TRAC, CD52 in donor "universal" CAR7 cells (2022).
- **Goals:** Pick multiplex knockout **edit sites** that fully disrupt the target protein while minimizing collateral/off-target; understand which residues a premature stop must reach.
- **Input:** BE tiling / knockout-efficiency screens of immune-modulator genes (TRAC, CD52, PDCD1, B2M, CD7). **BE3D today:** *single-protein* disruptive-residue mapping (which stop/missense edits kill function) **works today**; *cross-gene multiplex* prioritization is beyond current single-protein scope → **partial / needs_readiness_feature for the multi-gene layer**.
- **Why BE3D:** For each target, LFC3D shows which coding positions are structurally essential → the highest-confidence knockout window (base editing avoids DSB translocations — ~210-fold fewer than Cas9 multiplexing).
- **Alternatives:** Guide-efficiency predictors (BE-Hive) + manual stop-codon placement.
- **Decision criteria:** Reachability of a disruptive edit; multi-gene handling.
- **Sophistication:** High. **Adoption: MED**. [Alyssa/CAR7: UCL news 2022/dec; NCT05397184; multiplex translocation risk: PMC10123993]

---
**A17 — Microbiome / bacterial base-editing lab**
- **Role:** High-throughput **bacterial** base editing (E. coli, C. glutamicum; phage-delivered CBEs) for AMR studies and microbiome engineering.
- **Goals:** Map which residues of a bacterial enzyme (e.g. a resistance gene, essential enzyme) are functionally critical from genome-scale BE screens.
- **Input:** Bacterial BE tiling → per-residue LFC. **BE3D today:** bacterial proteins are covered by AlphaFold; custom FASTA+PDB path is native → **works today** (with standard bacterial codon table).
- **Why BE3D:** Structure-aware hotspots on essential/resistance enzymes; cross-strain conservation via built-in alignment.
- **Alternatives:** 1D essentiality maps; bespoke scripts.
- **Decision criteria:** Non-human protein support; alignment robustness.
- **Sophistication:** High. **Adoption: LOW-MED** — genome-wide bacterial BE screens still rare but rising. [C. glutamicum BE screen: NAR 52/7/4079; E. coli ScCas9 CBE: PMC11039988; phage-delivered CBE AMR: PNAS 2206744119]

---
### NEW INDUSTRY PERSONAS (therapeutic — highest patient impact)

---
**I10 — Hemoglobinopathy ex-vivo BE program (Beam BEAM-101)**
- **Role:** SCD/β-thalassemia programs reactivating fetal hemoglobin. **BEAM-101** base-edits the **HBG1/HBG2 promoter** to block BCL11A binding (BEACON: all 17 patients HbF >60%, HbS <40%; RMAT + orphan).
- **Goals:** Interpret/optimize the edit and its trans-factor consequences.
- **Input:** Promoter/regulatory edits. **BE3D today: LIMITED / needs_readiness_feature** — the *edit site is NON-CODING* (a BCL11A-binding motif), so protein-structure mapping does not apply to the edit itself. The mappable biology is the **BCL11A protein–DNA interface** (BCL11A zinc-finger domain has structure); a BCL11A-coding tiling screen would map today, but the promoter edit does not.
- **Why BE3D:** Only via a BCL11A protein-domain functional screen, not the therapeutic promoter edit.
- **Alternatives:** Regulatory/motif analysis tools; not BE3D's lane.
- **Decision criteria:** Would require non-coding/protein–DNA readiness.
- **Sophistication:** Very high. **Adoption: LOW** — flagship therapeutic but structurally off-scope. [BEAM-101 BEACON: investors.beamtx EHA 2025; Casgevy BCL11A enhancer comparator: crisprtx FDA approval Dec 2023]

---
**I11 — Cardiometabolic in-vivo BE program (Verve → Eli Lilly: PCSK9 / ANGPTL3 / LPA)**
- **Role:** In-vivo liver base editing to permanently lower LDL-C / remnant cholesterol / Lp(a). **VERVE-102** (PCSK9) gave mean 53%/max 69% LDL-C reduction (Heart-2 Ph1b); VERVE-201 (ANGPTL3), VERVE-301 (LPA); **Lilly acquired Verve for up to $1.3B (June 2025)**.
- **Goals:** Choose the **therapeutic edit residue** — which coding position, when knocked out, safely recapitulates a *protective human LOF variant* (PCSK9 Y142X/C679X → 88% lower CHD risk).
- **Input:** BE tiling of PCSK9/ANGPTL3 + structures. **BE3D today: WORKS TODAY** — these are coding, well-structured proteins; LFC3D over a target-tiling screen maps functional/knockout hotspots onto the structure for edit-site prioritization (native CBE/ABE).
- **Why BE3D:** Directly serves **U15 in-vivo therapeutic target-residue prioritization** — reduce wasted wet-lab validation by ranking candidate installable knockout edits by structural consequence.
- **Alternatives:** Internal validated platforms; human-genetics LOF catalogs.
- **Decision criteria:** Data security (local mode), regulatory reproducibility; likely borrows the *method*.
- **Sophistication:** Very high. **Adoption: MED** (high value, big-pharma internal tooling bias). [Lilly/Verve $1.3B: cnbc 2025/06/17; VERVE-102 Heart-2: vervetx.gcs-web; PCSK9 LOF protection: PMC5729040]

---
**I12 — Liver / AATD in-vivo BE program (Beam BEAM-302, SERPINA1 E342K)**
- **Role:** One-time A-to-G correction of the **PiZ (E342K)** SERPINA1 mutation. **BEAM-302** Ph1/2 (Mar 2025): first-ever clinical genetic correction of the causal mutation; 60 mg → 2.8× functional AAT, up to 78% Z-AAT reduction. Prime Medicine has a parallel prime-editing AATD program.
- **Goals:** Understand the single-residue structural consequence (E342K drives serpin polymerization/ER retention) and screen SERPINA1 for polymerization-modifying residues.
- **Input:** SERPINA1 tiling / polymerization screen + serpin structure. **BE3D today: WORKS TODAY** — SERPINA1 is a well-structured coding protein; the causal target is a *single defined residue* (Glu342Lys), ideal for structural hotspot interpretation.
- **Why BE3D:** Maps which residues govern the folding/polymerization pathology → prioritize edit and anticipate liability residues.
- **Alternatives:** Serpin biophysics/MD; single-variant assays.
- **Decision criteria:** Custom liganded/mutant structure; interpretability.
- **Sophistication:** Very high. **Adoption: MED** — clean single-residue serpin story; strong exploratory fit. [BEAM-302 data: globenewswire 2025/03/10; E342K polymer mechanism: PMC4957051; Prime Medicine AATD: investors.primemedicine 2025/03]

---
**I13 — Allogeneic CAR-T / multiplex-BE cell-therapy company**
- **Role:** Off-the-shelf cell therapy via multiplex editing. **Beam BEAM-201** = first **quadruplex** base-edited allogeneic anti-CD7 CAR-T (KO of TRAC, CD52, PDCD1, CD7; NCT05885464, ~102 R/R T-ALL/T-LL; first base-editing CAR-T dosed in US). Peers: Caribou CB-010 (PD-1 KO), Prime Medicine PASSIGE (B2M/TRAC + knock-in).
- **Goals:** Pick multiplex knockout edit sites that fully null each protein with minimal collateral; base editing chosen precisely to avoid DSB translocations.
- **Input:** Per-gene BE knockout screens. **BE3D today:** single-protein disruptive-residue mapping **works today**; **multi-gene cross-target prioritization → needs_readiness_feature** (U12).
- **Why BE3D:** For each of TRAC/CD52/PDCD1/CD7, LFC3D flags the structurally essential positions → most reliable stop-codon window.
- **Alternatives:** Guide-efficiency predictors; empirical KO screens.
- **Decision criteria:** Reachability, multi-gene batch mode, data security.
- **Sophistication:** Very high. **Adoption: MED**. [BEAM-201: investors.beamtx; ASH 2024 Blood 144/Suppl1/4838; Caribou CB-010: cgtlive; Prime PASSIGE: primemedicine ASH 2023]

---
**I14 — RNA-editing (ADAR, A-to-I) therapeutics company**
- **Role:** Programmable **A-to-I RNA editing** via endogenous/engineered ADAR — recodes a codon at the transcript level. Players: **Wave (WVE-006, AATD — first-in-human therapeutic RNA editing, Sept 2025)**, **Korro (KRRO-110, AATD)**, **AIRNA (AIR-001, $155M B)**, ProQR/Axiomer, **Shape**, Ascidian.
- **Goals:** Understand the protein-level consequence of the recoded residue (e.g. SERPINA1 Z→M restoration; ADAR recoding = a defined missense change).
- **Input:** RNA-editing recoding maps. **BE3D today: needs_readiness_feature** — requires an **ADAR editor-profile (A>I≈A>G recoding)** + RNA-level ingestion; but because recoding lands as a *protein missense change on a structured protein*, once ingested the LFC3D machinery applies.
- **Why BE3D:** Structural interpretation of which recoded residues restore/impair function (AATD programs are validated at protein level: M-AAT vs Z-AAT).
- **Alternatives:** Transcript/protein assays; no structural layer.
- **Decision criteria:** RNA-editing ingestion; recoding→residue mapping.
- **Sophistication:** Very high. **Adoption: LOW-MED** (gated on RNA readiness). [Wave WVE-006: ir.wavelifesciences RestorAATion-2; Korro KRRO-110: biospace; AIRNA $155M: genengnews; market ~$195M 2025: natlawreview]

---
**I15 — Mitochondrial-medicine biotech (DdCBE/TALED; Pretzel Therapeutics)**
- **Role:** Company translating organelle editing / heteroplasmy-shifting. **Pretzel Therapeutics** ($72.5M Series A, ARCH/Mubadala; founders incl. Minczuk, Gustafsson) — mtDNA correction platform.
- **Goals:** Same structural mapping as A12 but in a therapeutic pipeline (which mt-protein residues to correct; safety of dual-strand windows).
- **Input:** DdCBE/TALED screens of mt-encoded OXPHOS subunits. **BE3D today: needs_readiness_feature** (mitochondrial code + dual-strand editor-profile + heteroplasmy), then OXPHOS complex structures map natively.
- **Why BE3D:** Maps edit consequences onto respiratory-chain complexes; distinguishes assembly vs catalytic residues.
- **Alternatives:** In-house structural/biophysics.
- **Decision criteria:** Organelle codon/complex support; data security.
- **Sophistication:** Very high. **Adoption: LOW-MED**. [Pretzel launch: biospace 2022]

---
**I16 — Ag-bio crop prime/base-editing company (extends I4)**
- **Role:** Crop trait engineering via base/prime editing. **Pairwise** (base-editing license to **Tropic** for banana/coffee), **Inari** (multiplex editing, $144M Series G, 10–20% yield targets), prime editing demonstrated in rice/wheat.
- **Goals:** Identify functionally critical residues to edit/avoid in a trait gene; transfer hotspots across cultivars/species.
- **Input:** Trait/functional screens in non-model plants. **BE3D today: WORKS TODAY** — custom FASTA+PDB, AlphaFold non-model coverage, built-in cross-species alignment.
- **Why BE3D:** Cross-species hotspot transfer + structure-guided edit choice; prime editor-profile would broaden reachable substitutions (readiness upside, not blocker).
- **Alternatives:** Bespoke agronomic pipelines.
- **Decision criteria:** Non-human proteins; cross-cultivar conservation.
- **Sophistication:** Med-High. **Adoption: LOW-MED** — few such structured screens but distinctive. [Pairwise/Tropic: pairwise.com; Inari $144M: agtechnavigator; ag genome-editing market ~$11.5B(2025)→$48B(2034): fortunebusinessinsights]

---
**I17 — Structural-genomics / foundation-model (techbio) team (LFC3D as labels)**
- **Role:** Builds variant-effect / protein foundation models (ESM, AlphaMissense lineage); wants structure-anchored functional labels.
- **Goals:** Use **LFC3D scores as 3D-structural training/eval labels** — fine-tuning PLMs on functional-screen data measurably improves variant-effect prediction; structure improves it further.
- **Input:** Public BE/PE screens + AlphaFold. **BE3D today: WORKS TODAY** — export LFC3D per-residue tables as labels/features (MIT license).
- **Why BE3D:** A ready, structure-aware aggregation with a null model → cleaner labels than raw noisy per-guide LFC.
- **Alternatives:** Reimplement neighborhood aggregation; raw screens.
- **Decision criteria:** Code quality, label provenance, license.
- **Sophistication:** Very high. **Adoption: MED** (extends I8 with a concrete label-export use). [AlphaMissense: researchgate 374043093; DMS fine-tuning PLMs: arXiv 2405.06729; structure improves VEP: ScienceDirect S0959440X25000417]

---
### NEW REGULATORY / GOVERNANCE PERSONA

---
**R1 — Biosecurity / biocontainment / dual-use reviewer (IBC, DURC/PEPP oversight)**
- **Role:** Institutional Biosafety Committee / dual-use-research-of-concern reviewer assessing editing proposals (WHO frameworks, NASEM, BWC context).
- **Goals:** Flag whether a proposed/edited variant could confer **gain-of-function** in a concerning protein; document functional-hotspot rationale.
- **Input:** Existing BE3D outputs (functional hotspots, +/- split scores). **BE3D today: WORKS TODAY** as an *interpretive/governance* consumer — the split +/- (GOF vs LOF) channel and hotspot maps give a structured, auditable read on which edits are activating; a GOF-direction flag would sharpen it (readiness upside).
- **Why BE3D:** Turns a screen into an auditable "which residues activate function" map for oversight (U17).
- **Alternatives:** Manual expert review; no structural layer.
- **Decision criteria:** Auditability; GOF/LOF direction clarity.
- **Sophistication:** Med-High. **Adoption: LOW** — novel angle, not a primary market, but a credibility/safety differentiator. [WHO frameworks: link.springer 10.1007/s11673-024-10358-8; NASEM: PMC12222282; DURC/PEPP/BWC: PMC12183502]

---

## 2. NEW USE CASES (extend Q2's U1–U10)

- **U11** Map organelle-editor (mtDNA DdCBE/TALED) screen onto mitochondrial protein / OXPHOS-complex structure.
- **U12** Prioritize allogeneic-CAR-T / multiplex edit sites — pick coding knockout positions that fully null each target with minimal collateral.
- **U13** Prime-editing **saturation-SNV → 3D** spatial aggregation (VUS resolution + PS3 spatial support).
- **U14** Single-cell **multi-phenotype hotspots** — per-phenotype LFC3D, then cross-phenotype consensus.
- **U15** **In-vivo therapeutic target-residue prioritization** — rank installable knockout/correction edits by structural consequence (Verve/Beam-style).
- **U16** Export **LFC3D as structural labels** for ML / foundation models.
- **U17** **Biosecurity / dual-use screen review** — flag GOF-conferring functional hotspots (auditable direction).

---

## 3. UPDATED PERSONA × USE-CASE FIT MATRIX (NEW personas × key + new use cases)

Columns: selected existing high-signal use cases (U1 prioritize validation, U3 interpret tiling driver
screen, U8 meta-aggregate, U10 pub/G2P figures) + the seven NEW use cases (U11–U17). H/M/L = fit strength.
`RF` marks a cell that **needs a readiness feature** to reach full strength today.

| Persona | U1 | U3 | U8 | U10 | U11 | U12 | U13 | U14 | U15 | U16 | U17 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| A12 Mito DdCBE/TALED | M | M | M | H | **H·RF** | L | L | L | M | L | L |
| A13 Prime saturation-SNV | H | H | H | H | L | L | **H·RF** | L | M | M | L |
| A14 Single-cell/Perturb-BE | H | M | H | M | L | M | M | **H·RF** | M | M | L |
| A15 Neuroscience BE | H | H | M | H | L | L | M | L | **H** | L | L |
| A16 Immuno-onc T-cell (acad.) | H | M | M | H | L | **H·RF** | L | M | H | L | M |
| A17 Microbiome/bacterial BE | H | M | H | M | L | L | M | L | M | M | M |
| I10 Hemoglobinopathy (Beam) | L | L | L | M | L | L | L | L | **L·RF** | L | L |
| I11 Cardiometabolic (Verve/Lilly) | H | M | H | M | L | L | M | L | **H** | M | M |
| I12 Liver/AATD (Beam-302) | H | H | M | M | L | L | M | L | **H** | M | M |
| I13 Allogeneic CAR-T co. | H | M | H | M | L | **H·RF** | L | M | H | L | M |
| I14 RNA-editing (ADAR) co. | M | M | M | M | L | L | M | L | **M·RF** | M | L |
| I15 Mito biotech (Pretzel) | M | M | M | M | **H·RF** | L | L | L | M | L | L |
| I16 Ag-bio crop editing | M | M | M | M | L | L | M | L | M | M | L |
| I17 Foundation-model/techbio | M | M | H | L | L | L | H | M | L | **H** | L |
| R1 Biosecurity reviewer | M | M | M | M | L | L | L | L | M | L | **H** |

(Existing A1–A11 / I1–I9 map primarily onto U1–U10 as in Q2; their fit on new use cases: A4/A5 → U13
(M–H), A8/I8 → U16 (H), I1/I2 → U15 (H). `RF` cells are the readiness-gated opportunities.)

---

## 4. INCREMENTAL ADDRESSABLE PATIENT POPULATION (emerging-tech personas)

BE3D itself is an *interpretation tool*, but its patient impact is levered through the disease areas the
new **therapeutic** personas serve. Cited disease burdens the emerging personas add on top of Q2's mostly
academic/target-discovery framing:

**Tier 1 — large cardiometabolic surface (I11, I14-cardio):**
- Heterozygous familial hypercholesterolemia ~**1 in 250** (~1.3M US; ~30M+ global). [CDC blogs.cdc.gov/genomics 2021/01/25]
- Elevated Lp(a) ~**20% of the population** (~1 in 5; >1.4B people globally). [AHA JAHA 124.040361]
- These sit atop tens of millions with ASCVD — the largest patient pool any new persona adds.

**Tier 2 — severe monogenic disease (I10, I12, I14-AATD, A12/I15):**
- Sickle cell disease: GBD 2021 ~**5.68M** cases globally (+67% since 2000; ~300k affected births/yr); **US ~100,000**. Severe thalassemia ~**400,000**. [PMC10390339; CDC]
- AATD severe PiZZ: ~**1 in 2,000–4,000** of European descent (Z allele ~4% of N. Europeans). [Orphanet J 10.1186/1750-1172-3-16]
- Mitochondrial disease: diagnosed ~**1 in 5,000–10,000** births; mtDNA-mutation carriers ~**1 in 200**; TALED could in principle correct **~43% of 90** known pathogenic mtDNA mutations. [Springer 10.1007/s00415-015-7884-3; sciencedaily 2022/04]

**Tier 3 — oncology / cell therapy (A16, I13):**
- T-ALL ≈ **20%** of ~**6,100** US ALL cases/yr (worse-prognosis relapsed disease is the initial indication); base-edited allogeneic CAR-T extends to broader hematologic cancers. [Nature bcj201753]

**Tier 4 — clinical-interpretation surface (A13, A5-adjacent, I17):**
- Saturation prime/base editing targets the **>1,000,000 VUS** in ClinVar — not a patient count but a
  variant-interpretation surface touching millions of carriers across all Mendelian genes. [Cell Genomics S2666-979X(25)00070-9]

**Tier 5 — rare neuro (A15):** prion disease is rare (~1–2 per million) but the same in-vivo BE approach
generalizes to **SOD1-ALS, APOE/Alzheimer's** (tens of millions) — high per-persona optionality. [Nat Med s41591-024-03466-w]

**Net.** The emerging-tech personas expand BE3D's relevance from a few-hundred-lab academic substrate to
the **interpretation layer behind therapies aimed at a tier-1 cardiometabolic pool in the tens-of-millions
to >1B (Lp(a)), a tier-2 severe-monogenic pool of ~6–7M+ (SCD + thal + AATD + mito), and a >1M-VUS clinical
surface.** The caveat (patient-ethos honest): BE3D touches these **indirectly** (target-residue selection,
edit-site prioritization, VUS spatial support) — its lever is *faster, less-wasted validation and safer edit
choice*, most directly for the **works-today** therapeutic personas (I11 cardiometabolic, I12 AATD) and
the **readiness-gated** frontier (mito, RNA, multiplex CAR-T, saturation-PE).

---

## 5. NEW-PERSONA JSON (for dashboard)

```json
[
  {"id":"A12","name":"Mitochondrial-disease / organelle-editing researcher (DdCBE/TALED)","sector":"academic","role":"mtDNA base-editing / OXPHOS structure-function lab","fit_score":3,"primary_use_case":"Map mtDNA DdCBE/TALED screen onto OXPHOS-complex structure (U11)","needs_readiness_feature":true,"adoption_likelihood":"Medium"},
  {"id":"A13","name":"Prime-editing / saturation-SNV variant-mapping lab","sector":"academic","role":"Saturation prime-editing VUS mapper","fit_score":4,"primary_use_case":"Saturation-SNV -> 3D spatial aggregation for VUS (U13)","needs_readiness_feature":true,"adoption_likelihood":"Medium-High"},
  {"id":"A14","name":"Single-cell / Perturb-BE multi-phenotype screener","sector":"academic","role":"scBE / Perturb-seq multi-phenotype lab","fit_score":3,"primary_use_case":"Per-phenotype 3D hotspots + cross-phenotype consensus (U14)","needs_readiness_feature":true,"adoption_likelihood":"Medium"},
  {"id":"A15","name":"Neuroscience / neurodegeneration in-vivo BE lab","sector":"academic","role":"AAV/LNP base editing for PRNP/SOD1/APOE","fit_score":4,"primary_use_case":"Prioritize/interpret neuro coding edit residues (U15)","needs_readiness_feature":false,"adoption_likelihood":"Medium-High"},
  {"id":"A16","name":"Immuno-oncology / T-cell engineering academic lab","sector":"academic","role":"Multiplex base-edited CAR-T (UCL/GOSH lineage)","fit_score":3,"primary_use_case":"Prioritize multiplex knockout edit sites (U12)","needs_readiness_feature":true,"adoption_likelihood":"Medium"},
  {"id":"A17","name":"Microbiome / bacterial base-editing lab","sector":"academic","role":"High-throughput bacterial BE (AMR/microbiome)","fit_score":3,"primary_use_case":"3D hotspots on bacterial enzymes from genome-scale BE screens","needs_readiness_feature":false,"adoption_likelihood":"Low-Medium"},
  {"id":"I10","name":"Hemoglobinopathy ex-vivo BE program (Beam BEAM-101)","sector":"industry","role":"SCD/beta-thal HbF-reactivation program","fit_score":2,"primary_use_case":"BCL11A protein-DNA interface (edit site is non-coding promoter)","needs_readiness_feature":true,"adoption_likelihood":"Low"},
  {"id":"I11","name":"Cardiometabolic in-vivo BE program (Verve/Lilly: PCSK9/ANGPTL3/LPA)","sector":"industry","role":"In-vivo liver base-editing lipid program","fit_score":4,"primary_use_case":"In-vivo therapeutic target-residue prioritization (U15)","needs_readiness_feature":false,"adoption_likelihood":"Medium"},
  {"id":"I12","name":"Liver / AATD in-vivo BE program (Beam BEAM-302, SERPINA1 E342K)","sector":"industry","role":"In-vivo SERPINA1 correction program","fit_score":4,"primary_use_case":"Structural interpretation of E342K / polymerization residues (U15)","needs_readiness_feature":false,"adoption_likelihood":"Medium"},
  {"id":"I13","name":"Allogeneic CAR-T / multiplex-BE cell-therapy company","sector":"industry","role":"Multiplex base-edited allogeneic cell therapy (Beam/Caribou/Prime)","fit_score":3,"primary_use_case":"Prioritize multiplex knockout edit sites (U12)","needs_readiness_feature":true,"adoption_likelihood":"Medium"},
  {"id":"I14","name":"RNA-editing (ADAR, A-to-I) therapeutics company","sector":"industry","role":"Programmable RNA editing (Wave/Korro/AIRNA/Shape)","fit_score":2,"primary_use_case":"Structural read of recoded residue consequence (U15)","needs_readiness_feature":true,"adoption_likelihood":"Low-Medium"},
  {"id":"I15","name":"Mitochondrial-medicine biotech (Pretzel Therapeutics)","sector":"industry","role":"Therapeutic mtDNA editing / heteroplasmy shifting","fit_score":3,"primary_use_case":"Map mt-protein edit consequences on OXPHOS complexes (U11)","needs_readiness_feature":true,"adoption_likelihood":"Low-Medium"},
  {"id":"I16","name":"Ag-bio crop prime/base-editing company","sector":"industry","role":"Crop trait engineering (Pairwise/Inari/Tropic)","fit_score":2,"primary_use_case":"Cross-species hotspot transfer + edit choice in trait genes","needs_readiness_feature":false,"adoption_likelihood":"Low-Medium"},
  {"id":"I17","name":"Structural-genomics / foundation-model (techbio) team","sector":"industry","role":"Variant-effect / protein foundation-model builder","fit_score":3,"primary_use_case":"Export LFC3D as structural labels for ML (U16)","needs_readiness_feature":false,"adoption_likelihood":"Medium"},
  {"id":"R1","name":"Biosecurity / biocontainment / dual-use reviewer","sector":"regulatory","role":"IBC / DURC-PEPP oversight reviewer","fit_score":2,"primary_use_case":"Flag GOF-conferring functional hotspots (U17)","needs_readiness_feature":false,"adoption_likelihood":"Low"}
]
```

---

## 6. SOURCES (new, this round)

- **Mitochondrial:** DdCBE broadinstitute.org/news/new-molecular-tool-precisely-edits-mitochondrial-dna ; TALED Cell 2022 sciencedirect.com/science/article/pii/S0092867422003890 ; sciencedaily.com/releases/2022/04/220425121101.htm ; DdCBE cryo-EM nature.com/articles/s41467-023-42359-3 ; Pretzel biospace.com/pretzel-therapeutics-launches ; burden link.springer.com/article/10.1007/s00415-015-7884-3
- **RNA editing:** Wave ir.wavelifesciences.com/news-releases/.../wave-life-sciences-announces-positive-update-restoraation-2 ; Korro biospace.com (KRRO-110 Q3) ; AIRNA genengnews.com/topics/genome-editing/airna-closes-155m-series-b-round ; market natlawreview.com/press-releases/rna-editing-therapies-market
- **CAR-T / multiplex:** BEAM-201 investors.beamtx.com/news-releases/.../beam-therapeutics-names-first-car-t-base-editing-development ; ASH2024 ashpublications.org/blood/article/144/Supplement%201/4838 ; Alyssa ucl.ac.uk/news/2022/dec/world-first-use-base-edited-car-t-cells ; NCT05397184 ; translocations PMC10123993 ; Caribou cgtlive.com CB-010
- **Hemoglobinopathy:** BEAM-101 investors.beamtx.com (BEACON EHA2025) ; Casgevy ir.crisprtx.com FDA approval ; GBD SCD PMC10390339
- **Cardiometabolic:** Lilly/Verve cnbc.com/2025/06/17 ; VERVE-102 Heart-2 vervetx.gcs-web.com ; PCSK9 LOF PMC5729040 ; FH CDC blogs.cdc.gov/genomics 2021/01/25 ; Lp(a) AHA JAHA 124.040361
- **AATD:** BEAM-302 globenewswire.com 2025/03/10 ; E342K PMC4957051 ; Prime Medicine investors.primemedicine.com 2025/03 ; prevalence Orphanet 10.1186/1750-1172-3-16
- **Prime / saturation-SNV:** nature.com/articles/s41587-020-0455-x ; Cell Genomics S2666-979X(25)00070-9 ; PEER-seq PMC12008803 ; Nat Biotech s41587-024-02465-z
- **Single-cell BE:** Nat Methods s41592-024-02256-z ; Cell Systems S2405-4712(24)00366-1
- **Ag-bio:** pairwise.com (Tropic license) ; agtechnavigator.com (Inari $144M) ; fortunebusinessinsights.com/genome-editing-market
- **Microbiome/bacterial:** NAR academic.oup.com/nar/article/52/7/4079 ; PMC11039988 ; PNAS 2206744119
- **Neuroscience:** Nat Med nature.com/articles/s41591-024-03466-w ; broadinstitute.org/news/gene-editing-extends-lifespan-mouse-model-prion-disease
- **Foundation models:** AlphaMissense researchgate 374043093 ; PLM DMS fine-tuning arXiv 2405.06729 ; structure-VEP ScienceDirect S0959440X25000417
- **Biosecurity:** WHO link.springer.com 10.1007/s11673-024-10358-8 ; NASEM PMC12222282 ; DURC/PEPP/BWC PMC12183502 ; carnegieendowment.org 2024/10
</content>
</invoke>
