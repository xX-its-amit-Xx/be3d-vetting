# P2 — Scout: real public BE tiling / mutational-scanning datasets for BE3D

**Date:** 2026-08-23
**Author:** Scout Queen (delegated verification by 4 sub-agents; all download URLs curl-probed in-session)
**Goal:** 5–7 real, PUBLIC base-editor tiling / mutational-scanning screens (EXCLUDING shipped
examples KBTBD4, HDAC1, MORC2) that we can download and feed to BE3D, each paired with a
*genuinely open* 3D-hotspot question (not the source paper's headline).

BE3D input reminder (from the shipped example TSVs, e.g. `YeoKBTBD4…ABE-Screen.tsv`,
`MORC2/…TadA…Input.tsv`): columns `ID, Gene, Mutation_list, Mutation_type, sgRNA_score`.
`Mutation_list` = per-guide AA edits, `;`-joined (e.g. `M1V` or `A25T;G24K;`);
`Mutation_type` = matching categories (`Missense;Silent;`, or `Missense`, `Nonsense`, `Splice`,
`No Mutation`, `UTR`, `Flanking`); `sgRNA_score` = the numeric per-guide LFC/effect. So a usable
dataset needs, per guide: an AA-level edit + a category + a numeric score.

## Access notes learned during verification (important)
- **`ars.els-cdn.com` (Elsevier CDN)** and **`static-content.springer.com` (Springer)** both serve
  supplementary `.xlsx` openly to bare `curl` (HTTP 200, verified). **Use these, not** `cell.com`
  or `nature.com` article pages (Cloudflare / auth 303) **nor** PMC `.../bin/mmc*` links
  (JS proof-of-work → ~21 KB placeholder). This resolved an early false "paywalled" read.
- **MaveDB API** is fully curl-downloadable: `https://api.mavedb.org/api/v1/score-sets/<URN>/scores`
  → CSV with `accession,hgvs_nt,hgvs_splice,hgvs_pro,score,…` (verified on several URNs).
- **Negative result (useful):** no public *base-editor* tiling scan exists for **BTK/ibrutinib,
  ABL1/imatinib, or ALK** — only ENU/clinical substitution maps. Those clinically important
  pockets remain BE-data gaps.

---

## Ranked candidate table

Feasibility = data-accessibility (1–5) × structural-tractability (1–5) × novelty (1–5); max 125.

| # | Dataset (paper, year) | Data URL(s) — **verified curl 200** | Format & key columns | Gene(s) → UniProt | Editor / assay | Structure (AF / PDB / complex) | Open 3D gap → NEW BE3D output | Feasibility (A×S×N) | Reformat risk |
|---|---|---|---|---|---|---|---|---|---|
| **1** | **Coelho/Dincer et al. 2024, Nat Genet** "Base editing screens define the genetic landscape of cancer drug resistance" (Sanger) | `https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-024-01948-8/MediaObjects/41588_2024_1948_MOESM4_ESM.xlsx` (62 MB) · MaveDB subsets `urn:mavedb:00001204-a-1..6` | xlsx, sheet **`ST2 BE z-scores`** (45,129 rows): `Protein_Change` (AA), `variant_classification`/`most_severe_consequence` (type), per-drug **`L2FC_<cell>_<drug>_plasmid_average_zscore`** (score), + built-in `swissprot`,`Domain`,`is_disruptive_interface_*` | EGFR **P00533**, KRAS **P01116**, BRAF **P15056**, MAP2K1 **Q02750**, MAP2K2 P36507, AKT1 P31749, PIK3CA **P42336**, PARP1 **P09874**, PARP2 Q9UGN5 (11 genes) | ABE **and** CBE (Cas9-NG); proliferation under 10 drugs × 4 cell lines | EGFR-osimertinib 4ZAU / triple-mut+osi 6LUD; KRAS-G12C sotorasib 6OIM, adagrasib 6UT0; MEK1 allosteric 4MNE (no trametinib co-crystal); PARP1-olaparib 5DS3 | Paper flags per-variant interface disruption but never *spatially clusters* resistance scores on structure. **New:** a per-drug 3D resistance atlas separating **orthosteric (pocket-lining) vs allosteric/PPI-interface** resistance clusters; e.g. do sotorasib vs adagrasib resistance form distinct sub-clusters *within* the switch-II pocket? | 5×5×4 = **100** | Low (turnkey; 3-letter→1-letter HGVS conversion only) |
| **2** | **Hanna et al. 2021, Cell** "Massively parallel assessment of human variants with base editor screens" (Doench/GPP, Broad) | `https://ars.els-cdn.com/content/image/1-s2.0-S009286742100012X-mmc2.xlsx` (BRCA1/2, has Z-score) · `…-mmc4.xlsx` (MCL1/BCL2L1) · `…-mmc5.xlsx` (drug-target tiling) | xlsx "Library annotation" sheet: **`Amino acid edits`** (`Gly1015Lys;`), **`Mutation category`** (`Missense;`/`Silent;`/…). Data sheets = **raw counts** → compute LFC=log2(cond/pDNA); **mmc2** BRCA1/2 sheet ships a ready **`Z-score`** | MCL1 **Q07820**, BCL2L1 **Q07817**, BRCA2 **P51587**, BRCA1 P38398, PARP1 P09874 | CBE (BE3.9max/BE4max); dropout + drug arms (PARPi, cisplatin, S63845, A1331852) | **MCL1 BH3-peptide complexes 3PK1 (BAK), 3MK8 (BID)** + inhibitor 6QB4/5FDR; Bcl-xL 1BXL(BAK)/2BZW(BAD); BRCA2-DSS1-ssDNA 1MJE/1N0W | **New:** on the *same* MCL1 domain, does the dropout-essentiality hotspot map to the **BH3-groove PPI surface** (contacts to BAK/BID) while inhibitor-resistance maps to the **orthosteric drug pocket** — can BE3D separate a PPI-dependency interface from a drug-contact hotspot? For BRCA2, DSS1-interface vs ssDNA-OB-fold partition. | 3×5×5 = **75** | Medium (compute LFC + join on `sgRNA sequence`; BRCA1/2 subset already scored) |
| **3** | **Cuella-Martin et al. 2021, Cell** "Functional interrogation of DNA damage response variants with base editing screens" (Ciccia lab) | `https://ars.els-cdn.com/content/image/1-s2.0-S0092867421000842-mmc5.xlsx` (11 MB, Table S5, **LFC + AA + category**) · raw counts `…-mmc2.xlsx` | xlsx (sheets per cell line): **`AAChg`** (HGVS p., e.g. `p.T842I`), **`Function`** (`missense`/`synonymous`/nonsense/splice), per-condition **`T18_UNT_LFC` / `T18_CISP_LFC` / `T18_OLAP_LFC` / `T18_DOX_LFC` / `T18_CPT_LFC`** (+ p/FDR) | BRCA1 **P38398**, TP53 **P04637**, ATM Q13315, BARD1, CHEK2 O96017 (~86 DDR genes) | CBE (BE3-FNLS); competitive dropout ± cisplatin/olaparib/doxorubicin/camptothecin; MCF10A/MCF7/HAP1 | BRCA1-BARD1 RING **1JM7**; BRCA1-BRCT-pSer **1T29**; TP53-DNA 2AC0/3KMD tetramer; ATM 7NI4 | **New:** map *drug-specific* LFC (OLAP vs CISP) of BRCA1 RING missense onto the 1JM7 heterodimer — do condition-specific hotspots cluster on the **BARD1-dimerization interface vs the E2/UbcH5 ligase surface**, resolving heterodimerization- from ligase-dependent vulnerability? | 4×4×4 = **64** | Medium (use **mmc5/S5**, which carries `AAChg`; the processed mmc4 lacks residue-level AA — confirm S5 column on load) |
| **4** | **Lue et al. 2023, Nat Chem Biol** "Base editor scanning charts the DNMT3A activity landscape" (**Liau lab**) | `https://static-content.springer.com/esm/art%3A10.1038%2Fs41589-022-01167-4/MediaObjects/41589_2022_1167_MOESM3_ESM.xlsx` (Suppl. Data 1–7 in one book) | xlsx **Suppl. Data 3 ("Summary")**: **`Mutation_list`** (`E664K, V665M`), **`Mut_type`** (Missense/Silent/Nonsense/Splice/Non-exon), **`sgRNA_score_d9_citrine`** (numeric effect), `Hit_status` | DNMT3A **Q9Y6K1** (partner DNMT3L Q9UJW3) | CBE (BE3.9max/evoAPOBEC1-BE4max); **CpG-methylation→citrine reporter FACS** (activity, not dropout); MOLM13 | DNMT3A-DNMT3L-H3 active **4U7T**; autoinhibited 4U7P; DNMT3A-3L-DNA tetramer | Paper's PWES already finds a general interdomain cluster. **New:** partition LOF density across the *distinct* 4U7T surfaces — DNMT3A–DNMT3L (FF) interface vs homotetramer (RD) interface vs DNA face vs ADD–MTase autoinhibitory contact — to name the dominant **PPI interface** (not the catalytic pocket) vulnerability, and test CHIP-hotspot co-localization. | 5×4×3 = **60** | Low (one-table, scores shipped) |
| **5** | **Sangree et al. 2022, Nat Commun** "Benchmarking of SpCas9 variants enables deeper base editor screens of BRCA1 and BCL2" (Doench/GPP) | base `https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-022-28884-7/MediaObjects/` → `…MOESM9` (BRCA1 CBE), `…MOESM10` (BRCA1 ABE), `…MOESM11` (**BCL2 venetoclax**), `…MOESM13` (Source Data, z/LFC) | xlsx annotation sheet: **`Amino acid edits`**, **`Mutation category`**; counts per rep/arm → compute LFC, or pull z from MOESM13 | BCL2 **P10415**, BRCA1 **P38398** | CBE + ABE (WT/NG/SpG Cas9); BRCA1 dropout+cisplatin; **BCL2 venetoclax-resistance enrichment** (MOLM13) | BCL2-venetoclax **6O0K**; resistance mutants 6O0P (G101A), 6O0M (F104L) | **New:** do venetoclax-resistance residues form a *contiguous 3D shell* lining the BH3 groove **beyond** the known G101/D103/F104 rim (mapped on 6O0K) — an allosteric vs orthosteric resistance surface? | 3×4×4 = **48** | Medium (LFC from counts / Source Data; join on guide) |
| **6** | **Sánchez-Rivera et al. 2022, Nat Biotech** "Base editing sensor libraries … cancer-associated SNVs" | annotation `https://static-content.springer.com/esm/art%3A10.1038%2Fs41587-021-01172-3/MediaObjects/41587_2021_1172_MOESM4_ESM.xlsx`; scores `…MOESM9_ESM.xlsx` | xlsx: MOESM4 `Table_S2C` → **`HGVSp_Short`** (`p.R206H`), **`Variant_Classification`**; MOESM9 `Table_S7` → `sgrna` (`Trp53_E258K_4450`), **`LFC`**, `FDR` | mouse Trp53 **P02340** (human TP53 P04637) | 9 CBEs; proliferation/dropout in pancreatic cells | TP53-DNA **1TUP**/2AC0; tetramer 3KMD | **New:** do proliferation-scored Trp53 missense hotspots extend *beyond* canonical R175/R248/R273 across the DNA-contact surface + Zn-structural core (1TUP), yielding a BE-specific spatial LOF map? | 3×4×3 = **36** | Medium (AA + LFC in separate tables; join on parseable `sgrna` string; mouse numbering) |
| **7** | **Belli/Platt et al. 2024, Nat Biotech** "Multimodal scanning of genetic variants with base and prime editing" (EGFR) | `https://static-content.springer.com/esm/art%3A10.1038%2Fs41587-024-02439-1/MediaObjects/41587_2024_2439_MOESM5_ESM.xlsx` (screen scores) + `…MOESM3_ESM.xlsx` (annotation); GitHub `plattlab/multimodal_genetic_variants` | xlsx: scores sheet `Osimertinib_resistance_ABE8e/…` → `sgrna`,**`LFC`**,`FDR`; annotation sheet → **`Predicted amino acid change`**,**`Mutation category`** — **join on `sgRNA sequence`** | EGFR **P00533** | ABE8e + CBE(BE3.9max) + PE, full EGFR CDS; **gefitinib/osimertinib** selection | EGFR-osimertinib 4ZAU; L858R/T790M/C797S+osi 6LUD | **New:** do "blocks-the-drug" (C797S covalent-site) vs "reactivates-the-kinase" variants segregate into different EGFR compartments (ATP-pocket rim vs αC-helix/dimer interface); are 1st- vs 3rd-gen-TKI resistance hotspots spatially offset? | 3×4×3 = **36** | Medium (cross-sheet sgRNA join) |

### Bonus / fallback (NOT base editors — per-variant AA score + structure, higher reformat risk)
SGE/DMS on MaveDB (uniform curl CSV; **category must be derived** from `hgvs_pro`): BRCA1 SGE
`urn:mavedb:00000097-0-2` (P38398, 1JM7/1T15); MSH2 SGE `urn:mavedb:00000050-a-1` (P43246,
2O8B dimer-DNA); DDX3X SGE `urn:mavedb:00000658-0-1`; **LDLR DMS `urn:mavedb:00001269-a-1`**
(P01130; **3GCW LDLR–PCSK9 interface** = novel PPI-hotspot question). Also the published TP53
base-editor set on MaveDB `urn:mavedb:00001245-a-1` (ABE8e) / `-a-2` (CBE), score=LFC — real BE,
but multi-residue `hgvs_pro` (`p.[Met133Thr;Phe134Leu]`) raises the mapping question of one LFC →
several residues (reformat Medium-High).

---

## Top-3 run recipes

### #1 — Coelho/Dincer 2024 (drug-resistance, multi-target) — **the turnkey pick**
- **File:** `41588_2024_1948_MOESM4_ESM.xlsx`, sheet `ST2 BE z-scores`. Filter one gene at a time
  (e.g. `swissprot == "P01116"` for KRAS).
- **Column mapping:**
  - `Mutation_list` ← `Protein_Change` (convert 3-letter HGVS `p.Gly12Cys` → `G12C`; drop `p.`).
  - `Mutation_type` ← `variant_classification` (`Missense_Mutation`→`Missense`, `Silent`→`Silent`,
    `Nonsense_Mutation`→`Nonsense`, `Splice_Site`→`Splice`).
  - `sgRNA_score` ← the drug/cell z-score column, e.g. `L2FC_H23_Sotor_plasmid_average_zscore`
    (sotorasib) **and** `…_Adag_…` (adagrasib) as two separate BE3D screens for the same gene.
  - `Gene` ← `Hugo_Symbol`; UniProt ← `swissprot`.
- **Structure:** KRAS **P01116**, PDB **6OIM** (G12C+sotorasib) and **6UT0** (adagrasib). **Monomer**
  for the pocket question; optionally **PPI/complex** with RAF-RBD to test the RAS–RAF interface.
- **NEW question:** run sotorasib and adagrasib arms separately through LFC3D + 3D clustering; test
  whether their resistance hotspots form **spatially separable sub-clusters within the switch-II
  pocket** (steric-contact vs conformational-latch residues) — a drug-specific 3D resistance
  signature invisible to the pooled hit list. Repeat for EGFR (gefitinib vs osimertinib) and the
  MAPK axis (BRAF/MAP2K1) to build a "pocket vs interface" resistance atlas.

### #2 — Hanna 2021, MCL1 / BCL2L1 (apoptosis, PPI + drug) — **highest novelty**
- **File:** `1-s2.0-S009286742100012X-mmc4.xlsx` (MCL1/BCL2L1 tiling), "Library annotation" sheet
  for edits, data sheet for counts.
- **Column mapping:**
  - `Mutation_list` ← `Amino acid edits` (already `;`-joined 3-letter HGVS → convert to 1-letter).
  - `Mutation_type` ← `Mutation category`.
  - `sgRNA_score` ← compute per guide: `log2(mean(dropout reps)/pDNA)` (and separately the
    S63845/A1331852 drug arms), joining to the annotation on `sgRNA sequence`.
  - `Gene` ← MCL1; UniProt **Q07820**.
- **Structure / mode:** **PPI/complex mode** — MCL1 chain + BAK-BH3 peptide, PDB **3PK1** (and BID
  3MK8), plus inhibitor complex **6QB4** for the drug arm. Chain-aware LFC3D (as in the shipped
  KBTBD4–HDAC1 PPI example) gives interface hotspots.
- **NEW question:** does the **dropout (essentiality)** hotspot localize to the **BH3-groove PPI
  surface** while the **inhibitor-resistance** hotspot localizes to the **orthosteric pocket**? i.e.
  can BE3D dissociate a protein–protein dependency interface from a drug-contact surface on one
  domain — directly actionable for BH3-mimetic design.

### #3 — Cuella-Martin 2021, BRCA1 (DDR, drug-specific interface)
- **File:** `1-s2.0-S0092867421000842-mmc5.xlsx` (Table S5); pick a BRCA1 sheet. **Confirm the
  `AAChg` column is present on load** (this file carries residue-level AA; the processed mmc4 does
  not).
- **Column mapping:**
  - `Mutation_list` ← `AAChg` (HGVS p.; take the driver edit, convert to 1-letter).
  - `Mutation_type` ← `Function` (`missense`→`Missense`, `synonymous`→`Silent`, etc.).
  - `sgRNA_score` ← `T18_OLAP_LFC` (olaparib) and `T18_CISP_LFC` (cisplatin) as two screens; both
    already computed (with `_p-value`/`_FDR`).
  - `Gene` ← BRCA1; UniProt **P38398**.
- **Structure / mode:** **PPI/complex mode** on the BRCA1–BARD1 RING heterodimer **1JM7** (and BRCT
  cleft 1T29 as a second run).
- **NEW question:** does the *drug-specific* LFC (OLAP vs CISP) of BRCA1 RING missense partition
  between the **BARD1-dimerization interface** and the **E2/UbcH5 ubiquitin-ligase surface** — a 3D
  hotspot that separates heterodimerization-dependent from ligase-activity-dependent vulnerability?

---

## Shortlist (JSON)

```json
[
  {"rank":1,"name":"Coelho/Dincer 2024 Nat Genet — cancer drug-resistance BE screens","gene":"KRAS/EGFR/BRAF/MAP2K1/PIK3CA/PARP1 (11)","uniprot":"P01116 (KRAS); P00533; P15056; Q02750; P42336; P09874","url":"https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-024-01948-8/MediaObjects/41588_2024_1948_MOESM4_ESM.xlsx","editor":"ABE+CBE (Cas9-NG); drug-selection proliferation","open_gap":"resistance scores never spatially clustered on structure; orthosteric vs allosteric/PPI-interface resistance undefined","new_output":"per-drug 3D resistance atlas separating drug-pocket-lining hotspots from allosteric/interface clusters (e.g. sotorasib vs adagrasib sub-clusters in switch-II)","feasibility":100},
  {"rank":2,"name":"Hanna 2021 Cell — MCL1/BCL2L1 CBE tiling","gene":"MCL1 (also BCL2L1, BRCA2)","uniprot":"Q07820","url":"https://ars.els-cdn.com/content/image/1-s2.0-S009286742100012X-mmc4.xlsx","editor":"CBE (BE3.9max); dropout + BH3-mimetic arms","open_gap":"cannot tell if essentiality hotspot is the BH3-groove PPI surface vs the drug pocket","new_output":"BE3D PPI-mode separation of the BAK/BID BH3-groove dependency interface (3PK1) from the orthosteric inhibitor pocket (6QB4) on MCL1","feasibility":75},
  {"rank":3,"name":"Cuella-Martin 2021 Cell — DDR BE screens (BRCA1)","gene":"BRCA1 (also TP53, ATM; 86 DDR genes)","uniprot":"P38398","url":"https://ars.els-cdn.com/content/image/1-s2.0-S0092867421000842-mmc5.xlsx","editor":"CBE (BE3-FNLS); dropout ± olaparib/cisplatin/dox/CPT","open_gap":"drug-specific missense LFC not mapped to BRCA1-BARD1 heterodimer; dimerization vs ligase surface unresolved","new_output":"OLAP-vs-CISP 3D hotspot partition on 1JM7 separating BARD1-dimer interface from E2/ligase surface","feasibility":64},
  {"rank":4,"name":"Lue 2023 Nat Chem Biol — DNMT3A base editor scanning (Liau)","gene":"DNMT3A","uniprot":"Q9Y6K1","url":"https://static-content.springer.com/esm/art%3A10.1038%2Fs41589-022-01167-4/MediaObjects/41589_2022_1167_MOESM3_ESM.xlsx","editor":"CBE; CpG-methylation citrine-reporter activity (FACS)","open_gap":"PWES found a general interdomain cluster but did not attribute LOF to specific PPI interfaces","new_output":"partition LOF density across DNMT3L (FF) vs tetramer (RD) vs DNA-face vs ADD-autoinhibitory surfaces on 4U7T; name the dominant interface vulnerability + test CHIP hotspots","feasibility":60},
  {"rank":5,"name":"Sangree 2022 Nat Commun — BRCA1/BCL2 deeper BE screens","gene":"BCL2 (also BRCA1)","uniprot":"P10415","url":"https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-022-28884-7/MediaObjects/41467_2022_28884_MOESM11_ESM.xlsx","editor":"CBE+ABE (NG/SpG); venetoclax-resistance enrichment","open_gap":"venetoclax-resistance residues beyond the G101/D103/F104 rim not spatially organized","new_output":"contiguous 3D resistance shell lining the BH3 groove on 6O0K (orthosteric vs allosteric)","feasibility":48},
  {"rank":6,"name":"Sanchez-Rivera 2022 Nat Biotech — CBE sensor libraries (Trp53)","gene":"TP53/Trp53","uniprot":"P02340 (mouse); P04637 (human)","url":"https://static-content.springer.com/esm/art%3A10.1038%2Fs41587-021-01172-3/MediaObjects/41587_2021_1172_MOESM9_ESM.xlsx","editor":"9 CBEs; proliferation/dropout","open_gap":"BE proliferation LFC not spatially mapped beyond canonical hotspot codons","new_output":"BE-specific spatial LOF map across the p53 DNA-contact surface + Zn-structural core (1TUP)","feasibility":36},
  {"rank":7,"name":"Belli/Platt 2024 Nat Biotech — multimodal EGFR scan","gene":"EGFR","uniprot":"P00533","url":"https://static-content.springer.com/esm/art%3A10.1038%2Fs41587-024-02439-1/MediaObjects/41587_2024_2439_MOESM5_ESM.xlsx","editor":"ABE8e+CBE+PE; gefitinib/osimertinib selection","open_gap":"resistance vs reactivating variants not spatially separated on the kinase","new_output":"3D split of drug-blocking (C797S pocket) vs kinase-reactivating (aC-helix/dimer interface) hotspots; 1st- vs 3rd-gen TKI offset","feasibility":36}
]
```

## Caveats / reformatting-risk summary
- **Turnkey (score shipped in one table):** Coelho `ST2`, DNMT3A `Suppl Data 3`, Cuella-Martin
  `mmc5/S5` (verify `AAChg` present).
- **Compute LFC from raw counts + join on `sgRNA sequence`:** Hanna 2021, Sangree 2022.
- **AA and LFC in separate sheets (string-parse guide IDs to join):** Sánchez-Rivera, Belli/Platt,
  Kim H.H. 2022 (Nat Biotech; 82 MB / 29 sheets, thin per gene — a weaker BE3D fit, listed only as
  context).
- **HGVS conversion:** most BE files use 3-letter HGVS (`p.Gly12Cys`) — convert to BE3D's 1-letter
  `Mutation_list` style (`G12C`); MaveDB sets additionally need **category derived** from `hgvs_pro`.
- **Numbering:** confirm mouse Trp53 (Sánchez-Rivera) and any MaveDB `offset` map to the canonical
  UniProt/PDB residues before overlay.
- **Multi-residue guides:** several BE guides edit >1 codon (`p.[A;B]`); BE3D's `Mutation_list`
  already supports `;`-joined edits per guide, so keep all edits rather than collapsing.
