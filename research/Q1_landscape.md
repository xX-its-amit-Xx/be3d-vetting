# BE3D — Competitive Landscape / Field Meta-Search (Queen 1)

**Prepared:** 2026-08-23
**Subject:** BE3D (`broadinstitute/BE3D`; pip `beclust3d` v1.0.0; MIT). A Python package that takes a *completed* base-editor (BE) tiling-mutagenesis functional screen, maps per-residue log-fold-change (LFC) onto a 3D protein structure, computes a spatially-smoothed **LFC3D** score with randomization-based statistics, and **agglomeratively clusters** significant residues into 3D structural **hotspots** (single-screen and meta across screens), exporting to the Genomics 2 Proteins (G2P) Portal.

---

## 0. Executive orientation — where BE3D actually sits

BE3D is a **downstream interpretation layer**, not a screen-caller, not a guide designer, not a pathogenicity predictor, and not an outcome/efficiency predictor. Its pipeline consumes an already-scored screen (a TSV of mutation → LFC) and produces **3D structural hotspots**. That places it at the intersection of two mature-but-separate worlds:

1. **CRISPR / base-editing tiling-screen analysis** (sequence/linear-coordinate world), and
2. **3D spatial mutation clustering / hotspot detection** (historically a *cancer somatic-mutation* world).

BE3D's distinctive move is to **bridge** them: it brings the 3D-spatial-clustering paradigm (born in cancer genomics, where the input signal is *mutation recurrence/frequency*) to the *functional-screen* signal (LFC from a designed tiling library) — and does so with base-editing-specific handling, cross-species/isoform meta-aggregation, and a direct handoff to an interactive portal.

**Most important single finding of this search:** BE3D is *not* alone in its exact niche. A direct competitor, **ProTiler-Mut** (He et al., *Cell Systems*, 2026), was published essentially concurrently and does the same high-level thing — map CRISPR **tiling mutagenesis** screens onto 3D structure and call **hotspot substructures**. The two differ in method and emphasis (see §3). This is the single most consequential item for vetting BE3D's novelty.

---

## 1. The field, defined precisely, and adjacent fields BE3D touches

### 1a. Base-editing / CRISPR screen ANALYSIS & tiling-mutagenesis interpretation

This is BE3D's home turf. These tools convert raw sgRNA read counts (or already-computed guide scores) into per-guide / per-gene / per-region signal. Most operate in **linear (sequence) coordinates** and do NOT use 3D structure.

| Tool | What it does | Niche | 3D? | Gap BE3D fills vs it |
|---|---|---|---|---|
| **MAGeCK / MAGeCK-VISPR** | RRA and MLE gene-level essentiality calling from pooled CRISPR-KO screen counts; VISPR adds QC + web viz. The de-facto standard first-pass caller. | Genome-wide dropout/enrichment; gene-level hit calling. | No | MAGeCK gives you the LFC/gene call; BE3D consumes that residue-level signal and asks *where in 3D space* it concentrates. Purely upstream. |
| **BAGEL2** | Bayesian gene essentiality (Hart lab); log Bayes factor per gene. | Essentiality classification, fold-change modeling. | No | Same as above — upstream gene-level, no residue/structure axis. |
| **drugZ** | Detects synergistic/suppressor drug–gene interactions in CRISPR screens (Hart lab). | Chemogenomic modifier screens. | No | Different question (drug interaction), no structural mapping. |
| **CRISPRO** (Schoonenberg et al., Genome Biol 2018) | Maps guide functional scores to genome/transcript/protein coordinates **and onto structures**; produces linear + 3D genotype-phenotype maps; LOESS smoothing along sequence. | The original "guide-score → protein/structure" mapper. | Partial (renders on structure; not a spatial-neighborhood statistic) | CRISPRO visualizes signal on structure but does **not** compute a spatial-neighborhood aggregated score with a null model, nor cluster residues into statistically significant 3D hotspots, nor meta-aggregate across screens. BE3D adds the LFC3D statistic + clustering + meta. |
| **ProTiler** (He et al., Nat Commun 2019) | Detects "CRISPR-KO hyper-sensitive (CKHS)" **regions** along the protein via a two-step denoise + region-calling on tiling-sgRNA dropout; can predict CKHS from conservation/domains/PTMs. | Linear functional-domain discovery from Cas9 tiling KO. | No (1D regions) | Works in 1D sequence space; BE3D works in 3D. ProTiler finds contiguous linear stretches; BE3D finds spatially-proximal residues that are *distant in sequence*. |
| **CRISPR-SURF** (Hsu et al., Nat Methods 2019) | Deconvolves tiling screen signal into a stepwise function to call regulatory elements / functional stretches (works for non-coding & coding tiling). | Deconvolution of tiling signal into 1D intervals. | No | 1D deconvolution only; no structure, no clustering in 3D, no cross-screen meta. |
| **BEAN** (`pinellolab/crispr-bean`; Ryu et al., Nat Genet 2024) | Bayesian network that **calls** base-editing variant effect sizes, normalizing by per-guide editing outcome (reporter) and chromatin accessibility; can share info among neighboring guides in dense tiling. | State-of-the-art *upstream* BE screen variant-effect **quantification**. | No | This is the **best current upstream partner**, not a competitor: BEAN produces cleaner per-variant/per-residue effect sizes; BE3D could take BEAN output as its LFC input and add the 3D-hotspot layer BEAN lacks. |
| **ScreenPro / ScreenPro2** (Nuñez/Weissman ecosystem) | Flexible Python toolkit for processing pooled screens (incl. CRISPRi/a, Perturb-seq-adjacent). | General screen processing/QC. | No | Upstream processing; no structural axis. |
| **be-scan / BEscreen / CRISPR-BEasy** | Base-editing library **design** (guide tiling design, variant enumeration). | Library design *before* the experiment. | No | Entirely upstream (design phase); orthogonal to BE3D. |
| **CRISPRcleanR** | Corrects copy-number bias in CRISPR-KO fold-changes. | Bias correction of screen LFC. | No | Preprocessing/QC; feeds cleaner LFC into any downstream tool including BE3D. |

**Field gap in 1a:** Every mature tool here is either (i) an upstream *caller* of per-guide/per-gene/per-variant signal, or (ii) a *1D* (linear sequence) interpreter. Before BE3D/ProTiler-Mut, there was **no standard, statistically-grounded way to aggregate functional-screen signal over a 3D spatial neighborhood and cluster it into structural hotspots** — you either eyeballed it on a structure viewer (CRISPRO) or stayed in 1D (ProTiler, SURF).

### 1b. BE OUTCOME / EFFICIENCY prediction (adjacent — NOT the same field)

These predict, *from sequence, before/without a screen*, what edit a base editor will make and how efficiently. They are strictly **upstream** of BE3D and answer a different question (what will the editor do?) rather than BE3D's (what does the observed phenotype mean structurally?).

| Tool | What it does | Relationship to BE3D |
|---|---|---|
| **BE-Hive** (Arbab et al., Cell 2020) | ML prediction of CBE/ABE editing outcomes (bystander patterns) + efficiency from target sequence. | Upstream/orthogonal — helps *design/interpret which* edits a guide makes; BE3D assumes edits are already annotated. |
| **BE-FF** | Predicts synonymous-correction outcomes to expand BE targeting scope. | Upstream design aid. |
| **BEdeepon (ABEdeepon/CBEdeepon)** | Deep-learning prediction of BE efficiency + outcomes. | Upstream. |
| **DeepBE / DeepBaseEditor** | CNN prediction of BE outcomes & efficiency. | Upstream. |
| **CRISPRitz** | Enumerates off-targets / variant-aware guide search. | Upstream, off-target-focused. |

**Gap:** None to fill — BE3D deliberately does not do this. The value is in being clear that BE3D is *not* competing here (a common reviewer confusion).

### 1c. 3D spatial mutation-CLUSTERING / hotspot detection (the CLOSEST analogues)

This is the paradigm BE3D borrows. **Crucially, almost all of these take *somatic mutation frequency/recurrence* as input** (count how often a residue is mutated across tumors, test whether mutations cluster in 3D beyond chance). BE3D differs fundamentally by taking a **continuous functional-screen signal (LFC)** as input rather than discrete mutation counts.

| Tool | What it does | Input signal | Gap BE3D fills / how it differs |
|---|---|---|---|
| **HotMAPS** (Tokheim et al., Cancer Res 2016; KarchinLab) | Density-based detection of missense hotspot *regions* in 3D structures; permutation null; residue-level, variable region size. | Somatic mutation counts | Same *machinery idea* (spatial density + permutation), but built for tumor mutation counts, not screen LFC; no BE handling, no cross-species meta, no screen-QA. |
| **HotSpot3D** (Niu et al., Nat Genet 2016; web server Bioinformatics 2020) | Graph-based hierarchical clustering on inter-residue atomic distances; intra/inter-molecular + drug–mutation clusters. | Somatic mutation counts | Closest *clustering* analogue. BE3D likewise clusters proximal significant residues but its inputs are LFC3D significance calls from a designed screen (with a randomization null tuned to screen structure), and it adds meta-aggregation + G2P export. |
| **CLUMPS** (Kamburov et al., PNAS 2015) | Tests whether mutations cluster in 3D more than expected (WAP statistic); pan-cancer driver discovery. | Somatic mutation frequency | Statistic assumes recurrence-weighted mutations; not designed for per-residue continuous screen scores. |
| **Mutation3D** (Meyer et al., Hum Mutat 2016) | Atomic clustering of coding variants across the structural proteome to nominate driver genes; web server. | Somatic/coding variants | Driver-gene nomination, not functional-screen hotspot mapping. |
| **SpacePAC** (Bioconductor) | Spatial protein amino-acid clustering combining COSMIC + PDB. | COSMIC somatic counts | Same somatic-recurrence paradigm. |
| **e-Driver3D / e-Driver** | Detects clustering of mutations in protein regions/3D to find drivers. | Somatic mutations | Driver detection, recurrence-based. |
| **3DHotspots.org** (Gao et al.) | Web resource of pre-computed 3D hotspot mutations across cancers. | Somatic mutations (pan-cancer) | A *database of results*, not a tool you run on your own screen. |
| **Oncodrive3D** (Bertran et al., bioRxiv 2025) | Fast, accurate detection of 3D structural mutation clusters under positive selection (AlphaFold-based). | Somatic mutations | Newest-generation somatic 3D-cluster caller; again recurrence/selection input, not screen LFC. |
| **POSTAR / PhosphoSitePlus-3D (spatial PTM–mutation)** | Spatial co-localization of PTM sites and mutations. | PTM + mutation annotations | Annotation co-localization, not screen-signal aggregation. |

**The defining difference (whole category):** Somatic-mutation 3D clustering asks *"are mutations positionally non-random across a patient cohort?"* — the signal is **discrete recurrence**. BE3D asks *"where does a continuous, designed functional readout (dropout LFC) concentrate in 3D?"* — the signal is **a per-residue effect size from a controlled perturbation**, with **known negative/neutral controls** (silent / no-mutation guides) enabling a *screen-specific randomization null* (BE3D's `randomize_data`/`randomize_sequence`) rather than a generic mutation-density null. That control structure, plus the base-editing-specific parsing (bystander edit types aggregated per residue via mean/min/max/sum) and the **QA gate** (`hypothesis_test`: MWU + KS on knockout vs neutral guide LFC), is BE3D's methodological distinction from the entire cancer-hotspot lineage.

### 1d. Structure-function mapping / variant-to-structure portals & effect predictors

These map *variants or scores* onto structure, or predict variant effect. They are complements (BE3D's outputs are *consumed by* the portal in this list) or orthogonal predictors.

| Tool | What it does | Relationship to BE3D |
|---|---|---|
| **Genomics 2 Proteins Portal (G2P)** (Kwon et al., Nat Methods 2024; g2p.broadinstitute.org) | Human-proteome resource mapping ~20M variants onto ~42k sequences / ~78k structures; **interactive upload** of custom per-residue annotations/scores onto a structure. | **Downstream partner / output target** of BE3D — BE3D emits a G2P-formatted TSV (`g2p_formatted_hit_cluster`). G2P provides the interactive 3D viz BE3D itself does not build. Same Broad ecosystem. |
| **MuPIT / MutPanning / CRAVAT** (KarchinLab) | Interactive mapping of variants onto 3D structures; MutPanning finds positively-selected genes. | Viz/driver-detection; not screen-LFC spatial aggregation. |
| **ProtVar** (EMBL-EBI) | Maps & interprets human missense variation on structure with functional context. | Variant annotation/interpretation portal; complement. |
| **AlphaMissense** (Cheng et al., Science 2023) | Proteome-wide missense pathogenicity from AlphaFold-based model (auROC ~0.94 on ClinVar). | Orthogonal *predictor*; can be overlaid alongside BE3D hotspots but predicts pathogenicity from sequence/structure, not from a screen. |
| **EVE** (Frazer et al., Nature 2021) | Unsupervised evolutionary variant effect from deep generative model on MSAs. | Orthogonal predictor (evolutionary). |
| **ESM1b / ESM-based VEP** | Protein-LM zero-shot variant effect. | Orthogonal predictor. |
| **Missense3D / Missense3D-TM** | Structural damage prediction for a missense change (steric/burial/etc.). | Per-variant structural consequence, not spatial aggregation of screen signal. |
| **FoldX / ddG tools (Rosetta, ThermoMPNN)** | Predict ΔΔG stability change of a mutation on a structure. | Biophysical per-mutation prediction; orthogonal, could annotate BE3D hits. |
| **PhosphoSitePlus (mapping)** | PTM site annotations. | Annotation source for BE3D's enrichment/characterization step. |

**Gap:** These predict or display *per-variant* properties; none aggregate a *screen's* continuous signal over 3D neighborhoods and call significant spatial clusters. BE3D's `enrichment_test` (Fisher's exact for structural features among hits) and `pLDDT_RSA_scatter` actually *consume* the kind of structural annotations these tools/portals provide.

### 1e. Structure SOURCES (dependencies, not competitors)

| Source | Role for BE3D |
|---|---|
| **AlphaFold DB** | Default structure fetched by UniProt ID (isoform-aware); primary dependency. |
| **PDB** | Experimental structures for custom `user_pdb`; PPI chains for interaction-aware LFC3D. |
| **ESMFold** | Alternative single-sequence predicted structure (usable as custom PDB). |
| **Boltz / AlphaFold3-class** | Newer co-folding/complex predictors — could supply complex structures for PPI-aware LFC3D. |
| **DSSP** | Secondary structure / RSA annotation (used in characterization; optional). |
| **UniProt** | Canonical sequence + isoform resolution. |
| **MUSCLE / Clustal Omega** | Pairwise alignment for cross-species / cross-isoform conservation & meta-aggregation. |

These are inputs BE3D orchestrates; it is a *consumer*, so its quality is bounded by AlphaFold/PDB coverage and pLDDT for the target.

---

## 2. Feature-by-feature comparison matrix

Columns: **Input** (what signal it eats) · **3D req** (needs a 3D structure) · **3D-clust** (does spatial 3D clustering/hotspots) · **Meta** (multi-screen/cohort meta-aggregation) · **Consv** (cross-species/conservation) · **Null** (statistical significance via null/permutation model) · **Viz** (interactive visualization) · **BE-specific** (base-editing-aware) · **Interface** (GUI vs code) · **Open** (free/open source) · **Maint** (actively maintained, as of 2026) · **Curve** (learning curve).

| Tool | Input | 3D req | 3D-clust | Meta | Consv | Null | Viz | BE-specific | Interface | Open | Maint | Curve |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| **BE3D / beclust3d** | BE tiling screen LFC (TSV) | Yes | **Yes** (agglom.) | **Yes** | **Yes** | Yes (randomization) | Via G2P export | **Yes** | Code (+Colab) | Yes (MIT) | Yes (v1.0, 2026) | Medium |
| **ProTiler-Mut** | CRISPR tiling mutagenesis screens (multi-condition) | Yes | **Yes** (3D-RRA substructures) | Multi-condition | Partial | Yes (RRA-based) | Plots | Cas9-KO oriented | Code | Yes | Yes (2026) | Medium-high |
| **CRISPRO** | Guide functional scores | Renders on struct. | No (1D LOESS + render) | No | Yes | Limited | Static/PyMOL out | No (Cas9) | Code | Yes | Low (2018-era) | Medium |
| **ProTiler** | Cas9 tiling KO dropout | No | No (1D regions) | No | Uses domains/PTM | Yes (region calling) | Plots | No | Code | Yes | Low | Medium |
| **CRISPR-SURF** | Tiling screen signal | No | No (1D deconv) | No | No | Yes (empirical) | Web + code | No | Web + code | Yes | Low | Low-medium |
| **BEAN (crispr-bean)** | BE screen counts + reporter | No | No | Shares across guides | No | Yes (Bayesian) | Plots | **Yes** | Code | Yes | Yes | Medium-high |
| **MAGeCK / VISPR** | sgRNA read counts | No | No | Multi-condition (MLE) | No | Yes (RRA/MLE) | VISPR web | No | Code + web | Yes | Yes | Medium |
| **BAGEL2** | sgRNA fold-changes | No | No | No | No | Yes (Bayes factor) | Plots | No | Code | Yes | Yes | Low-medium |
| **drugZ** | sgRNA read counts | No | No | Paired conditions | No | Yes (normal Z) | Plots | No | Code | Yes | Yes | Low |
| **HotMAPS** | Somatic mutation counts | Yes | **Yes** (density) | Pan-cohort | No | Yes (permutation) | Web resource | No | Code + web | Yes | Moderate | Medium-high |
| **HotSpot3D** | Somatic mutation counts | Yes | **Yes** (graph cluster) | Pan-cohort | No | Yes | Web server | No | Code + web | Yes | Moderate | Medium-high |
| **CLUMPS** | Somatic mutation freq | Yes | **Yes** (WAP stat) | Pan-cancer | No | Yes (permutation) | No | No | Code | Yes | Low | High |
| **Mutation3D** | Coding variants | Yes | **Yes** (atomic cluster) | Cohort | No | Yes | Web server | No | Web + code | Yes | Low | Low (web) |
| **Oncodrive3D** | Somatic mutations | Yes (AlphaFold) | **Yes** (3D clusters) | Cohort | No | Yes (selection) | Plots | No | Code | Yes | Yes (2025) | Medium |
| **G2P Portal** | Any per-residue scores | Yes | No (displays) | No (viewer) | Shows features | No | **Yes (interactive)** | No | Web GUI | Yes (portal) | Yes (2024) | Low |
| **AlphaMissense** | Sequence (variant) | Model-internal | No | N/A | Model-internal | Calibrated | Via portals | No | Precomputed + code | Yes (data/code) | Yes | Low (lookup) |
| **Missense3D** | Single variant + struct. | Yes | No | No | No | Rule-based | Web | No | Web | Yes | Moderate | Low |
| **BE-Hive** | Target sequence | No | No | N/A | No | ML-calibrated | Web | **Yes (predict)** | Web + code | Yes | Low | Low |

```json
[
  {"tool":"BE3D / beclust3d","category":"BE-screen 3D hotspot interpretation","input":"BE tiling screen LFC (TSV)","needs_3d_structure":true,"spatial_3d_clustering":true,"multi_screen_meta":true,"cross_species_conservation":true,"statistical_null_model":"randomization/permutation","interactive_viz":"via G2P export","base_editing_specific":true,"interface":"code + Colab","free_open":true,"actively_maintained":true,"learning_curve":"medium","url":"https://github.com/broadinstitute/BE3D"},
  {"tool":"ProTiler-Mut","category":"tiling-screen 3D hotspot interpretation (DIRECT competitor)","input":"CRISPR tiling mutagenesis screens (multi-condition)","needs_3d_structure":true,"spatial_3d_clustering":true,"multi_screen_meta":"multi-condition","cross_species_conservation":"partial","statistical_null_model":"3D-RRA (robust rank aggregation)","interactive_viz":"plots","base_editing_specific":false,"interface":"code","free_open":true,"actively_maintained":true,"learning_curve":"medium-high","url":"https://www.cell.com/cell-systems/abstract/S2405-4712(26)00133-X"},
  {"tool":"CRISPRO","category":"guide-score to protein/structure mapping","input":"guide functional scores","needs_3d_structure":"renders on structure","spatial_3d_clustering":false,"multi_screen_meta":false,"cross_species_conservation":true,"statistical_null_model":"limited","interactive_viz":"static/PyMOL output","base_editing_specific":false,"interface":"code","free_open":true,"actively_maintained":false,"learning_curve":"medium","url":"https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1563-5"},
  {"tool":"ProTiler","category":"1D tiling functional-region calling","input":"Cas9 tiling KO dropout","needs_3d_structure":false,"spatial_3d_clustering":false,"multi_screen_meta":false,"cross_species_conservation":"domains/PTM features","statistical_null_model":"region-calling","interactive_viz":"plots","base_editing_specific":false,"interface":"code","free_open":true,"actively_maintained":false,"learning_curve":"medium","url":"https://github.com/MDhewei/protiler"},
  {"tool":"CRISPR-SURF","category":"1D tiling signal deconvolution","input":"tiling screen signal","needs_3d_structure":false,"spatial_3d_clustering":false,"multi_screen_meta":false,"cross_species_conservation":false,"statistical_null_model":"empirical/permutation","interactive_viz":"web + code","base_editing_specific":false,"interface":"web + code","free_open":true,"actively_maintained":false,"learning_curve":"low-medium","url":"https://www.nature.com/articles/s41592-018-0225-6"},
  {"tool":"BEAN (crispr-bean)","category":"upstream BE variant-effect quantification","input":"BE screen counts + reporter editing","needs_3d_structure":false,"spatial_3d_clustering":false,"multi_screen_meta":"info-sharing across guides","cross_species_conservation":false,"statistical_null_model":"Bayesian network","interactive_viz":"plots","base_editing_specific":true,"interface":"code","free_open":true,"actively_maintained":true,"learning_curve":"medium-high","url":"https://github.com/pinellolab/crispr-bean"},
  {"tool":"MAGeCK / MAGeCK-VISPR","category":"upstream screen caller","input":"sgRNA read counts","needs_3d_structure":false,"spatial_3d_clustering":false,"multi_screen_meta":"multi-condition MLE","cross_species_conservation":false,"statistical_null_model":"RRA / MLE","interactive_viz":"VISPR web","base_editing_specific":false,"interface":"code + web","free_open":true,"actively_maintained":true,"learning_curve":"medium","url":"https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0843-6"},
  {"tool":"BAGEL2","category":"upstream essentiality caller","input":"sgRNA fold-changes","needs_3d_structure":false,"spatial_3d_clustering":false,"multi_screen_meta":false,"cross_species_conservation":false,"statistical_null_model":"Bayes factor","interactive_viz":"plots","base_editing_specific":false,"interface":"code","free_open":true,"actively_maintained":true,"learning_curve":"low-medium","url":"https://github.com/hart-lab/bagel"},
  {"tool":"drugZ","category":"upstream drug-gene interaction caller","input":"sgRNA read counts","needs_3d_structure":false,"spatial_3d_clustering":false,"multi_screen_meta":"paired conditions","cross_species_conservation":false,"statistical_null_model":"normalized Z","interactive_viz":"plots","base_editing_specific":false,"interface":"code","free_open":true,"actively_maintained":true,"learning_curve":"low","url":"https://github.com/hart-lab/drugz"},
  {"tool":"HotMAPS","category":"somatic 3D mutation hotspot","input":"somatic mutation counts","needs_3d_structure":true,"spatial_3d_clustering":true,"multi_screen_meta":"pan-cohort","cross_species_conservation":false,"statistical_null_model":"permutation density","interactive_viz":"web resource","base_editing_specific":false,"interface":"code + web","free_open":true,"actively_maintained":"moderate","learning_curve":"medium-high","url":"https://github.com/KarchinLab/HotMAPS"},
  {"tool":"HotSpot3D","category":"somatic 3D mutation hotspot","input":"somatic mutation counts","needs_3d_structure":true,"spatial_3d_clustering":true,"multi_screen_meta":"pan-cohort","cross_species_conservation":false,"statistical_null_model":"clustering + permutation","interactive_viz":"web server","base_editing_specific":false,"interface":"code + web","free_open":true,"actively_maintained":"moderate","learning_curve":"medium-high","url":"https://academic.oup.com/bioinformatics/article/36/12/3944/5823299"},
  {"tool":"CLUMPS","category":"somatic 3D mutation clustering","input":"somatic mutation frequency","needs_3d_structure":true,"spatial_3d_clustering":true,"multi_screen_meta":"pan-cancer","cross_species_conservation":false,"statistical_null_model":"WAP permutation","interactive_viz":false,"base_editing_specific":false,"interface":"code","free_open":true,"actively_maintained":false,"learning_curve":"high","url":"https://www.pnas.org/doi/10.1073/pnas.1516373112"},
  {"tool":"Mutation3D","category":"somatic 3D atomic clustering / driver","input":"coding variants","needs_3d_structure":true,"spatial_3d_clustering":true,"multi_screen_meta":"cohort","cross_species_conservation":false,"statistical_null_model":"clustering significance","interactive_viz":"web server","base_editing_specific":false,"interface":"web + code","free_open":true,"actively_maintained":false,"learning_curve":"low","url":"https://pubmed.ncbi.nlm.nih.gov/26841357/"},
  {"tool":"Oncodrive3D","category":"somatic 3D structural cluster / selection","input":"somatic mutations","needs_3d_structure":true,"spatial_3d_clustering":true,"multi_screen_meta":"cohort","cross_species_conservation":false,"statistical_null_model":"positive-selection null","interactive_viz":"plots","base_editing_specific":false,"interface":"code","free_open":true,"actively_maintained":true,"learning_curve":"medium","url":"https://www.biorxiv.org/content/10.1101/2025.01.11.632354"},
  {"tool":"Genomics 2 Proteins Portal (G2P)","category":"variant-to-structure interactive portal (BE3D OUTPUT target)","input":"any per-residue scores","needs_3d_structure":true,"spatial_3d_clustering":false,"multi_screen_meta":false,"cross_species_conservation":"shows features","statistical_null_model":false,"interactive_viz":true,"base_editing_specific":false,"interface":"web GUI","free_open":true,"actively_maintained":true,"learning_curve":"low","url":"https://g2p.broadinstitute.org/"},
  {"tool":"AlphaMissense","category":"missense pathogenicity predictor (orthogonal)","input":"sequence variant","needs_3d_structure":"model-internal","spatial_3d_clustering":false,"multi_screen_meta":false,"cross_species_conservation":"model-internal","statistical_null_model":"calibrated","interactive_viz":"via portals","base_editing_specific":false,"interface":"precomputed + code","free_open":true,"actively_maintained":true,"learning_curve":"low","url":"https://www.science.org/doi/10.1126/science.adg7492"},
  {"tool":"Missense3D","category":"structural consequence of a missense variant","input":"single variant + structure","needs_3d_structure":true,"spatial_3d_clustering":false,"multi_screen_meta":false,"cross_species_conservation":false,"statistical_null_model":"rule-based","interactive_viz":"web","base_editing_specific":false,"interface":"web","free_open":true,"actively_maintained":"moderate","learning_curve":"low","url":"http://missense3d.bc.ic.ac.uk/"},
  {"tool":"BE-Hive","category":"BE outcome/efficiency predictor (upstream, orthogonal)","input":"target sequence","needs_3d_structure":false,"spatial_3d_clustering":false,"multi_screen_meta":false,"cross_species_conservation":false,"statistical_null_model":"ML-calibrated","interactive_viz":"web","base_editing_specific":true,"interface":"web + code","free_open":true,"actively_maintained":false,"learning_curve":"low","url":"https://www.crisprbehive.design/"}
]
```

---

## 3. Where BE3D WINS, LOSES, and what is genuinely UNIQUE

### Genuinely UNIQUE (defensible niche)
1. **Only openly-packaged pipeline purpose-built to turn a *base-editing* tiling screen into 3D structural hotspots end-to-end.** The combination of (a) BE-specific parsing of bystander/multi-edit guides aggregated per residue, (b) a **spatial-neighborhood LFC3D statistic**, (c) a **screen-aware randomization null** exploiting the screen's own silent/no-mutation controls, (d) **agglomerative 3D clustering into hotspots**, (e) **cross-screen and cross-species/isoform meta-aggregation** (BE-MetaClust3D), and (f) a **built-in QA gate** (`hypothesis_test`) plus (g) **direct G2P interactive-portal export** — that *whole stack* exists in no single other tool.
2. **Signal type is the deep differentiator.** Every 3D-clustering analogue (HotMAPS, HotSpot3D, CLUMPS, Mutation3D, Oncodrive3D) is built for *somatic mutation recurrence* under positive selection. BE3D is built for a *designed, controlled, continuous functional readout* with real negative controls — a fundamentally cleaner statistical setting the cancer tools cannot exploit.
3. **Meta-aggregation to rescue weak signal across noisy BE screens** (SUM across screens, with alignment for cross-species) is a distinctive contribution not present in the somatic-hotspot lineage or in CRISPRO/ProTiler/SURF.
4. **PPI-chain-aware LFC3D** (aggregation across interacting protein chains) extends hotspot detection to interfaces — beyond single-chain somatic tools.

### Where BE3D WINS
- **vs upstream callers (MAGeCK/BAGEL/drugZ/BEAN):** not competing — BE3D adds a 3D layer none of them have; it *complements* them (ideal: BEAN → BE3D).
- **vs 1D interpreters (ProTiler, CRISPR-SURF):** BE3D captures spatially-proximal-but-sequence-distant hotspots (allosteric sites, interfaces, folding cores) that 1D methods structurally cannot see.
- **vs CRISPRO:** CRISPRO renders on structure but has no spatial-neighborhood statistic, no significance-tested 3D clusters, no meta-aggregation. BE3D is the statistical/clustering successor.
- **vs somatic 3D tools:** BE3D uses the screen's controls for a tailored null and handles BE edit semantics; the cancer tools would need re-engineering to accept LFC.
- **Openness + accessibility:** MIT license, pip-installable, zero-install Colab path, YAML-driven local runs, and a clean handoff to a maintained interactive portal (G2P). Low barrier relative to CLUMPS/HotMAPS setup.

### Where BE3D LOSES / risks
- **Direct competition from ProTiler-Mut** (He et al., Cell Systems 2026): the *only* other tool doing tiling-screen→3D-hotspot. ProTiler-Mut uses a 3D-RRA (robust rank aggregation) statistic and emphasizes multi-condition phenotypes and mutation-disrupted PPIs. BE3D must be positioned against it explicitly (BE3D's edges: base-editing-native parsing, randomization null, cross-species/isoform meta, G2P integration; ProTiler-Mut's edges: rank-based robustness, multi-condition modeling, lineage from the established ProTiler brand). This head-to-head is the key vetting question.
- **No native interactive visualization** — depends on exporting to G2P; competitors like HotSpot3D/Mutation3D ship their own web viewers, and G2P dependency ties BE3D to one ecosystem.
- **Structure-quality ceiling:** results are bounded by AlphaFold/PDB availability and pLDDT; low-confidence or disordered regions and multi-domain/complex assemblies weaken spatial neighborhoods. Somatic tools face this too, but BE3D's per-screen scope makes a single bad structure decisive.
- **Statistical maturity/benchmarking:** the somatic-hotspot field has a decade of benchmarking, FDR conventions, and pan-cohort validation; BE3D's randomization null and single-linkage/euclidean clustering choices are reasonable but comparatively unvalidated across many genes, and single-linkage clustering can chain. Reviewers will want calibration/benchmark evidence.
- **Not a predictor:** produces hypotheses from an *existing* screen; offers nothing pre-experiment (unlike BE-Hive/AlphaMissense), so its audience is narrower (labs that already ran a BE tiling screen).
- **Maturity/adoption:** v1.0.0, new (2026); ecosystem tools (MAGeCK, BAGEL, AlphaMissense) have large user bases and citations. Adoption risk.

### Bottom line
BE3D occupies a **real, narrow, defensible niche**: the base-editing-native, control-aware, meta-aggregating 3D-hotspot interpreter for tiling screens, wired into an interactive portal. Its novelty is genuine relative to the cancer-hotspot lineage and to 1D screen tools, and it is complementary (not competitive) with upstream callers and with predictors. The **one serious novelty threat is ProTiler-Mut**, the concurrent tiling-screen→3D-hotspot method; the vetting should center on a direct methodological and empirical comparison there.

---

## Sources
- CRISPRO / ProTiler / CRISPR-SURF: https://www.nature.com/articles/s41467-019-12489-8 ; https://github.com/MDhewei/protiler ; https://www.researchgate.net/publication/329324948
- ProTiler-Mut / BE3D-adjacent Cell Systems paper: https://www.cell.com/cell-systems/abstract/S2405-4712(26)00133-X ; bioRxiv https://www.biorxiv.org/content/10.1101/2025.04.17.649336
- BEAN: https://github.com/pinellolab/crispr-bean ; https://www.nature.com/articles/s41588-024-01726-6
- HotMAPS: https://github.com/KarchinLab/HotMAPS ; https://aacrjournals.org/cancerres/article/76/13/3719/607886
- HotSpot3D: https://academic.oup.com/bioinformatics/article/36/12/3944/5823299
- CLUMPS / Mutation3D / SpacePAC / 3Dhotspots: https://pubmed.ncbi.nlm.nih.gov/26841357/ ; https://link.springer.com/article/10.1186/s13073-016-0393-x ; https://link.springer.com/article/10.1186/1471-2105-15-231
- Oncodrive3D: https://www.biorxiv.org/content/10.1101/2025.01.11.632354
- Multiscale functional map (structure + network): https://www.nature.com/articles/s41467-024-54176-3
- BE outcome predictors (BE-Hive/BEdeepon/DeepBE/BE-FF): https://www.crisprbehive.design/ ; https://www.nature.com/articles/s41467-021-25375-z ; https://www.biorxiv.org/content/10.1101/2021.03.14.435303v1.full ; https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7319459/
- MAGeCK-VISPR / BAGEL2 / drugZ / review: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0843-6 ; https://portlandpress.com/emergtoplifesci/article/5/6/779/230423 ; https://nf-co.re/crisprseq
- Genomics 2 Proteins Portal: https://www.nature.com/articles/s41592-024-02409-0 ; https://g2p.broadinstitute.org/ ; https://www.biorxiv.org/content/10.1101/2024.01.02.573913v1
- AlphaMissense / EVE / Missense3D: https://www.science.org/doi/10.1126/science.adg7492 ; https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7617522/
