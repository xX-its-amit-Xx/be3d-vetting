# BE3D — Personas, Use Cases & Market Impact (Queen 2)

**Tool:** BE3D (`broadinstitute/BE3D`, pip `beclust3d` v1.0.0, MIT). Python package that turns a *completed* base-editor (BE) tiling-mutagenesis functional screen into **3D structural hotspots** of functionally important residues via LFC3D spatial scoring + agglomerative clustering + multi-screen meta-aggregation, with output to the Genomics 2 Proteins (G2P) Portal.

**Key positioning:** downstream *interpretation/analysis* layer. NOT a BE outcome/efficiency predictor (BE-Hive/BE-FF), NOT a guide designer, NOT a variant pathogenicity predictor (AlphaMissense/ESM). It sits between functional-genomics screens and structural biology.

---

## 1. PERSONA CARDS

Each card: role, goals, input data, why BE3D helps, alternatives, decision criteria, sophistication, adoption likelihood.

### ACADEMIC

---
**A1 — Functional-genomics / CRISPR-screen lab PI (mammalian cell)**
- **Role:** PI running pooled BE tiling screens (proliferation/dropout, drug-resistance) in cancer cell lines.
- **Goals:** Convert noisy per-guide LFC into interpretable biology; publish structure-function stories; find the "why" behind hit residues.
- **Input:** Their own screen scores TSV (Mutation_type, Mutation_list, Gene, sgRNA_score) — exactly BE3D's native input. UniProt ID.
- **Why BE3D:** Purpose-built for their exact output; LFC3D rescues weak single-residue signal; QA module gates screen quality; publication-ready plots + G2P figures.
- **Alternatives:** Custom in-house scripts (PyMOL/ChimeraX manual mapping), MAGeCK for hit-calling only (no 3D), collaborate ad hoc with a structural biologist.
- **Decision criteria:** Fits their TSV without reformatting; reproducibility; figure quality; free/open.
- **Sophistication:** High (comp-savvy lab or has a bioinformatician).
- **Adoption: HIGH** — this is the bullseye persona; input is literally their screen output.

---
**A2 — Base-editing method developer (screen technology labs)**
- **Role:** Develops new BE screening methods / libraries (e.g. groups extending Hanna, Cuella-Martin, Sangree, Doench-style work).
- **Goals:** Demonstrate that a new screen design recovers known functional residues; benchmark signal; validate against structure.
- **Input:** Screens from many genes; often multiple editors (CBE/ABE) → meta-aggregation.
- **Why BE3D:** BE-QA quantifies screen quality; MetaClust3D aggregates across editors; provides a standard analysis to show method sensitivity.
- **Alternatives:** Bespoke benchmarking code; ROC vs known essential residues.
- **Decision criteria:** Statistical rigor (null/randomization, p-values), extensibility.
- **Sophistication:** Very high.
- **Adoption: HIGH** — natural early adopters; likely to cite/extend it. Broad GPP-adjacent.

---
**A3 — Structural biologist (interpreting functional data on structures)**
- **Role:** Crystallographer / cryo-EM / AlphaFold modeler wanting functional annotation on structures.
- **Goals:** Overlay functional screen data onto their structure; identify functional surfaces, interfaces, pockets.
- **Input:** Custom PDB + collaborator's screen data.
- **Why BE3D:** Accepts custom FASTA+PDB; DSSP secondary-structure + RSA enrichment; PPI-chain aware neighborhoods; maps function to 3D.
- **Alternatives:** Manual coloring in PyMOL/ChimeraX by B-factor; Consurf for conservation.
- **Decision criteria:** Custom structure support; control over radius; interpretability.
- **Sophistication:** High structurally, medium computationally.
- **Adoption: MED-HIGH** — strong fit when collaborating with a screen lab; less standalone.

---
**A4 — Cancer genomics / driver-gene lab**
- **Role:** Studies tumor suppressors/oncogenes; separates LOF vs GOF, allosteric vs active-site.
- **Goals:** Interpret saturation/tiling BE screens of a driver gene; find functional hotspots (e.g. degron, dimerization, catalytic sites); prioritize residues for mechanism.
- **Input:** BE screen of a driver (the KBTBD4/Yeo et al. example is exactly this).
- **Why BE3D:** Clusters significant residues into spatial hotspots that map to functional domains; split +/- scores distinguish activating vs inactivating.
- **Alternatives:** Hotspot mutation tools (cBioPortal 3D hotspots, HotMAPS) — but those use *patient mutation frequency*, not functional screens.
- **Decision criteria:** Handles saturation data; distinguishes effect direction; links to drug pockets.
- **Sophistication:** High.
- **Adoption: HIGH** — driver-gene structure-function is a major publication driver.

---
**A5 — Variant-effect / clinical-genetics researcher (MAVE/VUS)**
- **Role:** Builds functional evidence to reclassify variants of uncertain significance (VUS).
- **Goals:** Aggregate residue-level functional signal spatially; provide ACMG PS3 functional evidence; spatial context for VUS.
- **Input:** BE/SGE tiling screen of a clinically actionable gene (e.g. BRCA1, BARD1).
- **Why BE3D:** Adds 3D spatial corroboration to per-variant scores; hotspot membership as supporting evidence; conservation module.
- **Alternatives:** SGE per-variant function maps (Findlay), AlphaMissense/ESM1b, direct MAVE scores. BE3D complements rather than replaces these.
- **Decision criteria:** Interpretability for clinical reporting; statistical defensibility.
- **Sophistication:** Medium-high.
- **Adoption: MED** — BE screens are coarser than SGE for clinical VUS, but growing; complementary.

---
**A6 — Protein engineer / directed-evolution academic**
- **Role:** Engineers enzymes/proteins; defines functional/tolerant sites.
- **Goals:** Map which residues are functionally critical vs tolerant; define domains and hotspots for library design.
- **Input:** BE tiling (or DMS) functional screen of the protein.
- **Why BE3D:** Identifies spatial clusters of critical residues; RSA/secondary-structure enrichment guides where to mutate.
- **Alternatives:** DMS heatmaps, Rosetta, in-house tolerance maps.
- **Decision criteria:** Works on non-disease proteins; custom PDB.
- **Sophistication:** High.
- **Adoption: MED** — DMS (not BE) is their usual data; fit depends on data type.

---
**A7 — Bioinformatics core facility / staff scientist**
- **Role:** Services many PIs; runs standardized pipelines.
- **Goals:** Offer a repeatable "screen → 3D hotspots" service; reduce bespoke scripting.
- **Input:** Whatever screens clients bring.
- **Why BE3D:** Pip-installable, YAML-configured, Colab option, standard I/O — easy to templatize.
- **Alternatives:** Build/maintain their own; nf-core-style pipelines.
- **Decision criteria:** Install/maintenance burden, docs, deterministic outputs, support.
- **Sophistication:** Very high.
- **Adoption: MED-HIGH** — cores love turnkey packages if robust; they become force-multipliers.

---
**A8 — Comp-bio / method-development lab (build-on)**
- **Role:** Develops variant-effect/structure ML; wants features or benchmarks.
- **Goals:** Use LFC3D scores as features/labels; benchmark structure-aware aggregation.
- **Input:** Public BE screens + AlphaFold.
- **Why BE3D:** Ready spatial-aggregation implementation with a null model; MIT license.
- **Alternatives:** Reimplement graph/neighborhood aggregation themselves.
- **Decision criteria:** Code quality, extensibility, license.
- **Sophistication:** Very high.
- **Adoption: MED** — will fork/extend rather than use as-is.

---
**A9 — Grad student / postdoc (hands-on operator across A1–A6)**
- **Role:** The person who actually runs the pipeline.
- **Goals:** Get figures for thesis/paper fast; low setup friction.
- **Why BE3D:** Colab zero-install; example (KBTBD4); good plots.
- **Alternatives:** Ask a bioinformatician; hack scripts.
- **Decision criteria:** Tutorials, Colab, worked examples, error tolerance.
- **Sophistication:** Medium (variable).
- **Adoption: HIGH** — Colab + example lowers the barrier dramatically; students drive uptake.

---
**A10 — Cross-species / evolutionary functional-genomics lab**
- **Role:** Compares functional screens across species/orthologs/isoforms.
- **Goals:** Test conservation of functional hotspots across species.
- **Input:** Screens of orthologs; two sequences to align.
- **Why BE3D:** Built-in MUSCLE/CLUSTAL conservation + cross-species meta-aggregation.
- **Alternatives:** Manual alignment + custom mapping.
- **Decision criteria:** Alignment handling, robustness of API fallback.
- **Sophistication:** High.
- **Adoption: MED** — niche but distinctive; few tools do this natively.

---
**A11 — Drug-target-discovery academic (chemical biology / target validation)**
- **Role:** Validates targets, finds druggable/allosteric sites.
- **Goals:** Map screen hits onto known drug pockets; find allosteric functional sites.
- **Input:** BE screen of a candidate target + liganded PDB.
- **Why BE3D:** Hotspots on a custom liganded structure highlight functionally essential pocket-adjacent residues.
- **Alternatives:** Fpocket/SiteMap + manual overlay.
- **Decision criteria:** Custom PDB with ligand; interpretability.
- **Sophistication:** High.
- **Adoption: MED-HIGH** — strong story for target validation papers.

---
### INDUSTRY

---
**I1 — Biotech/pharma target-discovery / functional-genomics group**
- **Role:** Industrial screening teams (large pharma FG groups, DepMap-style).
- **Goals:** Prioritize residues/targets for programs; mechanistic support for target dossiers; find allosteric handles.
- **Input:** Proprietary BE tiling screens at scale.
- **Why BE3D:** Standardized 3D interpretation; MetaClust3D consolidates many internal screens; de-risks target selection.
- **Alternatives:** Internal comp platforms; CRO analysis; bespoke.
- **Decision criteria:** IP/data security (must run locally — BE3D supports local), validation, integration, support/SLA (concern for a small academic package).
- **Sophistication:** Very high.
- **Adoption: MED** — value is high, but enterprises favor supported/internal tooling; may adopt the *method*/reimplement. Local mode helps.

---
**I2 — Base-editing therapeutics company (Beam, Verve, Prime, Intellia, etc.)**
- **Role:** Companies whose core tech IS base editing.
- **Goals:** Understand target-protein structure-function; select edit sites; interpret internal BE tiling of a target; safety/mechanism.
- **Input:** In-house BE screens; targets like PCSK9, HBG promoter, SERPINA1.
- **Why BE3D:** Native to their assay type; from the Broad (credibility/adjacency); LFC3D gives mechanistic maps of edit consequences.
- **Alternatives:** Internal platforms (well-resourced), academic collaborations.
- **Decision criteria:** Data security, regulatory-grade reproducibility, support.
- **Sophistication:** Very high.
- **Adoption: MED** — likely to evaluate/borrow method; may not adopt an academic package wholesale for GxP work. Strong for exploratory research.

---
**I3 — Contract Research Organization (CRO) offering screening services**
- **Role:** Runs BE/CRISPR screens for clients; sells analysis as add-on.
- **Goals:** Differentiate deliverables with 3D hotspot reports.
- **Input:** Client screens.
- **Why BE3D:** Turnkey, standardized, presentable outputs; adds a premium analysis tier cheaply.
- **Alternatives:** Build own; deliver raw hit lists only.
- **Decision criteria:** Robustness at volume, white-labelable outputs, licensing (MIT is friendly).
- **Sophistication:** High.
- **Adoption: MED-HIGH** — cheap differentiation on top of existing services.

---
**I4 — Ag-bio / crop & livestock editing company**
- **Role:** Base editing in plants/animals for trait improvement.
- **Goals:** Identify functionally critical residues to edit or avoid; cross-species hotspot transfer.
- **Input:** Trait/functional screens; often non-model organisms (custom FASTA+PDB, AlphaFold).
- **Why BE3D:** Custom structure + cross-species alignment; AlphaFold coverage of non-model proteins.
- **Alternatives:** Bespoke; limited off-the-shelf options.
- **Decision criteria:** Non-human protein support, structure availability.
- **Sophistication:** Medium-high.
- **Adoption: LOW-MED** — fewer such screens; but distinctive niche fit.

---
**I5 — Synthetic biology / enzyme engineering firm**
- **Role:** Engineers enzymes/biologics; industrial protein optimization.
- **Goals:** Map critical vs tolerant residues; guide library design.
- **Input:** Functional screens (often DMS more than BE).
- **Why BE3D:** Spatial critical-residue maps; enrichment by structural feature.
- **Alternatives:** In-house DMS pipelines, Rosetta, ML models.
- **Decision criteria:** Data-type fit (BE vs DMS), custom PDB.
- **Sophistication:** Very high.
- **Adoption: LOW-MED** — data often DMS not BE; method transferable but not native.

---
**I6 — Antibody / protein-engineering firm**
- **Role:** Optimizes antibodies/binders.
- **Goals:** Epitope/paratope functional mapping; developability hotspots.
- **Input:** Functional/affinity screens + structures.
- **Why BE3D:** PPI-chain-aware 3D aggregation for interface residues.
- **Alternatives:** Specialized antibody platforms.
- **Adoption: LOW** — data type and workflow diverge; niche.

---
**I7 — Clinical dx / variant-interpretation company**
- **Role:** Sells VUS classification / functional evidence.
- **Goals:** Incorporate functional + spatial evidence into variant reports.
- **Input:** Public + proprietary MAVE/BE screens.
- **Why BE3D:** Spatial corroboration layer for functional calls.
- **Alternatives:** AlphaMissense/ESM, direct MAVE scores, in-house pipelines.
- **Decision criteria:** Clinical validity, auditability.
- **Sophistication:** High.
- **Adoption: LOW-MED** — regulated setting; would likely reimplement/validate internally.

---
**I8 — AI-bio / techbio startup (features & benchmarks)**
- **Role:** Builds ML models for variant effect / protein function.
- **Goals:** LFC3D as training features/labels; benchmark data.
- **Why BE3D:** Open, structure-aware aggregation; complements AlphaFold pipelines.
- **Alternatives:** Reimplement; use raw screens.
- **Adoption: LOW-MED** — will cherry-pick the method/data.

---
**I9 — Tool / SaaS / portal vendor (integration partner)**
- **Role:** Structure-viewer / genomics-portal vendors (incl. G2P itself, Benchling-like platforms).
- **Goals:** Offer BE screen → 3D hotspot as a feature.
- **Why BE3D:** MIT-licensed engine; G2P-formatted output already exists.
- **Adoption: MED** — natural for G2P ecosystem; others may integrate.

---
## 2. PERSONA × USE-CASE MATRIX

Use cases:
- **U1** Prioritize residues for follow-up validation
- **U2** QA / gate screen quality before interpretation
- **U3** Interpret saturation/tiling BE screen of a tumor suppressor/driver
- **U4** Find allosteric / functional-surface sites
- **U5** Map screen hits to a drug/ligand pocket
- **U6** Define domains/critical residues for protein engineering
- **U7** Cross-species / isoform conservation of hotspots
- **U8** Meta-aggregate many screens → consensus hotspots
- **U9** Provide functional+spatial evidence for VUS classification
- **U10** Publication-ready 3D figures / G2P interactive sharing

| Persona | U1 | U2 | U3 | U4 | U5 | U6 | U7 | U8 | U9 | U10 |
|---|---|---|---|---|---|---|---|---|---|---|
| A1 FG/CRISPR-screen PI | H | H | H | M | M | L | L | M | L | H |
| A2 BE method developer | H | H | M | M | L | L | M | H | M | H |
| A3 Structural biologist | M | L | M | H | H | M | M | L | L | H |
| A4 Cancer/driver-gene lab | H | M | H | H | H | L | L | M | M | H |
| A5 Variant-effect/clinical researcher | M | M | H | M | L | L | M | M | H | M |
| A6 Protein engineer (academic) | M | L | L | M | M | H | M | L | L | M |
| A7 Bioinformatics core | H | H | M | M | M | M | M | H | M | H |
| A8 Comp-bio method lab | M | M | M | M | L | M | M | H | M | L |
| A9 Grad student/postdoc (operator) | H | H | H | M | M | M | M | M | M | H |
| A10 Cross-species FG lab | M | M | M | M | L | M | H | H | L | H |
| A11 Drug-target-discovery academic | H | M | H | H | H | M | L | M | L | H |
| I1 Pharma target discovery | H | H | H | H | H | M | L | H | M | M |
| I2 BE therapeutics co. | H | H | M | H | H | M | L | H | M | M |
| I3 CRO screening services | H | H | M | M | M | M | M | H | L | H |
| I4 Ag-bio editing | M | M | M | M | L | H | H | M | L | M |
| I5 SynBio/enzyme eng. | M | L | L | M | M | H | M | L | L | M |
| I6 Antibody engineering | L | L | L | M | M | M | L | L | L | M |
| I7 Clinical dx / variant interp. | M | M | H | L | L | L | M | M | H | M |
| I8 AI-bio startup | M | M | M | M | L | M | M | H | M | L |
| I9 Tool/SaaS/portal vendor | M | M | M | M | M | M | M | M | M | H |

(H/M/L = fit strength. Best-fit cells cluster on U1, U2, U3, U8, U10 for the core screen-lab personas.)

---
## 3. IMPACT / MARKET QUANTIFICATION

*Numbers are cited where available; clearly flagged as estimates otherwise.*

### 3a. Screening / functional-genomics market (BE3D's substrate)
- **CRISPR screening market:** ~**$1.52B (2025) → $1.79B (2026), $3.43B by 2030, ~17.7% CAGR.** [towardshealthcare; globenewswire]
- **Functional-genomics (CRISPR/Cas) segment:** **$691.6M (2024) → $1,618.7M (2030), 15.5% CAGR.** [Grand View Research]
- Growth drivers cited: precision drug discovery, gene-therapy R&D, bioinformatics-driven screening — i.e. exactly BE3D's demand pool.
- Broader CRISPR gene-editing market: **$4.76B (2025) → $18.89B (2035)**; another source: CRISPR market ~$3.27B in 2025. [towardshealthcare]

### 3b. Base editing as a fast-growing modality (the specific tailwind)
- Base-editor screens have seen "an explosion of studies" since the 2020 landmark pooled-screen papers (Hanna et al.; Cuella-Martin et al.) — the technique BE3D consumes is <6 years old and scaling. [Nature Genetics 2025; ScienceDirect Mol Cell 2023]
- Dedicated tooling now emerging (e.g. CRISPR-BEasy for BE tiling-library design, PMC 2025) — signals a maturing ecosystem BE3D can anchor the *analysis* end of.
- **Estimate (flagged):** Plausibly **hundreds of labs worldwide** run BE tiling/mutational-scanning screens today (low thousands run some form of pooled CRISPR screen). If ~200–500 labs run BE tiling screens and each does 1–5 screens/yr, that's a realistic **serviceable population of a few hundred groups / ~500–2,000 screens/yr** for BE3D. *(Order-of-magnitude estimate; no single authoritative count exists.)*

### 3c. Base-editing therapeutics industry (industry personas I2)
- **~332 gene-editing companies across 23 countries (Feb 2025); ~217 in the US (~two-thirds).** [SAMPS 2025] Base editing is a major sub-platform within this.
- **Beam Therapeutics:** launched a **~$500M** financing on positive BEAM-302 (AATD) Phase 1/2 in-vivo base-editing data; ex-vivo SCD program dosing (2024). [pharmaphorum; BioSpace]
- **Verve Therapeutics:** in-vivo base editors VERVE-101/-102 (cardiovascular, PCSK9/LDL) — VERVE-102 halved LDL at top dose in Ph1b; well tolerated. [BioSpace; CGTlive]
- **Prime Medicine:** prime-editing (adjacent), incl. an AATD program (2025) — arbitration dispute with Beam. [SEC filings; pharmaphorum]
- **Funding climate:** avg gene-editing financing round **~$68M in 2024**; only ~1% of biotech VC went to gene editing in 2024, but 2025 rebound (Tune $175M B, Light Horse $62M A). [Nature 2025]
- **Implication for BE3D:** a well-funded cluster of ~dozens of base-editing companies are the highest-value (if hardest-to-convert) industry targets; several are Broad/Beam-adjacent, aiding credibility.

### 3d. AlphaFold as an enabling tailwind (BE3D auto-fetches AlphaFold structures)
- **AlphaFold DB: >214M predicted structures** covering essentially all UniProt (2024). [PMC 10767828]
- **1.6M unique visitors from >190 countries; archive downloaded >23,000 times** (by Jan 2024); founding paper **>10,000 citations**. [PNAS 2024; arXiv 2403.02124]
- **Why it matters:** BE3D needs a 3D structure for every screened protein. AlphaFold's near-complete coverage means BE3D works out-of-the-box on almost any gene a screen targets — removing the historical bottleneck (no crystal structure). This directly expands BE3D's addressable target space to ~any protein-coding gene.

### 3e. Clinical-genetics / VUS tailwind (personas A5, I7)
- Saturation editing is proving decisive for VUS: **BARD1 SGE assessed 8,818 SNVs + 2,097 indels, discriminated pathogenic/benign at AUC 0.99, resolved 95.4% of VUS** (2025). [medRxiv 2025.11.03] BRCA1 SGE (Findlay) classified thousands of variants.
- BE3D adds the *spatial* interpretation layer atop such per-variant functional maps — a complementary, growing clinical use case.

### 3f. Addressable-population summary (estimates flagged)
- **Core academic bullseye (A1, A2, A4, A9):** a few hundred BE-screen / functional-genomics labs globally — **HIGH conversion**, near-zero switching cost (native input format, free, Colab).
- **Adjacent academic (A3, A5, A7, A10, A11):** another several hundred labs — **MED**, adoption when collaborating with a screen lab.
- **Industry (I1–I3, I9):** dozens–hundreds of pharma FG groups, ~30+ base-editing companies, CROs, portal vendors — **MED**, high value but favor internal/validated tools; likely to borrow the *method*.
- **Net:** primary near-term market is academic (open-source adoption + citations); industry is high-value but slower and partly method-transfer rather than package adoption. BE3D's growth is levered to two fast-rising curves: BE-screen adoption and AlphaFold ubiquity.

---
## 4. PERSONA JSON (for dashboard)

```json
[
  {"id":"A1","name":"Functional-genomics / CRISPR-screen lab PI","sector":"academic","role":"PI running pooled BE tiling screens","fit_score":5,"primary_use_case":"Interpret screen: per-residue LFC -> 3D hotspots + prioritize validation","would_use_instead":"Custom PyMOL/ChimeraX scripts, MAGeCK (hit-calling only)","adoption_likelihood":"High"},
  {"id":"A2","name":"Base-editing method developer","sector":"academic","role":"Screen-technology / library developer","fit_score":5,"primary_use_case":"Benchmark new screen designs; QA; meta-aggregate editors","would_use_instead":"Bespoke benchmarking code, ROC vs known residues","adoption_likelihood":"High"},
  {"id":"A3","name":"Structural biologist","sector":"academic","role":"Crystallography/cryo-EM/AlphaFold modeler","fit_score":4,"primary_use_case":"Overlay functional data on structure; find functional surfaces/interfaces","would_use_instead":"Manual B-factor coloring in PyMOL, ConSurf","adoption_likelihood":"Medium"},
  {"id":"A4","name":"Cancer genomics / driver-gene lab","sector":"academic","role":"Tumor-suppressor/oncogene mechanism lab","fit_score":5,"primary_use_case":"Interpret saturation BE screen of a driver; LOF/GOF hotspots","would_use_instead":"cBioPortal 3D hotspots, HotMAPS (mutation-frequency based)","adoption_likelihood":"High"},
  {"id":"A5","name":"Variant-effect / clinical-genetics researcher","sector":"clinical","role":"MAVE/VUS functional-evidence researcher","fit_score":3,"primary_use_case":"Spatial corroboration of functional scores for VUS (PS3)","would_use_instead":"SGE per-variant maps, AlphaMissense/ESM1b","adoption_likelihood":"Medium"},
  {"id":"A6","name":"Protein engineer / directed evolution (academic)","sector":"academic","role":"Enzyme/protein engineering PI","fit_score":3,"primary_use_case":"Map critical vs tolerant residues; define engineering hotspots","would_use_instead":"DMS heatmaps, Rosetta, in-house tolerance maps","adoption_likelihood":"Medium"},
  {"id":"A7","name":"Bioinformatics core facility / staff scientist","sector":"academic","role":"Core services bioinformatician","fit_score":4,"primary_use_case":"Standardized screen->3D-hotspot pipeline as a service","would_use_instead":"Build/maintain in-house pipeline","adoption_likelihood":"Medium"},
  {"id":"A8","name":"Computational biology / method-development lab","sector":"academic","role":"Variant-effect/structure ML developer","fit_score":3,"primary_use_case":"LFC3D as features/labels; benchmark spatial aggregation","would_use_instead":"Reimplement neighborhood aggregation","adoption_likelihood":"Medium"},
  {"id":"A9","name":"Grad student / postdoc operator","sector":"academic","role":"Hands-on pipeline operator","fit_score":5,"primary_use_case":"Run Colab pipeline; generate thesis/paper 3D figures","would_use_instead":"Ask a bioinformatician; hack scripts","adoption_likelihood":"High"},
  {"id":"A10","name":"Cross-species / evolutionary FG lab","sector":"academic","role":"Comparative functional-genomics PI","fit_score":3,"primary_use_case":"Cross-species/isoform conservation of functional hotspots","would_use_instead":"Manual alignment + custom mapping","adoption_likelihood":"Medium"},
  {"id":"A11","name":"Drug-target-discovery academic","sector":"academic","role":"Chemical-biology / target-validation PI","fit_score":4,"primary_use_case":"Map screen hits to drug/allosteric pockets on liganded PDB","would_use_instead":"Fpocket/SiteMap + manual overlay","adoption_likelihood":"Medium"},
  {"id":"I1","name":"Pharma target-discovery / functional-genomics group","sector":"industry","role":"Industrial FG / target-ID team","fit_score":4,"primary_use_case":"Prioritize target residues; consolidate internal screens; allosteric handles","would_use_instead":"Internal comp platforms, CRO analysis","adoption_likelihood":"Medium"},
  {"id":"I2","name":"Base-editing therapeutics company","sector":"industry","role":"BE therapeutics R&D (Beam/Verve/Prime/Intellia)","fit_score":4,"primary_use_case":"Interpret internal BE tiling of a target; edit-site mechanism","would_use_instead":"Internal validated platforms, academic collaboration","adoption_likelihood":"Medium"},
  {"id":"I3","name":"CRO offering screening services","sector":"industry","role":"Contract screening/analysis provider","fit_score":4,"primary_use_case":"Premium 3D-hotspot report tier on top of screens","would_use_instead":"Build own; deliver raw hit lists","adoption_likelihood":"Medium"},
  {"id":"I4","name":"Ag-bio / crop & livestock editing company","sector":"industry","role":"Agricultural base-editing R&D","fit_score":2,"primary_use_case":"Cross-species hotspot transfer; residues to edit/avoid","would_use_instead":"Bespoke scripts (few off-the-shelf options)","adoption_likelihood":"Low"},
  {"id":"I5","name":"Synthetic biology / enzyme engineering firm","sector":"industry","role":"Industrial protein optimization","fit_score":2,"primary_use_case":"Critical-residue maps to guide library design","would_use_instead":"In-house DMS pipelines, Rosetta, ML","adoption_likelihood":"Low"},
  {"id":"I6","name":"Antibody / protein-engineering firm","sector":"industry","role":"Antibody/binder optimization","fit_score":2,"primary_use_case":"Interface/paratope functional mapping","would_use_instead":"Specialized antibody platforms","adoption_likelihood":"Low"},
  {"id":"I7","name":"Clinical dx / variant-interpretation company","sector":"clinical","role":"VUS classification / functional-evidence provider","fit_score":2,"primary_use_case":"Spatial functional evidence in variant reports","would_use_instead":"AlphaMissense/ESM, direct MAVE scores, in-house","adoption_likelihood":"Low"},
  {"id":"I8","name":"AI-bio / techbio startup","sector":"industry","role":"ML variant-effect/protein-function developer","fit_score":3,"primary_use_case":"LFC3D features/benchmarks for models","would_use_instead":"Reimplement; use raw screens","adoption_likelihood":"Low"},
  {"id":"I9","name":"Tool / SaaS / portal vendor","sector":"industry","role":"Structure-viewer / genomics-portal vendor","fit_score":3,"primary_use_case":"Integrate BE-screen->3D-hotspot as a product feature","would_use_instead":"Build own engine","adoption_likelihood":"Medium"}
]
```

---
## Sources
- CRISPR screening market: https://www.towardshealthcare.com/insights/crispr-market-sizing ; https://www.globenewswire.com/news-release/2026/04/20/3277327/28124/en/Arrayed-CRISPR-Screening-Market-Presents-Lucrative-Opportunities-Through-2026-2030-Total-Revenue-to-Grow-by-1-64-Billion-at-17-7-CAGR.html
- Functional-genomics segment: https://www.grandviewresearch.com/horizon/statistics/crispr-and-cas-genes-market/biomedical/functional-genomics/global
- Gene-editing company count & funding: https://www.samps.org/blog/gene-editing-landscape-2025 ; https://www.nature.com/articles/d43747-025-00019-z ; https://crisprmedicinenews.com/companies/
- Beam / Verve / Prime clinical & financing: https://pharmaphorum.com/news/beam-launches-500m-financing-base-editing-trial-data ; https://www.biospace.com/drug-development/safer-crispr-base-editing-breaks-through-in-the-clinic-as-beam-verve-advance ; https://www.genengnews.com/topics/genome-editing/base-editors-and-prime-editors-begin-to-realize-their-clinical-promise/
- Base-editor screens growth: https://www.nature.com/articles/s41588-025-02366-0 ; https://www.sciencedirect.com/science/article/pii/S1097276523004318 ; https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12230738/ (CRISPR-BEasy)
- AlphaFold adoption: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10767828/ ; https://www.pnas.org/doi/10.1073/pnas.2315002121 ; https://arxiv.org/pdf/2403.02124
- SGE / VUS: https://www.medrxiv.org/content/10.1101/2025.11.03.25339440v2.full ; https://www.biorxiv.org/content/10.1101/294520.full.pdf
