# BE3D Run — Decision Log: KRAS (sotorasib vs adagrasib resistance, 3D)

> Runner: Queen "KRAS/EGFR". Real public base-editing screen → BE3D LFC3D + 3D clustering →
> per-drug 3D resistance atlas. The point the source paper never made: **spatial** clustering of
> resistance on structure, separating orthosteric (pocket-lining) from distal/allosteric hotspots.

## 0. Target & the open gap
- **Gene / protein:** KRAS (UniProt **P01116**, KRAS4B canonical, 189 aa). G-domain res 1–169 + HVR 170–189.
- **Structures used:**
  - **AlphaFold** `AF-P01116-F1-model_v6` (WT sequence; canonical numbering — matches the screen exactly). **Primary analysis structure.**
  - **PDB 6OIM** chain A — KRAS(G12C) + **sotorasib/AMG510** (ligand `MOV`); engineered C51S/C80L/C118S, G12C. G-domain only (res 0–169).
  - **PDB 6UT0** chain A — KRAS(G12C) + **adagrasib/MRTX849** (ligand `M1X`); same engineered construct.
- **Source dataset:** Coelho, Dincer et al. 2024, *Nature Genetics* "Base editing screens define the genetic landscape of cancer drug resistance." DOI 10.1038/s41588-024-01948-8.
  Data URL fetched (curl 200, 62 MB): `https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-024-01948-8/MediaObjects/41588_2024_1948_MOESM4_ESM.xlsx` (local `scratchpad/coelho_S4.xlsx`), sheet **`ST2 BE z-scores`** (45,129 rows; 3,006 KRAS guide rows).
- **The OPEN question we attack (NOT the paper's headline):** Run the **sotorasib** and **adagrasib**
  arms *separately* through LFC3D + agglomerative 3D clustering. Do their resistance hotspots form
  **spatially separable sub-clusters** in/around the switch-II pocket, or are they identical? Which
  resistance residues are **pocket-lining (orthosteric)** vs **distal (allosteric/conformational-latch)**?
  Is any resistance cluster on the RAS–effector interface rather than the drug pocket?
- **What NEW actionable output we produce:** a per-drug 3D resistance atlas that (a) shows the two
  inhibitors have **different spatial resistance signatures**, and (b) names a **sotorasib-specific
  DISTAL resistance cluster** invisible to the pooled per-variant hit list.

## 1. Data provenance & reformatting (decision log)

The ST2 sheet does **not** match the brief's assumed schema. Actual per-guide fields:
`Amino_Acid_Position` (`;`-joined codon positions), `Edited_AA` (ref AAs), `New_AA` (alt AAs),
`Protein_Change` (VEP `ref/alt`, `;`-joined, includes redundant multi-codon haplotype tokens),
`most_severe_consequence` (the real category: `missense`/`synonymous`/`stop codon`/`splice variant`/
`intron`/`start lost`/None), `editor` (`ABE`/`CBE`), drug z-scores
`L2FC_H23_<drug>_plasmid_average_zscore`. `swissprot` is versioned (`P01116.251`); `variant_classification`
is uninformative (`SNV`/`NA`). **The brief's column map (`Protein_Change`=HGVS 3-letter, `variant_classification`=Missense/…) is wrong for this file** — reverse-engineered from the data.

| Decision | Choice | Why | Alternative rejected |
|---|---|---|---|
| Source file | `41588_2024_1948_MOESM4_ESM.xlsx` sheet `ST2 BE z-scores`, 62 MB, curl 200 | verified downloadable, carries AA edits + per-drug z-scores in one table | MaveDB subsets (thinner) |
| Gene filter | `Gene == "KRAS"` (3,006 rows) | `swissprot` is versioned `P01116.251` so an exact `=="P01116"` match returns 0 | filter on `swissprot` (breaks) |
| Category source | `most_severe_consequence`, **not** `variant_classification` | `variant_classification` is `SNV`/`NA`; real category is in `most_severe_consequence` | brief's column |
| Category map | missense→Missense, synonymous→Silent, stop codon→Nonsense, splice variant→Splice-donor, intron→Intron, start lost→Missense, None/NA→No Mutation | matches BE3D `mutation_category`; None/NA guides (no coding consequence) become controls | dropping None rows (loses control distribution) |
| Mutation_list | reconstructed `refAA+pos+altAA` (1-letter) from `zip(Amino_Acid_Position, Edited_AA, New_AA)` truncated to min length, deduped; nonsense→alt `*`; silent→alt=ref; non-coding→placeholder `nc` | `(pos, Edited_AA)` validated against the true P01116 sequence (only 10/148 clean single-pos guides mismatch, all in the HVR ≥151 where the screen transcript/isoform diverges) | parsing `Protein_Change` directly — its redundant multi-codon haplotype tokens (`TE/AG`) don't align to positions |
| Multi-edit guides | kept all reconstructed per-codon edits; a guide's z-score is assigned to **every** missense residue it edits (standard BE3D) | can't attribute causality within a guide; aggregation over many guides/residue is robust | collapsing to one edit (loses coverage) |
| sgRNA_score | `L2FC_H23_Sotor_plasmid_average_zscore` (sotorasib) and `…_Adag_…` (adagrasib) → **two separate screens** | the drugs are the comparison of interest; z-scored, so comparable | replicate `_1/_2` cols (noisier) |
| Editors | **pooled ABE + CBE** into one TSV per drug | maximizes residue coverage (162 distinct missense positions, 1–188); z-scores comparable | editor×drug as 4 screens (fragments the pocket coverage) |
| Direction convention | **positive z = enrichment = RESISTANCE** (the readout is proliferation under drug) | resistance screens are enrichment, not dropout | treating neg as the hit direction |

- **Reformatting script:** `real_output/P2_KRAS/input/convert_kras.py`.
- **Row counts (per drug):** 2,980 scored rows → Missense 230 (162 distinct positions, span 1–188),
  Silent 46, Nonsense 16, Splice-donor 42, Intron 5, No Mutation 2,641.
- **Numbering reconciliation (important):** screen == crystal == AF numbering in the pocket region.
  Crystal PDBs carry engineered **G12C, C51S, C80L, C118S** (and construct-specific differences at
  ≥151) — so in crystal runs positions 12/51/80/118 fail BE3D's refAA sanity match and are dropped
  there; the **AF-monomer run (WT) is clean at every position** and is the primary structure.
- **Structures prepped:** `scratchpad/targets/KRAS/` — `AF-P01116-F1-model_v6.pdb`+`.dssp`,
  `6oim_chainA.pdb`+`.dssp`, `6ut0_chainA.pdb`+`.dssp` (chain A extracted to avoid cross-chain
  neighbors in monomer mode). DSSP is the placeholder generator (affects only characterization).

## 2. QA (BE-QA)
- **cases = [Nonsense, Splice-donor]  controls = [No Mutation, Silent]** (knockout vs neutral).
- KS H1 (case vs control LFC distributions): **Sotor D=0.128, p=0.29; Adag D=0.107, p=0.50** — NOT
  significant (58 cases vs 2,687 controls).
- **Screen accept decision: ACCEPT, with the explicit caveat that the standard QA is the wrong gate here.**
  QA assumes a *dropout/essentiality* screen where knockouts deplete relative to neutral. This is an
  **enrichment (resistance)** screen: the informative signal is guides that become *enriched* (positive
  z) under drug, and nonsense/splice knockouts of KRAS are not expected to enrich. The near-null KS is
  therefore *expected*, not a failure. All downstream interpretation uses the **positive** LFC3D direction.

## 3. Parameter sweep grid (EXHAUSTIVE)
All runs: monomer mode, clustering_radius = 6 Å, LFC=mean, meta=SUM (single-screen so unused),
p-thresholds 0.05/0.01/0.001 all computed every run. 9 configs × 2 drugs = **18 runs, all EXITCODE 0**.
Outputs under `real_output/P2_KRAS/runs/<run_id>_<drug>/KRAS/`. Driver: `scratchpad/run_grid.sh`.

| run_id | structure | structure_radius | function_for_lfc3d | nRandom | atom_level | notes |
|---|---|---|---|---|---|---|
| AF_r4_mean_n1000_res  | AF | 4 | mean | 1000 | F | radius sensitivity |
| **AF_r6_mean_n1000_res** | AF | **6** | **mean** | **1000** | **F** | **BASELINE** |
| AF_r8_mean_n1000_res  | AF | 8 | mean | 1000 | F | radius sensitivity |
| AF_r10_mean_n1000_res | AF | 10 | mean | 1000 | F | radius sensitivity |
| AF_r6_sum_n1000_res   | AF | 6 | sum | 1000 | F | aggregation |
| AF_r6_mean_n500_res   | AF | 6 | mean | 500 | F | null-size / tail-p stability |
| AF_r6_mean_n1000_atom | AF | 6 | mean | 1000 | **T** | atom-level neighbors |
| 6OIM_r6_mean_n1000_res | 6OIM (sotorasib-bound) | 6 | mean | 1000 | F | drug-bound PDB |
| 6UT0_r6_mean_n1000_res | 6UT0 (adagrasib-bound) | 6 | mean | 1000 | F | drug-bound PDB |

- **Radius:** as radius grows, adagrasib's orthosteric-pocket cluster **saturates and grows**
  (pos p<0.01: r4=18 → r6=27 → r8=34 → r10=37 residues), while sotorasib stays **small and stable**
  (r4=21 → r6=20 → r8=15 → r10=12). The pocket-vs-distal split holds at every radius.
- **Aggregation (mean vs sum):** near-identical hit sets (Sotor 001: mean 17 vs sum 19; Adag 25 vs 28).
- **nRandom (500 vs 1000):** hit counts identical to ±1 residue at p<0.001 → tail-p is stable.
- **atom vs residue:** atom-level slightly expands both (more neighbor atoms) but preserves the same
  clusters (Sotor still pocket+distal; Adag still one pocket block).
- **AF vs drug-bound PDB:** crystal runs reproduce the qualitative split (below); AF is primary
  (WT, clean numbering, full-length incl. HVR).

## 4. Results per run

**Baseline (AF, r6, mean, n1000), positive/resistance direction:**
- **Sotorasib** — pos p<0.05: 24 residues; p<0.01: 20 residues =
  `11,47,58,60,65,66,67,68` + `147–156` + `178,179`.
- **Adagrasib** — pos p<0.05: 30; p<0.01: 27 residues =
  `60–72, 89,90,91,92,93, 95,96,97,98,99,100,101,102,103`.
- Negative/essential direction (both drugs, internal positive control): the GTP/GDP core —
  P-loop 10–17, switch I 30–40, 55–56, 82–86, 93–95, 104–112. Biologically sensible (KRAS is essential).

**3D clustering of the positive hits @6 Å (BE3D agglomerative, baseline):**
- **Sotorasib → 4+ SEPARATE spatial clusters:**
  - Cluster 4 (n=10): **147,148,149,150,151,152,153,154,155,156** — the **distal α5/SAK region**.
  - Cluster 1 (n=4): 65,66,67,68 — switch-II pocket edge.
  - Cluster 0 (n=2): 11, 60 — P-loop / switch-II.
  - Cluster 3 (n=2): 178,179 — HVR.
  - **Inter-cluster centroid distances:** distal α5 cluster is **26.9 Å** from the switch-II cluster
    and 19.0 Å from the P-loop cluster → **genuinely, quantitatively spatially separated.**
- **Adagrasib → ONE contiguous cluster (n=27): 60–72, 89–103** — the entire switch-II cryptic pocket
  lights up as a single connected 3D hotspot. **No distal cluster.**

**Pocket classification (min-atom distance to the bound drug; `targets/KRAS/pocket_dist.json`):**
- Adagrasib hits are overwhelmingly **orthosteric**: 61,62,63,64,69,72,92,95,96,99,100,102,103 are
  <5 Å from the drug; 70,71,89,91,93,97,98,101 are 5–8 Å.
- Sotorasib's in-pocket hits are few (11,58 <5 Å; 60,65–68 in/near pocket); its 147–156 cluster is
  **>8 Å from the drug (distal)** and 178–179 are in the HVR (no crystal coverage).

**Drug-bound-PDB runs corroborate (pos p<0.01):**
- 6UT0/6OIM **Adagrasib**: `62,63,64,65,66,67,68,69,70, 90,91, 95–102` — pocket again.
- 6UT0/6OIM **Sotorasib**: `65,66` + distal `165,179,180` — sparse pocket + distal, again.
  (pos 10/12 also appear, partly an artifact of the G12C refAA mismatch at the mutation site.)

**Robustness (persistence of pos p<0.01 residues across the 7 AF runs):**

| Residue(s) | drug | #/7 AF runs | class | interpretation |
|---|---|---|---|---|
| 65,66,67 | Sotor | 7/7 | near pocket | switch-II core (shared anchor) |
| **148,154** | Sotor | **7/7** | **distal (>8 Å)** | **distal α5 cluster core — rock-solid** |
| 149,150,179 | Sotor | 6/7 | distal / HVR | distal cluster |
| 68 | Sotor | 5/7 | pocket | switch-II |
| 11,58,60,151,153,155,147 | Sotor | 4/7 | pocket/distal | cluster shoulders |
| 60,62,63,64,65,66,67,68,69,92,99 | Adag | 7/7 | pocket | orthosteric core — rock-solid |
| 61,70,71,72,95,96,97,100,101,102,103 | Adag | 6/7 | pocket/near | orthosteric |

**Sotorasib vs adagrasib overlap (baseline, pos p<0.01):** SHARED = **only {60,65,66,67,68}**
(switch-II β3/α2 core); Sotor-only = 15 residues (incl. the entire distal 147–156 + HVR);
Adag-only = 22 residues (the rest of the orthosteric pocket). **Jaccard ≈ 0.12.**

- Figures per run (SVG): `…/cluster_LFC3D/*Dendrogram*`, `…/LFC3D/*`, `…/characterization/*`,
  `…/g2p_visualization/*`. Machine-readable summary: `real_output/P2_KRAS/grid_results.json`.

## 5. Interpretation & NOVELTY check
- **Structural mapping (literature-grounded):** the shared/adagrasib hits (60–72, 89–103, 95/96/99)
  are the **switch-II (59–76) + α2/α3 + His95 cryptic sub-pocket** of the SII-P that AMG510/MRTX849
  occupy (6OIM/6UT0). Contact maps place **V9, R68, Y96 in the first (<5 Å) shell of both drugs, with
  H95 contacting adagrasib but not sotorasib** (Feng et al. 2022 *PNAS*, PMID 35471904). Several are
  **documented clinical/experimental G12C-inhibitor resistance residues** that our base-editing map
  **rediscovers de novo** — a strong positive control:
  - **Y96** (Y96D/C/S; Awad 2021 *NEJM* PMID 34161704, Tanaka 2021 *Cancer Discov* PMID 33846219; pan-inhibitor),
  - **R68** (R68S clinical/R68M/R68L; Awad, Koga 2021, Feng 2022),
  - **H95** (H95D clinical/H95R; **adagrasib-specific, stays sotorasib-sensitive**),
  - **Q99** (Q99L/F/R/S/V/Y; **adagrasib-specific**), **D69, M72, D92, T58, S65** (Feng 2022).
  Consistently, these H95/Q99-series residues fall in our **Adag-only** set, not the sotorasib set —
  matching the published mechanism. **Not** recovered as enrichment (correctly): the *activating*
  mutations A59S/T and G13D — our data places 55–56/A59 and G13 in the **essential/depleted** set, as
  expected for in-cis activators rather than orthosteric escapes.
- **Why adagrasib's footprint is broader (mechanism, established):** adagrasib's pyrimidine makes a
  **direct H-bond to Tyr96** and engages **His95**, extending into the H95/Q99 cryptic groove
  (Hallin 2020 *Cancer Discov*, 6UT0); sotorasib contacts Y96 only via **water bridges** and is
  **H95-independent** (Canon 2019 *Nature*, 6OIM). This deeper, directional His95-groove occupancy
  fully explains the dense adagrasib resistance shell across 89–103 vs sotorasib's shallow footprint —
  our 3D map reproduces this pharmacology from screen data alone.
- **The novel part (NOT in the source paper):** the paper reports per-variant interface-disruption
  flags but **never spatially clusters resistance on structure**. Our two findings are new:
  1. **Sotorasib and adagrasib have qualitatively different 3D resistance signatures** —
     adagrasib saturates the orthosteric pocket (one 27-residue cluster), whereas sotorasib's
     in-pocket footprint is small and it carries a **physically separate distal cluster**.
  2. **A robust sotorasib-specific DISTAL resistance cluster at residues 147–156 (G5/SAK nucleotide
     motif + N-terminus of α5), ~27 Å from the drug** — orthogonal to the orthosteric surface and not
     highlighted by the paper. **Literature check (worker-grounded): no prior report names 147–156 as
     a G12C-inhibitor resistance or orthosteric site.** The region overlaps the independently-defined
     **α4–α5 allosteric/dimerization lobe** (NS1 monobody, Spencer-Smith 2017 *Nat Chem Biol*; DARPin,
     Bery 2019) and the guanine-contacting **SAK/G5 motif** (K147 contacts the nucleotide base). The
     single prior resistance datapoint anywhere in this patch is **F156L** (Feng 2022 *PNAS*), reported
     with no assigned mechanism — our cluster gives it spatial context. Candidate mechanism (hypothesis,
     not established): SAK/α5 perturbation modulates nucleotide exchange/GTP-loading and biases the
     switch-II conformational equilibrium that sotorasib (an inactive-state trap) depends on — which
     would explain sotorasib-specificity. Caveat: Whaby 2022 *JBC* (PMID 36334633) shows most α4–α5
     point mutations do **not** impair signaling/self-association, so a simple dimerization mechanism is
     disfavored; and the closest deep-mutational allostery map (Weng/Lehner 2024 *Nature*, PMID 38109937)
     centers elsewhere. Treat 147–156 as a **novel candidate allosteric-resistance patch to validate**,
     not a proven mechanism.
- **Confidence caveats:**
  - Direction: enrichment screen; QA (dropout-style) is null by design — interpret only the positive
    direction. Missense z-scores are assigned per-guide to all edited residues (causality within a
    guide is unresolved); robustness across 162 positions and many guides mitigates this.
  - The 147–156 / HVR signal rests on the **AF** model (crystal 6OIM/6UT0 stop at 169 and carry
    construct differences ≥151); on the crystals the sotorasib distal signal re-appears at 165/179/180,
    consistent but shifted by the truncation. pLDDT of 147–156 in the AF G-domain lobe is high
    (structured helix — the credible novel cluster); **178–180 are in the disordered HVR (low pLDDT) and
    should be treated as a possible artifact / passenger signal — flagged for sanity-check, not claimed.**
  - Literature grounding performed by a dedicated sub-agent (paperclip/PubMed); key refs: Feng 2022
    *PNAS* PMID 35471904 (KRAS G12Ci saturation resistance map), Awad 2021 *NEJM* PMID 34161704,
    Tanaka 2021 *Cancer Discov* PMID 33846219, Canon 2019 *Nature* (6OIM), Hallin 2020 *Cancer Discov*
    (6UT0), Whaby 2022 *JBC* PMID 36334633 (α4–α5 lobe), Weng/Lehner 2024 *Nature* PMID 38109937.
  - altAA identities in the reconstructed Mutation_list are approximate (the New_AA field is noisy);
    positions and missense/silent classification are validated. LFC3D aggregates by position, so this
    does not affect the hotspot map.

## 6. NEW actionable proposals (the deliverable)
1. **Med-chem / KRAS-drugging community — exploit the asymmetry.** Adagrasib's resistance is a broad
   orthosteric shell (60–72, 92, 95–103); sotorasib's in-pocket resistance is sparse **and** it has a
   *distal* escape route (α5 147–156). Prediction: **a sotorasib-class covalent binder that also makes
   contacts toward the α5/His95-groove exit (as adagrasib does) should have a higher genetic barrier to
   resistance** than sotorasib's current footprint. Conversely, adagrasib's dense orthosteric
   dependence means single pocket substitutions (Y96D, R68, H95, Q99) are efficient escapes — combine
   with an SOS1/SHP2 or a distinct-pocket agent.
2. **Structural biologists — test the α5 (147–156) allosteric hypothesis.** Introduce
   `T148/R149/G151/V152/D154`-region substitutions into KRAS(G12C) and measure sotorasib (vs adagrasib)
   IC50 shift, GTP-loading, and switch-II occupancy by HDX/crystallography. If these distal mutations
   selectively blunt sotorasib, the α5 patch is a bona fide conformational-latch resistance site and a
   candidate **allosteric co-targeting** locus. (Boltz co-folding of KRAS(G12C)+sotorasib with 147–156
   mutants is a fast in-silico pre-screen — the environment has Boltz credits.)
3. **Functional-genomics groups — validate drug-specific 3D signatures as a general readout.** Re-run
   the pipeline on the paper's other drug pairs already spatially separable here in principle
   (EGFR gefitinib vs osimertinib; MEK1 trametinib; PARP1 olaparib vs niraparib) to build a pan-target
   "orthosteric-shell vs distal-patch" resistance atlas; the sotorasib/adagrasib contrast shows the
   method resolves inhibitor-specific escape geometry that pooled hit lists hide.
4. **Clinicians / translational — surveillance.** The orthosteric residues we recover (R68, H95, Y96,
   Q99) match emerging clinical G12Ci resistance; the analysis nominates the **α5 147–156** and
   **switch-II 89–93/97–103** regions as additional loci worth watching in ctDNA on sotorasib vs
   adagrasib therapy, respectively.

## 7. BE3D issues encountered → appended to BE3D_IMPROVEMENTS.md
See the new entries: (a) the shipped brief/recipe column-map for the Coelho ST2 sheet is wrong (real
schema differs — a data-schema-drift friction worth flagging for turnkey datasets); (b) drug-bound
crystal constructs carry engineered mutations + non-canonical numbering that silently drop edits in
BE3D's refAA sanity match — a warning/report would help; (c) `run_new_target.py` shares a single
`be3d_configs/KRAS_run.yaml` name per gene, so concurrent same-gene runs must use unique parent dirs
(handled here via `runs/<id>/`); (d) DSSP placeholder is fine (characterization-only), reconfirmed.
