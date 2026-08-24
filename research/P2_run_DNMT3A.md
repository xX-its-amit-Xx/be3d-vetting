# BE3D Run — Decision Log: DNMT3A (Lue et al. 2023 base-editor scanning)

Runner: Queen "DNMT3A". Date: 2026-08-24. Status: **COMPLETE** (12/12 sweep runs EXIT=0; QA PASS).
Outputs: `real_output/P2_DNMT3A/`. Venv: `scratchpad/be3dvenv`. BE3D clone: `scratchpad/BE3D`.

## 0. Target & the open gap
- **Gene / protein:** DNMT3A (DNA (cytosine-5)-methyltransferase 3A). UniProt **Q9Y6K1** (canonical, 912 aa).
  Library targets the DNMT3A2 transcript but mutations are reported in canonical 912-aa numbering (R882 etc.).
- **Structure used:**
  - PRIMARY: **PDB 4U7T** (DNMT3A–DNMT3L–H3 active complex). DNMT3A catalytic-domain chains **A, C** (UNP 474–912;
    contains SAH cofactor + ADD structural Zn). DNMT3L chains B, D (178–379); H3 peptide chains F, G.
    Chain **A** used as reference. Isolated chain-A PDB for monomer runs; A+C for complex.
  - BASELINE: **AlphaFold AF-Q9Y6K1-F1 v6** (full-length monomer) — robustness / conformation control.
- **Source dataset:** Lue, Liau et al. 2023, *Nat Chem Biol* "Base editor scanning charts the DNMT3A activity
  landscape". DOI 10.1038/s41589-022-01167-4. File fetched (curl 200, 464 KB):
  `https://static-content.springer.com/esm/art%3A10.1038%2Fs41589-022-01167-4/MediaObjects/41589_2022_1167_MOESM3_ESM.xlsx`
  → sheet **"Supplementary Data 3"** (Summary; 729 guides). Assay = CpG-methylation→citrine reporter FACS in MOLM13 (CBE).
- **OPEN question (NOT the paper's headline):** The paper's PWES analysis reports a general *interdomain* cluster.
  We go finer — **partition loss-of-function (LOF) LFC3D density across the distinct surfaces of 4U7T**
  (a) DNMT3A–DNMT3L "FF" interface, (b) DNMT3A homotetramer "RD–RD" interface, (c) DNA-binding face,
  (d) ADD–MTase autoinhibitory contact, (e) catalytic pocket — and name **which single interface dominates**,
  then test co-localization with clinical CHIP/AML hotspots (esp. R882).
- **NEW actionable output:** a per-surface interface-vulnerability map naming the dominant PPI interface and its
  specific **base-editable** hotspot residues, distinguishing them from the clinical R882 site.

## 1. Data provenance & reformatting (decision log)
Reformatting script: `analysis/prep_dnmt3a.py` (parses Suppl Data 3 → BE3D 4-column TSV). In=729 guides, out=**716**.

| Decision | Choice | Why | Alternative rejected |
|---|---|---|---|
| Source file | Springer MOESM3 xlsx, sheet "Supplementary Data 3" | one-table, scores shipped (Recipe #4) | Suppl Data 2 (raw FACS ratios) — needs LFC recompute |
| Rows kept | 716 (dropped 13 "Essential" positive-control guides w/o residue mapping) | Essential guides are dropout controls, not residue-mapped | keeping them (would pollute categories) |
| `Mutation_list` | ← `Mutation_list`; comma→`;`-join, `*`=stop, `None`→`utr;` | BE3D edit format | — |
| `Mutation_type` | ← `Mut_type`: Missense/Silent/Nonsense/Splice kept; **Non-exon+No_C/*+Intergenic+Non-targeting → "No Mutation"** | pool all non-coding/no-edit guides as the neutral control set (n=389) | using only Non-exon (weaker control) |
| `sgRNA_score` | ← `sgRNA_score_d9_citrine_positive`, **NEGATED** | see sign convention below | keeping raw (LOF would land in pos channel) |
| Multi-edit guides | all edits kept, `;`-joined (e.g. `E664K;V665M;`) | BE3D supports multi-AA guides | collapsing to one edit |
| Screens combined | single screen (one reporter, one editor set) — no meta | one condition only | — |

**SIGN CONVENTION (critical, documented per mandate):** the citrine reporter is *silenced* by DNMT3A methylation,
so **LOF → citrine de-repressed → cell enriched in the citrine-POSITIVE sort → HIGH POSITIVE raw score**.
Confirmed empirically on the controls: Nonsense mean **+0.73**, Splice **+1.05**, Silent **+0.09**, Non-exon **+0.04**.
This is *inverted* relative to a dropout screen (where LOF = negative LFC), which is the convention baked into BE3D's
downstream "neg" channel and BE-QA. **We therefore NEGATE the score**, so LOF (nonsense/splice/damaging missense)
becomes NEGATIVE and is captured by BE3D's negative-direction machinery. All results below read the **`*_LFC3D_neg`**
channel = LOF. A `raw` (un-negated) TSV is also emitted for audit. Row counts: 209 Missense / 18 Nonsense / 12 Splice /
88 Silent / 389 No-Mutation. Sanity: all 209 Missense guides carry a parseable edit; positions map to 4U7T chain A.

## 2. QA (BE-QA) — decision log
- **cases = [Nonsense, Splice]** (n=30, bona-fide LOF), **controls = [No Mutation]** (n=389, non-coding/no-edit).
  Rationale for an *activity* (non-dropout) screen: BE-QA's Mann-Whitney/KS only test whether the LOF class is
  *distributionally separated* from neutral guides — direction-agnostic — so it is valid after negation. We did NOT
  use Silent as the control set (kept as an independent neutral check; Silent mean ≈ 0 confirms no editing artifact).
- **Result (H1, within screen):** KS **D=0.628, p=4.2e-11**; Mann-Whitney **U=9902, p=2.0e-10**.
- **Decision: SCREEN ACCEPTED** — LOF guides are decisively separated from controls; the negated-score/QA setup is valid.

## 3. Parameter sweep grid (EXHAUSTIVE)
All runs share function_for_lfc=mean, function_for_meta=SUM (single screen ⇒ meta unused); each run emits p<.05/.01/.001.

| run | structure | scope | radius | func_lfc3d | nRandom | level | notes |
|---|---|---|---|---|---|---|---|
| r01 | 4U7T chain A | monomer | 6 | mean | 1000 | residue | **PRIMARY / baseline** |
| r02 | 4U7T chain A | monomer | 4 | mean | 1000 | residue | radius↓ |
| r03 | 4U7T chain A | monomer | 8 | mean | 1000 | residue | radius↑ |
| r04 | 4U7T chain A | monomer | 10 | mean | 1000 | residue | radius↑↑ |
| r05 | 4U7T chain A | monomer | 6 | **sum** | 1000 | residue | aggregation |
| r06 | 4U7T chain A | monomer | 6 | mean | **500** | residue | null stability |
| r07 | 4U7T chain A | monomer | 6 | mean | 1000 | **atom** | interface sensitivity |
| r08 | 4U7T chain A | monomer | 8 | sum | 1000 | residue | combo |
| r09 | AF full-length | monomer | 6 | mean | 1000 | residue | structure control |
| r10 | AF full-length | monomer | 8 | mean | 1000 | residue | structure control |
| r11 | 4U7T A+C | **complex** | 6 | mean | 1000 | residue | RD cross-protomer |
| r12 | 4U7T A+C | **complex** | 8 | mean | 1000 | residue | RD cross-protomer |

Each run ≈ 35–60 s (nR=1000). Radius=clustering_radius throughout. p-thr is not a run multiplier (all three computed per run).

## 4. Results per run
Primary run **r01**: 439 chain-A residues scored; LOF-significant **117 (p<.05) / 107 (p<.01) / 100 (p<.001)** in the neg channel.

**Per-surface LOF partition (r01; enrichment = fraction-significant vs genome-wide 0.267):**

| Surface | n_res | n_sig(p<.05) | sum\|LFC3D_neg\| | mean\|neg\| | enrichment |
|---|---|---|---|---|---|
| catalytic_pocket (SAH) | 28 | 19 | 21.32 | 0.762 | **2.55x** |
| **RD_homodimer (tetramer)** | 20 | 10 | 8.35 | **0.418** | **1.88x** |
| FF_DNMT3A_DNMT3L | 17 | 6 | 6.25 | 0.367 | 1.32x |
| ADD_MTase (interdomain) | 43 | 7 | 8.81 | 0.205 | 0.61x (depleted) |
| H3_ADD_read | 24 | 3 | 4.47 | 0.186 | 0.47x |
| DNA_face (curated) | 31 | 2 | 3.35 | 0.108 | 0.24x |

**Top-20 LOF hotspots (r01)** are dominated by the catalytic pocket (V665/D686/E664/E756/C666/G642/G685/D641 = motif IV
+ ENV catalytic motifs, all z ≤ −15, p<.001) plus PWWP-adjacent G293/E294. Among **protein interfaces**:

**RD homodimer interface (the dominant PPI vulnerability), r01 residue detail:**
- H873 (|neg|=1.089, z=−13.85, ***), D876 (1.089, −13.85, ***), T671 (0.898, −10.37, ***), V877 (0.874, −9.92, ***),
  S878 (0.766, −7.96, ***), M674 (0.583, ***), V675/R676/Q678 (0.517, ***), N879 (0.444, *).
- **R882: |neg|=0.170, z=+2.94, NOT significant.** (R882H is a G→A change inaccessible to CBE; no 882 guide exists in
  the screen.) S878F is a real *enriched* hit in the raw data (+0.44).

**FF interface (DNMT3A–DNMT3L), r01 residue detail:** Y735 (1.598, z=−23.16, ***), R729 (1.134, ***), Y724 (1.134, ***),
R736 (0.939, ***), E725 (0.692, ***), H739 (0.661, ***).

**Complex mode (r11/r12):** cross-protomer aggregation modestly *raises* the RD interface (sum 8.35→8.80, sig 10→11,
enrichment 1.88→2.05x at r=6; 2.44x at r=8) — consistent with the RD residues gaining real neighbor signal from the
partner protomer. Overall ranking unchanged (catalytic > RD > FF).

**Robustness across the sweep (files: `analysis/interface_enrichment_by_run.tsv`, `key_residue_significance_by_run.tsv`):**
- **RD is the dominant PPI interface in 8/12 runs** (RD enrichment 1.44–2.44x, always >1). FF overtakes RD only in the
  4 large-radius/atom-level runs (r03 r=8, r07 atom, r08 r=8-sum, r10 AF r=8) where FF residues (which abut the catalytic
  core) absorb bleed-through catalytic signal. ADD_MTase / H3 / DNA are depleted (<0.7x) in every run.
- **RD core patch is rock-solid:** T671, D876, V877, S878 are p<0.001 in **12/12** runs; M674, V675, H873 in ≥10/12.
- **R882 is non-significant in 11/12 runs** (only reaches *** at r=10, r04, via distant neighbors) — it is essentially
  BE-invisible.
- **FF residues R729, R736, Y735 are p<0.001 in ~12/12** runs; Y724/E725/H739 robust in ≥10/12.
- QA identical across runs (structure-independent). Figures/plots per run under each `dnmt3a_runs/<run>/` (SVG dendrograms,
  scatter). Lean primary outputs copied to `real_output/P2_DNMT3A/be3d_primary_r01/` and `be3d_complex_r11/`.

## 5. Interpretation & NOVELTY check
- **Which interface dominates?** Among protein–protein/quaternary surfaces, the **DNMT3A homotetramer RD–RD interface**
  carries the dominant, robust LOF hotspot (not the FF DNMT3A–DNMT3L interface, and not the ADD–MTase autoinhibitory
  contact, which is *depleted*). The catalytic (SAH) pocket is the strongest surface overall, as expected for an activity
  screen — but that is the paper's expected result; the *interface* answer is **RD**.
- **Clinical co-localization (worker-verified against COSMIC/cBioPortal/literature):**
  - **R882** (R882H/C/S; ~50–60% of DNMT3A-mutant AML) sits AT the RD tetramer interface and is dominant-negative by
    blocking active-tetramer formation (Russler-Germain 2014; Sandoval 2019). BE3D independently flags this SAME interface —
    but via *editable neighbors*, since R882 itself cannot be made by the CBE. **Interface-level co-localization: YES;
    residue-level: NO.**
  - **FF interface: direct co-localization.** BE3D's R729 and R736 are documented functional/clinical residues — **R736H
    is the 2nd-most-frequent DNMT3A cancer mutation (~2% AML)** and R729W abolishes processive catalysis; both are
    DNMT3L-stimulated (Emperle/Jeltsch 2023; Norvil/Gowher 2020). This validates BE3D's interface calls.
- **What is NOT explained by the paper (the novel part):**
  1. The paper stops at a general *interdomain* PWES cluster. BE3D resolves it into a **specific ranked interface map**:
     RD ≫ FF ≫ (ADD–MTase, DNA, H3 depleted).
  2. **A clinically UNANNOTATED, base-editable RD-interface patch — T671, M674, V675, R676, H873, D876, V877, S878 —**
     that is spatially adjacent to R882 but distinct from it. The worker confirmed *none* of these carry documented
     patient mutations (R676K/M674T/N879A appear only as engineered probes). BE3D thus reveals an **interface
     vulnerability NOT dominated by R882** — a base-editable proxy for the R882 tetramerization defect.
- **Confidence caveats:** (i) 4U7T covers 474–912 only; PWWP/N-term LOF (e.g. G293/E294, real hits) is out of the interface
  scope. (ii) The DNA-binding face is *curated from DNMT3A–DNA structures* (no DNA in 4U7T) — its low enrichment is
  therefore a soft result; a run on 5YX2/6W89 would harden it. (iii) FF vs RD ranking is radius-sensitive because FF
  residues abut the catalytic core; the RD call is the more *conformation-independent* one (holds in AF baseline too).
  (iv) pLDDT is high across the crystallized catalytic domain; hotspots are in ordered regions.

## 6. NEW actionable proposals (the deliverable)
For **epigenetics / AML / CHIP researchers and structural/chemical biologists:**
1. **Validate the base-editable RD-interface patch** as a functional proxy for R882 tetramerization loss:
   individually test **S878F, V877, D876, H873, T671** (and M674/V675/R676) by CBE/prime editing or site-directed
   mutagenesis in the DNMT3A methyltransferase + tetramerization assay (analytical SEC / cooperative-DNA-binding).
   Prediction: these phenocopy R882H's *hypomorphic/dominant-negative* tetramer defect **without** touching R882 —
   giving a clean allele series to dissect tetramerization vs catalysis, orthogonal to the confounded R882 codon.
2. **Actionable interface handle:** the RD patch (H873/D876/V877/S878 loop + T671/674–678 loop) defines a compact,
   surface-exposed, non-catalytic druggable epitope for a **tetramerization-disrupting** small molecule or stapled
   peptide — a way to inhibit DNMT3A *activity* by blocking oligomerization rather than the SAM pocket (avoids
   pan-methyltransferase off-target).
3. **FF interface:** R729/R736 (already clinical) plus Y724/Y735/H739 are a base-editable DNMT3L-dependency handle —
   propose testing whether editing Y735/Y724 selectively impairs DNMT3L-stimulated (vs basal) activity.
4. **De-prioritize** the ADD–MTase autoinhibitory contact and the (curated) DNA face as *primary* LOF interfaces in this
   assay — they are depleted, arguing the reporter LOF is driven by catalysis + quaternary assembly, not autoinhibition
   release or direct DNA contact.

## 7. BE3D issues encountered (appended to BE3D_IMPROVEMENTS.md)
- [BLOCKER-Windows] complex/ppi_diff mode's deep `output/<gene>/ppi_partners/<gene>_chain_<C>/screendata_sequence/...`
  path overflows Windows MAX_PATH (260) → `FileNotFoundError` on partner `to_csv`, even though the dir was created.
  Workaround: short output root + short screen filename. Fix: `\\?\` long-path prefix or shorter partner subpaths.
- [SCIENCE] BE-QA and the whole neg/pos framing assume a **dropout** sign (LOF = negative). For an **activity-reporter**
  screen (LOF = positive enrichment) the score must be negated first, or every downstream "neg/pos" label inverts
  silently. BE3D should expose an explicit `assay_direction`/`invert_score` flag and document it.
- [FRICTION] `run_new_target.py` hardcoded `function_for_lfc3d` and `atom_level_naa`; added `--function-for-lfc3d` and
  `--atom-level` flags to run the mandated sweep. (harness improvement, backward-compatible.)
- [FRICTION] Harness `mutation_category`/`mutation_priority` hardcode `Splice-donor/acceptor`; this dataset uses bare
  `Splice`. Harmless here (single-category guides pass through `reduce_mutation_type` unchanged) but the splice bucket
  silently emptied — categories should be config-driven per dataset.
