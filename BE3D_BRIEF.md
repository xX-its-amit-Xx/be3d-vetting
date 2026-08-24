# BE3D — Shared technical brief (for delegated research agents)

**BE3D** (`broadinstitute/BE3D`; pip package `beclust3d` v1.0.0; MIT) is a Python package for
**interpreting structure-function relationships in base-editor (BE) tiling-mutagenesis screens** by
mapping mutation readouts (log fold change, LFC) onto 3D protein structures.
Authors: Calvin XiaoYang Hu (Harvard), Yoochan Myung (Broad Institute). Pairs with the
**Genomics 2 Proteins Portal (G2P)** at g2p.broadinstitute.org for interactive visualization.
Repo cloned locally at (this session): scratchpad/BE3D — has README.md, README_FUNCTIONS.md, and
beclust3d/ source (~6100 LOC: qa/, lfc3d/, aggregate/, helpers/).

## Scientific context
- **Base editing (BE)**: CRISPR-derived tools (cytosine/adenine base editors, CBE/ABE) that install
  point mutations without double-strand breaks. **Tiling mutagenesis screens** systematically install
  many variants across a gene's coding region, then use a functional readout (e.g. proliferation/dropout
  LFC) to find which residues matter.
- BE3D's insight: single-residue screen signal is noisy; **aggregating signal over 3D spatial
  neighborhoods** (LFC3D) surfaces functionally important structural hotspots that per-residue analysis misses.

## The 3 modules & pipeline (function-level)
1. **BE-QA** — `hypothesis_test`: Mann-Whitney U + Kolmogorov-Smirnov on LFC distributions, knockout
   (nonsense/splice) vs neutral (silent/no-mutation) guides. H1 = within one screen; H2 = case in one
   screen vs pooled controls across screens. Screen quality gate.
2. **BE-Clust3D** — the core:
   - `sequence_structural_features`: UniProt seq + AlphaFold/PDB structure + DSSP secondary structure,
     builds per-residue neighbor lists within `radius` (default 6Å; residue- or atom-level).
   - `conservation`: MUSCLE/CLUSTAL align two seqs (cross-species/isoform) → per-residue conservation.
   - `parse_be_data` / `plot_rawdata` / `randomize_data`: raw guide→mutation tables, plots, null baseline.
   - `prioritize_by_sequence` / `randomize_sequence` / `plot_screendata_sequence`: aggregate multiple edit
     types per residue (mean/min/max/sum), with p-values.
   - `calculate_lfc3d`: **the LFC3D score** — aggregate per-residue LFC over the 3D spatial neighborhood
     (supports protein-protein interaction chains). Randomization gives a null.
   - `average_split_score` / `bin_score` / `znorm_score`: split +/- effects, percentile bins,
     z-normalize vs null, assign significance at p<0.05/0.01/0.001.
   - `clustering` / `plot_clustering`: agglomerative (single-linkage, euclidean) spatial clustering of
     significant residues over distance thresholds → hotspots; dendrograms.
   - Characterization: `enrichment_test` (Fisher's exact on structural features among hits),
     `plot_enrichment_test`, `lfc_lfc3d_scatter`, `pLDDT_RSA_scatter`, `hits_feature_barplot`.
3. **BE-MetaClust3D** — `average_split_meta`, `bin_meta`, `znorm_meta`: meta-aggregate LFC3D across
   many screens (SUM etc.) → consensus hotspots, rescue weak signals; supports cross-species via alignment.
- `g2p_formatted_hit_cluster`: export hits/clusters for the G2P portal.

## Inputs / outputs / run modes
- Inputs: BE screen scores TSV (cols: Mutation_type, Mutation_list, Gene, sgRNA_score) + UniProt ID
  (auto-fetch AlphaFold, isoform-aware) OR custom FASTA + PDB (+ optional DSSP, alignment).
- Outputs: per-residue LFC/LFC3D tables, significance labels, clusters, characterization plots,
  G2P-compatible TSV for interactive 3D viewing.
- Run: Google Colab (zero-install) or local conda env + YAML config (`be3d_local.py`). Example: KBTBD4 (Yeo et al.).

## Where it sits
It is NOT a base-editor outcome/efficiency predictor (that's BE-Hive/BE-FF/BEdeepon), NOT a guide designer,
NOT a variant pathogenicity predictor (AlphaMissense/ESM). It is a **downstream analysis / interpretation**
layer: turning a completed BE tiling screen into 3D structural hotspots with statistics — bridging
functional genomics screens and protein structural biology.
