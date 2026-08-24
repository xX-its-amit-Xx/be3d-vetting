# BE3D (`beclust3d` v1.0.0) — Design Deep-Dive (Queen 4)

Vetting target: `broadinstitute/BE3D`, pip package `beclust3d` v1.0.0, MIT. A Python
package for interpreting base-editor (BE) tiling-mutagenesis screens by mapping per-residue
log-fold-change (LFC) readouts onto 3D protein structure (the **LFC3D** score) and finding
statistically significant structural hotspots. ~6,100 LOC in `beclust3d/`.

All line references are to the locally cloned repo at
`.../scratchpad/BE3D/beclust3d/...` and `.../scratchpad/BE3D/examples/...`.

---

## 1. ARCHITECTURE MAP

### 1.1 Package layout

```
beclust3d/
  __init__.py                 # flat re-export of the whole public API (25 names)
  qa/
    hypothesis_tests.py        # hypothesis_test() — screen QA gate
    hypothesis_tests_helpers.py# hypothesis_one/two, MW/KS, plotting
  lfc3d/                       # the core BE-Clust3D module
    structure.py               # sequence_structural_features(), _lite variant
    structure_helpers.py       # 903 LOC: UniProt/AF/DSSP fetch+parse, neighbor lists, burial
    conservation.py            # conservation() — MUSCLE/Clustal cross-species/isoform align
    preprocess_data.py         # parse_be_data(), sanitary_check()
    preprocess_data_helpers.py # mutation parsing (identify_mutations, reduce_mutation_type, aa_map)
    preprocess_data_plot.py    # plot_rawdata()
    prioritize_sequence.py     # prioritize_by_sequence() — per-residue LFC + z/p
    prioritize_sequence_plot.py# plot_screendata_sequence()
    randomize_data.py          # randomize_data() — per-guide LFC permutation null
    randomize_sequence.py      # randomize_sequence() — per-residue null
    calculate_lfc3d.py         # calculate_lfc3d() — THE LFC3D SCORE (6A neighborhood aggr.)
    clustering.py              # clustering(), cluster_residues_from_pdb()
    clustering_plot.py         # plot_clustering()
    characterization.py        # enrichment_test() — Fisher's exact on structural features
    characterization_plot.py   # enrichment/scatter/barplot
  aggregate/
    nonaggregate.py            # average_split_score, bin_score, znorm_score (per-screen)
    metaaggregate.py           # average_split_meta, bin_meta, znorm_meta (across screens)
    aggregate_helpers.py       # binning, calculate_stats (z/p), pooled_mean_std, mu_sigma_screens
    aggregate_plot.py          # average_split_bin_plots()
  helpers/
    align/muscle-osx-arm64.v5.3   # vendored MUSCLE binary (macOS arm64 only!)
    visualization/g2p.py       # g2p_formatted_hit_cluster() — Genomics2Proteins portal export
examples/
  be3d_local.py (1608 LOC)     # YAML-driven end-to-end driver (main(), modes, PPI, blind, ppi_diff)
  be3d_local_helper.py, be3d_plotly.py
  BE3D_Colab.ipynb, BE3D_local.ipynb
  yaml/*.yaml                  # 8 example configs (KBTBD4/HDAC1 monomer, PPI, complex, cross-species MORC2)
  data/, pdb/                  # KBTBD4-HDAC1 ABE/CBE screens, MORC2 screens, example PDBs
```

There is **no `beclust3d` package-level CLI or `main()`**; the orchestration lives entirely in
`examples/be3d_local.py`, which imports the library functions and wires them together from a YAML config.

### 1.2 Data-flow pipeline (DataFrames + on-disk artifacts)

Working directory subdirs are created lazily by each stage (`os.mkdir` guarded by `os.path.exists`):

```
raw screen TSV(s)  ── parse_be_data ──►  screendata/{gene}_{screen}_{Missense,Silent,Nonsense,No_Mutation}.tsv
  (Mutation cat,                          (per-guide, exploded to one row per amino-acid edit:
   LFC, gene, edits)                       this_edit, edit_pos, refAA, altAA, LFC)
                                          + randomize_data ──► screendata_rand/..._Missense_rand.tsv.gz
                                                                (adds LFCr1..LFCrN permutation cols)

UniProt+AlphaFold+DSSP ── sequence_structural_features ──► sequence_structure/*_coord_struc_features.tsv
                                          (unipos, unires, chain, x/y/z_coord, bfactor_pLDDT,
                                           SS9/SS3, ACC, RSA, exposure, PHI/PSI, dBurial,
                                           normSumdBurial, pLDDT_dis, Naa_count, Naa, Naa_pos, Naa_chain)

df_dict{mut:df} + df_struc + df_consrv + df_control
   ── prioritize_by_sequence ──► screendata_sequence/{gene}_{screen}_protein_edits.tsv
        (one row per structure residue: {func}_{mut}_LFC, _stdev, _Z, _p, _plab, _plab_thr,
         all_{mut}_edits, conservation)
   ── randomize_sequence ──► screendata_sequence_rand/..._protein_edits_rand.tsv.gz
        (per-residue {func}_Missense_LFCr1..N, from averaging per-guide rand over edit_pos)

df_struc + [protein_edits...] + [protein_edits_rand...]
   ── calculate_lfc3d ──► LFC3D/{gene_type}_{gene}_LFC_LFC3D_LFC3Dr.tsv.gz
        (per residue per screen: {screen}_LFC, _LFC_Z, _LFC3D, LFCr1..N, LFC3Dr1..N,
         AVG_LFCr[_neg/_pos], AVG_LFC3Dr[_neg/_pos])

   ── average_split_score / bin_score / znorm_score  (per screen; nonaggregate.py)
        ──► LFC3D/{gene_type}_{gene}_{LFC3D}_bidirectional.tsv, _dis_wght.tsv, NonAggr_LFC3D.tsv
             (_neg/_pos split; percentile bins+weights; *_z, *_p, *_psig at p<.05/.01/.001)
   ── average_split_meta / bin_meta / znorm_meta  (across screens; metaaggregate.py)
        ──► meta-aggregate/{gene}_{LFC3D}_bidirectional.tsv, _dis_wght.tsv, MetaAggr_LFC3D.tsv
             (SUM_LFC3D[_neg/_pos], SUM_LFC3Dr..., SUM_LFC3D_{neg,pos}_{05,01,001}_psig)

df_struc + df_pvals(_psig cols) ── clustering ──► cluster_{LFC,LFC3D,union}/{gene}_{screen}_Aggr_Hits.tsv
        (adds {col}_Clust_{d}A cluster labels for d=4..max_distances)

everything ── g2p_formatted_hit_cluster ──► g2p_visualization/{,POS_,NEG_}lfc_lfc3d_union_cluster_for_g2p.tsv
characterization ── enrichment_test ──► characterization/{gene}_enrichment_test.pickle
```

### 1.3 Core data model (columns)

- **Identity/coord**: `unipos` (int, 1-based UniProt position; note: PDB residue numbering is
  assumed == unipos throughout), `unires` (1-letter AA), `chain`, `x_coord/y_coord/z_coord`
  (CA coordinates; `'-'` sentinel when missing), `bfactor_pLDDT` (B-factor = AlphaFold pLDDT).
- **Neighborhood**: `Naa_count`, `Naa` (`;`-joined neighbor 1-letter AAs), `Naa_pos`
  (`;`-joined neighbor positions), `Naa_chain` (`;`-joined neighbor chains). This is the
  precomputed 6A adjacency list that LFC3D smoothing consumes.
- **DSSP/exposure**: `SS9`/`SS3` (secondary structure), `ACC`, `RSA`, `exposure`
  (core/buried/medburied/medexposed/exposed), `PHI/normPHI/PSI/normPSI`, `dBurial`,
  `normSumdBurial`, `pLDDT_dis` (very low/low/confident/high).
- **Signal**: `{screen}_LFC`, `{screen}_LFC_Z`, `{screen}_LFC3D`; randomized `..._LFCr{r}`,
  `..._LFC3Dr{r}`; averages `..._AVG_LFCr[_neg/_pos]`, `..._AVG_LFC3Dr[_neg/_pos]`.
- **Significance**: `*_z`, `*_p`, `*_psig` (values `p<0.05`/`p>=0.05` etc.), `*_plab`
  (`-p=0.01`, `+p=0.001`, …), `*_dis` (percentile bin label), `*_wght` (bin weight −0.95..+0.95).
- **Conservation**: `conservation` (conserved/similar/weakly_similar/not_conserved),
  `v_score` (3/2/1/−1), `alternate_res`, `alternate_res_pos`.

### 1.4 Dependency / call graph (text)

```
examples/be3d_local.py : main()
  ├─ conservation()                       [optional cross-species/isoform]
  ├─ sequence_structural_features()
  │     ├─ query_uniprot / parse_uniprot
  │     ├─ query_af / parse_af / update_pdb_element_symbols / parse_coord
  │     ├─ query_domains|parse_domains
  │     ├─ run_dssp / parse_dssp
  │     ├─ count_aa_within_radius            (residue-level, CA-CA)      ── OR ──
  │     │   count_residue_contacts_all_atoms_single (atom-level, min atom-atom; Biopython NeighborSearch)
  │     └─ degree_of_burial
  ├─ hypothesis_test()  ── hypothesis_one/two ── ks_2samp / mannwhitneyu ── hypothesis_plot
  │        (qa_only → sys.exit; qa_passed_only → filter screens by H2 KS p<0.05)
  ├─ parse_be_data()  ── identify_mutations / reduce_mutation_type / aa_map
  ├─ plot_rawdata()
  ├─ randomize_data()                        [per-guide permutation]
  ├─ for each screen:
  │     prioritize_by_sequence()  ── get_plabel ; scipy.stats.norm ; statistics.NormalDist
  │     randomize_sequence()                 [per-residue null via groupby(edit_pos).mean()]
  │     plot_screendata_sequence()
  └─ run_Clust3D_per_species():
        clustering()  (LFC single-screen)  ── sklearn.AgglomerativeClustering / cluster_residues_from_pdb
        calculate_lfc3d()  ── _resolve_neighbor_sources / _gather_values
        average_split_score → bin_score → znorm_score            (nonaggregate)
        [if >1 screen] average_split_meta → bin_meta → znorm_meta (metaaggregate)
        clustering() (LFC & LFC3D, neg/pos × 3 pthr) → plot_clustering()
        enrichment_test() / plot_enrichment_test / scatter / barplot
        g2p_formatted_hit_cluster()
```

`be3d_local.py` also implements advanced run modes beyond a plain monomer: `complex`
(`run_complex_mode`, each chain in turn is target, others lite-preprocessed as PPI neighbors),
`ppi_diff` (`run_ppi_diff_mode`, runs PPI and no-PPI legs then paints ppi/noppi/delta onto a PDB
B-factor), and a `blind_target` mode (`run_blind_target` / `compute_blind_lfc3d`) that derives a
purely partner-contributed LFC3D for a scaffold chain with no screen of its own — explicitly noted
in-code as descriptive only, since a purely cross-chain value has zero variance in the null and
therefore no valid p-value (be3d_local.py:193-217).

---

## 2. ALGORITHM INTERNALS

### 2.1 The 3D neighborhood / LFC3D score (`calculate_lfc3d.py`)

**Neighbor topology is precomputed upstream** in `structure_helpers.py`, then LFC3D just aggregates
LFC values over that fixed adjacency list.

- **Residue-level neighbors** — `count_aa_within_radius` (structure_helpers.py:439-519). Brute-force
  O(N²) double loop over residues; for each target residue on `target_chain` with valid CA coords,
  computes Euclidean CA-CA distance (`math.sqrt(dx²+dy²+dz²)`, line 497) to every other residue and
  keeps those with `pairwise_dist <= radius` (default 6.0A, line 500). Writes `Naa/Naa_pos/Naa_chain`.
  Missing coords (`'-'`) are skipped and flag an "Incomplete Structure" warning.
- **Atom-level neighbors** — `count_residue_contacts_all_atoms_single` (structure_helpers.py:650-903).
  Two neighbors are within radius if their **minimum heavy-atom-to-atom distance <= radius**. Fast path
  uses Biopython `PDBParser` + `NeighborSearch` (KD-tree, line 783/790); a pure-Python brute-force
  fallback parser exists if Biopython import fails (line 820-875). Selected via
  `atom_level_naa=True` on `sequence_structural_features` (structure.py:129-132).

**LFC3D computation** — `calculate_lfc3d()` (calculate_lfc3d.py:21-239):
1. Builds `naa_pos_chain_dict` mapping each residue to its `chain_pos` neighbor strings (lines 111-122).
2. Per residue, `_resolve_neighbor_sources` (lines 241-286) builds a fixed list of value sources =
   `[('local', aa)]` (the residue itself) **plus** each neighbor as either `('local', naa_idx)` for
   same-chain neighbors or `('cross', gene_identifier, naa_idx)` for PPI cross-chain neighbors
   (lines 263-284). Cross-chain PPI neighbors are never conservation-gated (line 277); chains present
   in the PDB but not listed in `ppi_chain_gene_dict` are silently ignored (line 279-280).
   This is computed **once per screen** because topology does not depend on the randomization index.
3. `_gather_values` (lines 288-305) fetches the LFC value for each source (local from `taa_LFC_dict`,
   cross from `ppi_edits_dict`), drops `'-'`, and returns floats.
4. The per-residue LFC3D = `function_aggr_lfc3d(taa_naa_LFC_vals)` (line 172), where the aggregation
   function is tunable (`function_type_lfc3d`, default `'mean'`; `func_map` supports
   mean/median/sum/min/max, line 38). So **LFC3D of residue i = mean of {LFC(i)} ∪ {LFC(j) : j within
   6A}** — i.e. the residue's own LFC smoothed together with its spatial neighbors' LFCs.
5. The identical aggregation is repeated for all `nRandom` (default 1000) randomized LFC columns
   (lines 180-202) to build a per-residue null. Then `AVG_LFC3Dr`, `AVG_LFC3Dr_neg`, `AVG_LFC3Dr_pos`
   are the mean over the 1000 nulls (all / negatives-only / positives-only) (lines 220-229). Note the
   two-mean design flagged in-code (lines 18-20): the LFC3D aggregator is tunable, but the average over
   randomizations is hard-coded mean.

### 2.2 Randomization null model

Two levels:
- **`randomize_data`** (randomize_data.py:14-87) — the actual null. Takes the per-guide Missense LFC
  list and produces `nRandom` (default 1000) **permutations of the LFC values across guides**
  (`np.random.permutation(LFC_list)`, line 79). With `seed=True` uses
  `np.random.default_rng(i)` seeded by iteration index i (line 76-77) — reproducible. The driver
  `main()` calls it with `seed=True` (be3d_local.py:742). So the null breaks the guide→position
  association while preserving the LFC value distribution (a label-permutation / positional null).
- **`randomize_sequence`** (randomize_sequence.py:18-129) — collapses per-guide randomized columns to
  per-residue by `df_rand.groupby('edit_pos')[rand_columns].mean()` (line 102), i.e. the same mean-over-
  edits-at-a-position aggregation used for the real per-residue LFC, then maps back to every structure
  row. Non-matched / non-conserved positions get `'-'` (lines 110-114).

So the null is: permute guide LFCs → re-aggregate to residues → re-smooth over 3D neighbors →
compare observed LFC3D to this null distribution. `nRandom` default 1000 (per-screen), 500 in
`average_split_meta`.

### 2.3 z-normalization & p-value assignment

- **Per-residue sequence-level (LFC)** — `prioritize_by_sequence` (prioritize_sequence.py:193-224).
  Background μ/σ come from the **No_Mutation control guides**, computed *separately for the negative
  and positive halves* of the control LFC distribution (lines 126-136: `mu_neg,sigma_neg` from
  controls <0; `mu_pos,sigma_pos` from controls >0). For each residue, if LFC<0 z-scored vs the neg
  control model, if >0 vs the pos model (lines 204-213). p-value is **one-sided**
  `norm.sf(abs(z))` (line 206/211). `plab_thr` = `p<pthr` vs `p>=pthr`. `get_plabel` (lines 230-252)
  maps z to hard-coded one-tailed cutoffs: |z|>1.64→p=0.05, >2.32→p=0.01, >3.09→p=0.001 (a two-tailed
  variant is commented out).
- **LFC3D per-screen** — `znorm_score` (nonaggregate.py:180-304). The null μ/σ are the mean/std of the
  **randomized LFC3D averages** (`AVG_LFC3Dr_neg`/`_pos`) across residues (lines 261-270), computed
  separately for neg and pos. Then `calculate_stats` (aggregate_helpers.py:39-54) gives
  `z = NormalDist(mu,sigma).zscore(signal)`, one-sided `p = norm.sf(abs(z))`, and label
  `p<pthr` iff `p<pthr AND abs(signal)>abs(mu)` (line 53). Done for each pthr in `[0.05,0.01,0.001]`
  producing `_neg/_pos_{05,01,001}_psig`.
- **Meta-aggregated** — `znorm_meta` (metaaggregate.py:232-356) is identical but the null μ/σ come from
  the SUM-aggregated randomized meta columns (`SUM_LFC3Dr_neg/_pos`).

Note: LFC z-scoring uses control-guide μ/σ; LFC3D z-scoring uses the randomization-null μ/σ. Both are
**one-tailed by direction** (neg and pos scored against separate models), which is the intended design
but does mean significance is direction-conditional, not a symmetric two-tailed test.

### 2.4 Agglomerative clustering (`clustering.py`)

- `sklearn.cluster.AgglomerativeClustering` with `clustering_kwargs =
  {"n_clusters": None, "metric": "euclidean", "linkage": "single"}` (clustering.py:31) run on the
  **CA x/y/z coordinates** of only those residues whose `_psig` column equals the requested cutoff
  (e.g. `p<0.05`) (lines 123-147).
- It sweeps a range of **distance thresholds**: `distances = [4,5,...,max_distances]`
  (line 117; `int(i+1) for i in range(3,max_distances)` — starts at 4, not 3; driver passes
  `max_distances=20`). At each threshold it refits and records the number of clusters (elbow-style
  curve, `yvalue_lists`). "Hotspots" = the spatial clusters of significant residues at a chosen
  distance (default 6A used by the G2P export).
- Atom-level option `cluster_residues_from_pdb` (lines 164-218) builds a **precomputed min heavy-atom
  distance matrix** (cKDTree per residue) and clusters with `metric="precomputed"`, `linkage` default
  `average` when called from the driver (line 144) — note this diverges from the residue-level `single`
  linkage. There is also an implicit "union" clustering (cluster_union dir) combining LFC and LFC3D hits.

### 2.5 Meta-aggregation (`metaaggregate.py`)

`average_split_meta` (lines 20-147): for each residue, gathers each screen's LFC3D, splits into neg
and pos lists **independently**, and applies the meta aggregation function (default
`aggr_func_name='SUM'`, `func_map` supports MEAN/MIN/MAX/MEDIAN/SUM, line 28) separately to the neg
list, the pos list, and the combined list (lines 87-105). Because neg and pos are aggregated
independently, a residue can carry **both** a `SUM_LFC3D_neg` and a `SUM_LFC3D_pos` if screens disagree
in sign. Randomized meta nulls are built the same way per randomization (lines 114-133) then averaged
(lines 135-142). SUM is the default because summing per-residue LFC3D across screens reinforces
consistent signal and can rescue weak-but-concordant hits — the stated purpose of meta-aggregation.

### 2.6 QA hypothesis tests (`qa/`)

`hypothesis_test` (hypothesis_tests.py:16-145) runs both **Mann-Whitney U** (`mannwhitneyu(...,
method="asymptotic")`) and **Kolmogorov-Smirnov** (`ks_2samp`) two-sample tests comparing case vs
control LFC distributions (hypothesis_tests_helpers.py:134-147):
- **Hypothesis 1** (`hypothesis_one`): within a single screen, per gene — do case (e.g. Nonsense/Splice)
  LFCs differ from control (e.g. No Mutation/Silent) LFCs? (lines 22-67)
- **Hypothesis 2** (`hypothesis_two`): a gene's case LFCs in one screen vs that screen's control LFCs
  **pooled across all genes** in the screen (lines 72-132).
Empty groups return sentinel `[-999,-999]` (line 147). The driver uses H2 KS p<0.05 as the screen QA
gate (`qa_passed_only`, be3d_local.py:660-694): screens failing are dropped, and a `QA_STATUS.txt`
marker is written; a fully failed gene exits non-zero.

### 2.7 Characterization (`characterization.py`)

`enrichment_test` (lines 21-114) runs **Fisher's exact test** (scipy `fisher_exact`, two-sided) on a
2×2 of {below/above hit threshold} × {in/out feature set} (lines 133-140), reports **log2 odds ratio**
and its log2 CI via `scipy.stats.contingency.odds_ratio` (lines 141-147). Used to ask whether hits are
enriched in domains, low-pLDDT, secondary-structure categories, etc. Results are pickled (line 111-113).

---

## 3. ENGINEERING QUALITY ASSESSMENT

### Strengths
- **Clear modular decomposition** by pipeline stage (qa/lfc3d/aggregate/helpers) and a clean flat public
  API in `__init__.py`. Every public function has a full NumPy-style docstring with typed params.
- **Reproducibility of the null is available**: `randomize_data(seed=True)` uses per-iteration seeded
  `np.random.default_rng(i)` (randomize_data.py:76), and the driver enables it (be3d_local.py:742).
- **Defensive input contracts**: nearly every function opens with `assert` checks on required columns
  and length consistency (e.g. calculate_lfc3d.py:96-106, clustering.py:99-111, prioritize_sequence.py:99-111).
- **Real performance thought in the hot path**: `calculate_lfc3d` precomputes neighbor sources once per
  screen and reuses them across all 1000 randomizations (calculate_lfc3d.py:142-155, docstring at 249-258);
  `randomize_sequence` replaces a per-row re-filter with a single `groupby().mean()` (randomize_sequence.py:100-107);
  atom-level neighbors use a KD-tree. These are documented, deliberate optimizations.
- **Graceful degradation**: pure-Python PDB parser fallback if Biopython is unavailable
  (structure_helpers.py:820); remote MUSCLE API fallback if local MUSCLE/Clustal cannot run
  (conservation.py:99-118); user-supplied FASTA/PDB/DSSP/alignment all bypass the network fetches.
- **Intermediate artifacts persisted at every stage** (TSV / gzipped TSV), so runs are resumable and
  auditable, and `ppi_diff` uses `RUN_COMPLETED.txt` markers to skip completed legs.
- **Portability is documented**: README explicitly warns arm/M-series Macs cannot run Clustal locally
  and to use `mode='query'`; DSSP/MUSCLE install paths are given.

### Weaknesses / concrete smells
- **No tests and no CI.** The only files matching `*test*` are the QA *hypothesis*-test modules; there
  is no `tests/` directory, no `.github/` workflows, no pytest config. For a statistical method this is
  the single biggest gap — there is nothing verifying the LFC3D math, the null, or the z/p assignment.
- **`warnings.filterwarnings('ignore')` at import in many modules** (calculate_lfc3d.py:16,
  randomize_sequence.py:16, nonaggregate.py:16, metaaggregate.py:16, prioritize_sequence.py:18). This
  globally silences all warnings for any process that imports the package — including the pandas
  chained-assignment and dtype warnings the code itself triggers, hiding real correctness signals.
- **`'-'` string sentinel for missing numeric data everywhere**, mixed into otherwise-float columns,
  forcing repeated `replace('-', np.nan)` / `astype(float)` round-trips and `object` dtypes
  (e.g. nonaggregate.py:155-166, calculate_lfc3d.py:206-207). Fragile and error-prone; a proper NaN would
  be safer. Type coercion via `pd.to_numeric(errors='coerce')` can silently turn bad values into NaN.
- **Chained-assignment / SettingWithCopy patterns**: `prioritize_by_sequence` builds `df_protein` from a
  slice `df_struc[required_columns]` then assigns new columns to it (prioritize_sequence.py:114-122) —
  exactly what the silenced warnings would flag.
- **Hard PDB-numbering assumption**: throughout, `unipos` (UniProt position) is assumed equal to PDB
  `residue_number` (parse_coord matches on `residue_number == unipos`, structure_helpers.py:187;
  neighbor `naa_idx = int(naa_pos)-1`, calculate_lfc3d.py:270; B-factor painting keyed by residue_number).
  There is an AA-identity check with a warning on mismatch (structure_helpers.py:202-212) but no
  alignment/renumbering — user PDBs whose residue numbers don't match UniProt will silently drop residues.
- **Statistical model caveats worth noting for the science review**: LFC z-scores assume the control
  (No_Mutation) guide LFCs are Normal, split into independent neg/pos half-Normals; LFC3D z-scores treat
  the randomization-null mean/std as Normal. Significance is one-tailed per direction. The `bin_score`
  percentile weights (aggregate_helpers.py:78-90, fixed −0.95..+0.95) are heuristic. None of this is wrong
  per se, but it is undocumented modeling choice, untested, and the two-tailed alternatives are only
  present as commented-out code (prioritize_sequence.py:241-247, nonaggregate.py:253-259).
- **External-API fragility**: hard dependence on live `rest.uniprot.org` (FASTA + domains JSON),
  `alphafold.ebi.ac.uk` (PDB by `AF-{uniprot}-F1-model_v6` name — version pinned into the filename,
  structure.py/be3d_local.py), and EBI MUSCLE REST. `query_uniprot` uses the `wget` PyPI package
  (unmaintained) rather than `requests`. `query_domains` prints and returns silently on non-200
  (structure_helpers.py:384-386) rather than raising. No retry/backoff, no caching beyond "file exists".
- **Packaging inconsistency**: `pyproject.toml` uses **hatchling** as the build backend (lines 1-3) but
  then configures `[tool.setuptools.packages.find]` (lines 29-30) — setuptools config under a hatchling
  build is ignored; the two are mismatched. Dependencies are unpinned in `pyproject` but exhaustively
  pinned in `environment.yml` (a full conda freeze, Linux-specific), with a separate `environment_arm.yml`.
  README also references `pip install git+https://github.com/broadinstitute/beclust3d-public.git` — a
  *different* repo name than the vetted `broadinstitute/BE3D`.
- **A vendored macOS-arm64-only MUSCLE binary** ships in the wheel path
  (`helpers/align/muscle-osx-arm64.v5.3`) — dead weight / non-portable on Linux/Windows and not wired to
  `muscle_path` (which defaults to a bare `'muscle'` on PATH).
- **The real orchestration is a 1608-line untested script** (`examples/be3d_local.py`) with a giant
  `main(**kwargs)` that unpacks ~35 kwargs positionally (lines 530-571), `sys.exit()` calls inside library
  flow (lines 658, 694), a `# TODO: this should be moved to helper` (line 582), and tab-indented source.
  Much real logic (PPI, blind, ppi_diff modes) lives here rather than in the installable package.
- **O(N²) residue neighbor search** (`count_aa_within_radius`) is pure-Python double loop — fine for a
  single protein, slow for large complexes; the atom-level path is KD-tree accelerated but the default
  residue path is not.

### Config approach
YAML-driven via `examples/yaml/*.yaml` loaded by `load_config` → `yaml.safe_load` (be3d_local.py:10-14);
8 example configs cover monomer, chain-B, PPI-diff, complex, and cross-species MORC2 (each with a
`_colab` variant). The config is copied into the output dir for provenance (be3d_local.py:591).

### Reproducibility summary
Randomization is seedable and seeded by default; UniProt/AlphaFold versions are pinned into filenames;
the config is archived per run; intermediate artifacts are all written. But: no environment lock beyond
the conda freeze, no test suite to catch regressions, external DBs can drift, and global warning
suppression hides drift.

---

## 4. PUBLICATION & LAB

_(Filled in from the parallel literature-resolution agent — see summary below. Author line embedded in
every source file is: **Calvin XiaoYang Hu, Yoochan Myung, Surya Kiran Mani, Sumaiya Iqbal**;
`pyproject.toml` lists Calvin XiaoYang Hu (xiaohu@g.harvard.edu, Harvard) and Yoochan Myung
(ymyung@broadinstitute.org, Broad Institute). Pairs with the Genomics 2 Proteins Portal at
g2p.broadinstitute.org. Example data: KBTBD4/HDAC1 ABE+CBE screens (Yeo et al.), and MORC2 cross-species
screens.)_

### 4.1 Authorship (from the code itself)
Every source file header credits **Calvin XiaoYang Hu, Yoochan Myung, Surya Kiran Mani, Sumaiya Iqbal**.
`pyproject.toml` lists two authors: Calvin XiaoYang Hu (xiaohu@g.harvard.edu, Harvard) and Yoochan Myung
(ymyung@broadinstitute.org, Broad Institute). Output is "compatible with the Genomics 2 Proteins Portal"
(g2p.broadinstitute.org).

### 4.2 The lab / PI — "APAR" is INCORRECT; the lab is the Iqbal Lab
- The correct lab is the **Iqbal Lab**, PI **Sumaiya Iqbal**, in the **Ladders 2 Cures (L2C) Accelerator
  at the Broad Institute of MIT and Harvard** (Iqbal is also affiliated with ATGU at Mass General /
  Harvard Medical School). Sources: https://www.iqballab.org/ , https://www.broadinstitute.org/bios/sumaiya-iqbal
- **"APAR lab" is NOT confirmed anywhere.** I grepped the entire repo: the only 3 matches of "APAR" are
  coincidental substrings inside base64-encoded image blobs in `examples/BE3D_local.ipynb` — not a real
  term. It is not a lab name, module name, or acronym in the code or on any lab page. **The user's "APAR
  lab" label should be corrected to "Iqbal Lab (Broad)."**
- BE3D is effectively an **Iqbal Lab (Broad, computational) × Liau Lab (Brian B. Liau, Harvard Chemistry &
  Chemical Biology, experimental base-editing screens)** collaboration. Calvin Hu is affiliated with the
  Liau Lab; Yoochan Myung is the Iqbal Lab postdoc who presents BE3D. Source: http://www.liaulab.com/publications

### 4.3 The paper — NO standalone peer-reviewed paper / DOI yet (as of Aug 2026)
- There is **no dedicated BE3D/beclust3d journal article or bioRxiv preprint** that could be confirmed. The
  README has no "how to cite" and no DOI. Treat any claim of a formal "BE3D paper" cautiously.
- BE3D so far exists as: (a) the GitHub tool (repo `broadinstitute/BE3D`; pip source is actually named
  **`beclust3d-public`**, `github.com/broadinstitute/beclust3d-public`), and (b) a conference talk —
  "BE3D: A Computational Workflow for Integrative Structure-Function Analysis of Base-Editor Tiling
  Mutagenesis Data," presented by Yoochan Myung at **ISMB/ECCB 2025** (Liverpool, Jul 2025). Source:
  https://www.iqballab.org/news
- The likely **application/companion paper** (same author set, base-editing tiling + 3D structure-function)
  is a PRC2 sequence-function atlas paper from the Liau + Iqbal labs, listed as submitted/2026 with no DOI:
  "Base editing charts a sequence-function atlas of Polycomb Repressive Complex 2 (PRC2) at amino acid
  resolution." Source: http://www.liaulab.com/publications
- ⚠️ **Do NOT cite** the ProTiler-Mut / He et al. Cell Systems paper (DOI 10.1016/j.cels.2026.101651,
  bioRxiv 10.1101/2025.04.17.649336) as BE3D — it is a different tool that surfaces in searches.

### 4.4 Companion: Genomics 2 Proteins (G2P) Portal paper — CONFIRMED (separate)
- "Genomics 2 Proteins portal: a resource and discovery tool for linking genetic screening outputs to
  protein sequences and structures," Kwon, Safer, Nguyen, Hoksza, May, Arbesfeld, Rubin, Campbell, Burgin,
  Iqbal. **Nature Methods, 2024. DOI: 10.1038/s41592-024-02409-0** (PMID 39294369; preprint bioRxiv
  10.1101/2024.01.02.573913). Portal: g2p.broadinstitute.org.

### 4.5 Canonical example datasets
- **KBTBD4 (confirmed the documented worked example).** Ships as `examples/data/YeoKBTBD4-HDAC12025-{ABE,CBE}-Screen[-PPI].tsv`
  and `examples/pdb/8voj_KBTBD4_HDAC1_chain_ABC.pdb`, with YAMLs for monomer, chain-B, PPI-diff, and complex.
  Ties to **Yeo et al.**, "Asymmetric Engagement of Dimeric CRL3^KBTBD4 by the Molecular Glue UM171 Licenses
  Degradation of HDAC1/2 Complexes," bioRxiv 2024, **DOI 10.1101/2024.05.14.593897** (PMID 38798619; Liau
  lab; companion code `liaulab/HDAC1_KBTBD4_base_editing_scanning_2024`). Part of the KBTBD4/UM171
  molecular-glue story published in Nature 2025 (10.1038/s41586-024-08532-4; 10.1038/s41586-024-08533-3).
- **MORC2 — a genuine second example dataset.** Ships as `examples/data/MORC2/*.tsv` (human + mouse
  in-vitro/in-vivo Apobec/TadA screens) with `examples/pdb/MORC2.pdb` and cross-species YAMLs
  (`cross_species_MORC2.yaml`) — this is the worked example for the conservation / cross-species path.
- **DNMT3A and MEN1 are NOT example datasets** — they appear only as generic docstring placeholders
  ("e.g., 'DNMT3A', 'MEN1'") throughout the source; there is no DNMT3A/MEN1 data or config in the repo.

### 4.6 Sources
- https://github.com/broadinstitute/BE3D ; https://github.com/broadinstitute/beclust3d-public
- https://www.iqballab.org/ ; https://www.iqballab.org/news ; https://www.iqballab.org/publications
- https://www.broadinstitute.org/bios/sumaiya-iqbal ; https://www.atgu.mgh.harvard.edu/people/sumaiya-iqbal/
- http://www.liaulab.com/publications
- https://www.nature.com/articles/s41592-024-02409-0 (G2P, Nat Methods 2024)
- https://www.biorxiv.org/content/10.1101/2024.05.14.593897v1 (Yeo et al., KBTBD4)
- https://www.nature.com/articles/s41586-024-08532-4 ; https://www.nature.com/articles/s41586-024-08533-3 (KBTBD4/UM171 Nature 2025)
