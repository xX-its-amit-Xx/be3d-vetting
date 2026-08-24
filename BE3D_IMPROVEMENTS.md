# BE3D — Improvement Log (from real usage)

Running log of concrete, actionable improvements to BE3D discovered while *actually using it*.
Append as you go. Each entry: **what broke / was awkward → why it matters → concrete fix.**
Tag severity: [BLOCKER] stops a run · [FRICTION] costs time/errors · [SCIENCE] affects correctness/interpretation · [ENH] nice-to-have.

## Portability / packaging (found running on Windows + Python 3.14 + pandas 3.x)
- [BLOCKER] `sequence_structural_features` unconditionally shells out to `mkdssp`; if the binary
  isn't installed the whole run dies at the structure stage — even though DSSP is only needed for the
  final *characterization* step. **Fix:** wrap `run_dssp` in try/except and continue with SS/RSA marked
  unavailable; or ship a pure-Python DSSP fallback (biotite/pydssp) so a missing binary can never abort
  LFC3D/clustering. Today the only workaround is to hand-generate a `user_dssp` file.
- [BLOCKER] `to_csv(filename, "\t", ...)` uses a positional `sep` — removed in pandas 3.x → TypeError.
  **Fix:** use keyword `sep="\t"` (nonaggregate.py, metaaggregate.py).
- [BLOCKER] Significance labels like `p<0.05` are baked into **output filenames**; `<`/`>=` are illegal
  in Windows filenames → run dies late (after long compute). **Fix:** sanitize path components (`<`→`lt`,
  `>=`→`ge`) while keeping the in-data labels intact (clustering_plot.py).
- [FRICTION] Positional `Series[0]` indexing in the runner's `find_union` breaks under modern pandas
  (label-based) → `KeyError: 0`. **Fix:** `list(series)[0]` (be3d_local.py).
- [FRICTION] `pip install` target in README is `beclust3d-public` but the repo is `BE3D`; pyproject uses
  hatchling yet configures setuptools packages.find — mismatch. **Fix:** align repo/package names & backend.
- [ENH] A macOS-arm64-only MUSCLE binary is vendored in the wheel path — dead weight/non-portable
  elsewhere and not wired to `muscle_path`. **Fix:** drop the vendored binary; document `muscle_path`.

## To be extended by workers during real-problem runs
(new entries below)

### DNMT3A run (Lue et al. 2023 activity-reporter screen; PDB 4U7T) — Queen "DNMT3A"
- [SCIENCE] **Dropout assumption baked into BE-QA and the neg/pos channels.** BE3D treats LOF as a *negative* score
  (dropout convention). The Lue DNMT3A screen is an **activity reporter** (CpG-methylation→citrine FACS) where LOF =
  citrine de-repression = **positive** enrichment (verified: Nonsense +0.73, Splice +1.05, Silent/Non-exon ≈ 0). If fed
  raw, every downstream "neg" hotspot is actually GOF/neutral and the biology inverts *silently* — QA still "passes"
  because MWU/KS are direction-agnostic, hiding the error. **Fix:** add an explicit `assay_direction` / `invert_score`
  config flag (dropout vs enrichment) and record it in outputs. Workaround used: negate the score during reformatting so
  LOF lands in the `*_LFC3D_neg` channel; emit a raw TSV for audit.
- [BLOCKER][Windows] **complex / ppi_diff mode overflows Windows MAX_PATH (260).** The partner tree
  `output/<gene>/ppi_partners/<gene>_chain_<C>/screendata_sequence/<gene>_<screen>_protein_edits.tsv` is deeply nested;
  with a normal-length output root + screen filename the file path exceeds 260 chars and `prioritize_by_sequence`'s
  `to_csv` dies with `FileNotFoundError` **after** the directory was created (misleading). Monomer mode is unaffected
  (shallower tree). **Fix:** open with the `\\?\` long-path prefix (or `os.makedirs`+`open` via extended path), or
  shorten the partner subpath. Workaround: run complex mode under a short root (`C:\Temp\d3`) with a short screen name.
- [FRICTION] `run_new_target.py` hardcoded `function_for_lfc3d='mean'` and `atom_level_naa=False`, so the mandated sweep
  over aggregation (mean/sum) and residue-vs-atom neighbors was impossible without editing the harness. Added
  `--function-for-lfc3d` and `--atom-level` flags (backward-compatible defaults).
- [FRICTION] `mutation_category`/`mutation_priority` are hardcoded to `Splice-donor`/`Splice-acceptor`; datasets that use
  a bare `Splice` category get an empty splice bucket with no warning. `reduce_mutation_type` leaves single-category
  guides intact so it was harmless here, but category vocab should be dataset-configurable and warn on unmatched values.
- [ENH] complex mode with a homo-oligomer of the SAME gene (chains A & C both DNMT3A, sharing one screen) works and is a
  clean way to add cross-protomer neighbor signal at a homodimer interface — but it required a hand-written YAML
  (`ppi_chain_gene_dict:{A:DNMT3A,C:DNMT3A}` + `partner_uniprot:{DNMT3A:...}`); worth a documented recipe/example, since
  the shipped examples only cover hetero-complexes (KBTBD4–HDAC1).

## Plug-and-play harness for arbitrary targets (P2; verified on LYZ / UniProt P61626)
- [BLOCKER] `be3d_local.py:main()` does `shutil.copy2(config_yaml, output_dir/basename(config_yaml))`
  to archive the config. If the config YAML you pass in already lives *inside* `output_dir`
  (an obvious place to put it), copy2 hits src==dst → `PermissionError: [WinError 32] file in use`
  on Windows (and `shutil.SameFileError` on POSIX), killing the run immediately at the structure
  stage. **Fix:** either guard the copy (`if os.path.abspath(src) != os.path.abspath(dst)`), or
  always keep the source YAML outside `output_dir`. The harness (`run_new_target.py`) writes its
  generated YAML to a sibling `be3d_configs/` dir for exactly this reason.
- [FRICTION] AlphaFold DB only hosts the canonical F1 fragment keyed by the BASE accession, so a
  target given as an *isoform* id (e.g. `Q9NVX7-2`) has no AF file to download — the URL 404s.
  `prep_target.py` must be handed the base accession; isoform-specific numbering needs a user PDB.
  Worth documenting prominently since the flagship KBTBD4 example uses an isoform uniprot with a
  hand-supplied `user_pdb`, which can mislead new users into passing an isoform id to the fetcher.
- [FRICTION] `sequence_structural_features` generates DSSP from the *parse_af-processed* PDB, but a
  hand-generated `user_dssp` (from gen_dssp.py) is built from the raw AF monomer. For a single-chain
  AF F1 model this is fine — raw and processed residue numbers both equal unipos, so parse_dssp's
  inscode match succeeds. It would break for multi-chain/renumbered inputs; DSSP still only affects
  the final characterization step, so it never blocks LFC3D/clustering regardless.
- [ENH] A single, self-contained monomer YAML (conservation off, PPI dict = null, user_pdb +
  user_dssp wired in) is all that's needed to run any new gene end-to-end. The `mutation_priority`
  list and `qa.cases/controls` are the only fields that realistically need per-screen tuning; the
  harness exposes those as flags and hardcodes the rest to sensible defaults.
- [SCIENCE] Smoke test used a *fabricated* screen (refAA taken from the real P61626 sequence so
  Missense/Silent/Nonsense edits validate; 40/20/12 coding guides + 25 `No Mutation` controls).
  It proves plumbing only — the LFC/LFC3D/cluster values are meaningless. A real run needs a real
  reformatted screen with the four columns (gene_col, mut_col, mut_list_col, val_col).

## PPI / complex mode & control-category handling (found running MCL1 / Hanna 2021, P2)
- [BLOCKER] The neutral-control category MUST be literally named "No Mutation". `parse_be_data`
  writes one file per raw category token (`screendata/{gene}_{screen}_{token}.tsv`), but both
  `main()` and `preprocess_ppi_partner` hardcode reading `{gene}_{screen}_No_Mutation.tsv`. A screen
  whose controls are tokened `UTR`/`Intron` (as the Hanna MCL1 library is) dies with
  `FileNotFoundError` at `prioritize_by_sequence`. **Fix:** either let the config name the control
  category, or map the `no_mutation` bucket → the control filename. Workaround used: remap control
  tokens to "No Mutation" in the reformatted TSV.
- [BLOCKER, PPI] `mode: complex` cannot accept a **structural-only partner chain** (a bound
  peptide / ligand that is not itself base-editing-tiled). `run_complex_mode` → `preprocess_ppi_partner`
  filters the shared screen to the partner gene and unconditionally reads that partner's
  `No_Mutation.tsv`; with zero partner rows it crashes (same FileNotFoundError, now under
  `ppi_partners/<PARTNER>_chain_<c>/`). This is the *common* PPI case — the interacting BH3 peptide
  (BAX in 3PK1, the SAHB in 3MK8) or a small molecule is never tiled — so the biologically relevant
  interface run is impossible. `blind_target` mode covers a data-less *target* chain but there is no
  symmetric option for a data-less *partner*. **Fix:** allow a partner chain with no screen data to be
  registered as a pure structural neighbor (its cross-chain LFC contribution is simply empty), or
  auto-synthesize an all-`-` partner edits table. Verified: complex mode DOES run fine when every
  partner chain has real screen data (MCL1 A/C homodimer in 3PK1 → exit 0).
- [SCIENCE] Complex/PPI mode with a data-less partner is numerically identical to monomer:
  cross-chain neighbor values are `-` and `_gather_values` skips them, so the target's LFC3D is
  unchanged. A complex-mode run on 3PK1 A/C reproduced the monomer run's top hotspots exactly. Users
  should not expect a bound-complex PDB alone to add interface signal — the partner needs its own
  screen. Document this, because it is easy to assume "complex mode on a co-crystal = interface-aware".
- [FRICTION] In monomer mode a multi-chain (complex) PDB is accepted and cross-chain neighbors are
  silently dropped (safe, no crash), but then BE3D emits **no interface annotation at all** — you must
  compute which target residues contact the partner with an external distance script. An optional
  "annotate cross-chain contacts" flag in monomer mode would make interface-vs-core questions
  first-class without needing partner screen data.
- [ENH] `znorm_score` always computes p<0.05/0.01/0.001 regardless of `pthr.single_screen`, so a
  p-threshold sweep is free (no re-run needed) — but only 0.05 propagates to clustering /
  characterization. Exposing the clustering threshold per-p-level would let one run yield the full
  significance-vs-cluster sweep.

## KRAS run (Coelho/Dincer 2024 drug-resistance BE screen; PDB 6OIM/6UT0 + AF) — Queen "KRAS/EGFR"
- [SCIENCE][data] **The turnkey recipe's column map for the Coelho ST2 sheet is WRONG** (schema drift).
  The brief says `Mutation_list ← Protein_Change` (3-letter HGVS `p.Gly12Cys`), `Mutation_type ← variant_classification`
  (Missense/Silent/…). In the actual file: `variant_classification` is `SNV`/`SNV;substitution`/`NA` (useless);
  the real category lives in **`most_severe_consequence`** (`missense`/`synonymous`/`stop codon`/`splice variant`/
  `intron`/None); `Protein_Change` is VEP `ref/alt` (`E/G;T/A;TE/AG`) with redundant multi-codon haplotype tokens,
  NOT HGVS. `swissprot` is versioned (`P01116.251`) so `=="P01116"` returns 0 rows. Edits must be reconstructed from
  `zip(Amino_Acid_Position, Edited_AA, New_AA)`. **Fix/lesson:** a "turnkey" tag needs the column map verified against
  a live header dump before a worker trusts it; the recipe should carry the exact observed column names + a 3-row sample.
- [SCIENCE][enrichment] Reconfirms the dropout-convention issue from the DNMT3A entry, from the other side:
  this is a **resistance ENRICHMENT** screen, so BE-QA (cases=Nonsense/Splice vs controls) is null by design
  (KS p=0.29/0.50) and the meaningful signal is the **positive** `*_LFC3D_pos` channel. QA "passing/failing" is
  not the right accept gate for enrichment assays — an `assay_direction` flag (as already proposed) would also let QA
  test the correct tail.
- [FRICTION][structure] **Drug-bound crystal constructs silently drop edits.** 6OIM/6UT0 carry engineered
  G12C + C51S/C80L/C118S and construct-specific residues ≥151 with non-canonical numbering in the HVR. BE3D maps edits
  by `unipos==residue_number` and sanity-matches refAA; mismatched positions are dropped from the mapped set with only a
  buried print, so a user can lose pocket residues (here pos 12/51/80/118) without noticing. **Fix:** emit an explicit
  "N edits unmapped due to refAA mismatch: [list]" summary to `RUN_COMPLETED.txt` / a QA file. Workaround: use the WT
  AlphaFold model as the primary structure and treat the liganded PDB as a corroborating overlay.
- [FRICTION][concurrency] `run_new_target.py` writes its config to a sibling `be3d_configs/<gene>_run.yaml` whose name
  depends only on the gene. Running the SAME gene many times in parallel (a parameter sweep) races on that one filename
  and on the config-copy into each output_dir. **Workaround that works:** give every run a unique *parent* dir
  (`runs/<run_id>/<gene>`) so each gets its own `be3d_configs/`. A `--config-name`/unique-suffix option would be cleaner.
