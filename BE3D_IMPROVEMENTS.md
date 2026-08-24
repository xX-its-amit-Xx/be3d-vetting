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

## Round-2 EGFR/BRAF resistance benchmark (Queen "EGFR/BRAF")

- [HIGH-VALUE][interpretation] **No base-editing REACHABILITY report.** BE3D silently scores only the
  amino-acid substitutions a CBE/ABE can install. On EGFR ~18/25 and on BRAF 12/15 curated clinical
  resistance/driver residues receive NO signal because their alleles need transversions/indels (T790M,
  C797S, L858R, exon19del; BRAF V600E, R509 — the latter two had ZERO guide coverage). A per-curated-
  residue "reachable by CBE/ABE? / #guides covering" annotation (or a warning when a requested hotspot is
  unreachable) would stop users misreading intrinsic BE limits as BE3D false negatives. This is the single
  biggest interpretability gap for oncogene-resistance use.
- [HIGH-VALUE][stats] **High positive base rate, no warning.** In both genes ~30–40% of scored residues are
  called significant in the enrichment direction. The A0 discrimination test only clears chance in 1 of 3
  drug arms (EGFR gefitinib). BE3D should print the base rate (#sig/#scored) and, ideally, an enrichment-
  vs-null-set statistic, so "recovered residue X" isn't over-trusted.
- [MEDIUM][structure] **Disordered/low-pLDDT residues dominate top positive hits.** EGFR's strongest
  positive LFC3D signals are the disordered C-terminal tail (res >1000, pLDDT<50); precision@10 for the
  functional set was ~0 because of it. A pLDDT gate or a low-confidence flag on hotspot calls would help.
- [MEDIUM][artifact] **Contiguous single-neighborhood runs unflagged.** Robust hits often come as adjacent
  runs (EGFR 806–812, BRAF 521–550) that may reflect a few high-LFC overlapping guides / bystander edits
  rather than independent per-residue signal. Flagging "contiguous stretch, possible bystander" would aid
  triage.
- [CONFIRMED-OK] The KRAS-derived Coelho ST2 reformatting recipe (category from most_severe_consequence;
  edits from zip(Amino_Acid_Position,Edited_AA,New_AA); ABE+CBE pooled; positive=resistance) ported cleanly
  to EGFR and BRAF (BRAF refAA-vs-UniProt mismatch only 1.9%). z-score columns are NOT gene-specific — every
  gene carries all cell-line×drug arms; pick the arm by the clinically matched cell line (PC9→EGFR, HT29→BRAF).
- [CONFIRMED-OK] Unlike KRAS (resistance-only, null QA), EGFR/BRAF KS-H1 is strongly significant
  (D=0.37 p=4e-52; D=0.26 p=6e-14) because the target is essential in the matched line — QA gate passes.

## Run: MAP2K1 / MAP2K2 (MEK1/MEK2, trametinib resistance) — Queen "MEK1/MEK2"
- [POSITIVE] The Coelho ST2 schema + `convert_kras.py` template transferred cleanly to two new genes
  (MEK1/MEK2): only 0.9% of guides dropped on the refAA-vs-structure sanity check. Reconstructing the
  reference sequence FROM the AlphaFold PDB (not a separately-downloaded FASTA) guarantees screen numbering
  == structure numbering and removes a whole class of off-by-one/isoform bugs — recommend baking this into
  `prep_target.py` (emit a `<uniprot>.seq` from the PDB it just downloaded).
- [FRICTION] Output-tree shape differs between harness invocations: the KRAS runs nested outputs under
  `runs/<id>/<GENE>/...` but these MEK runs put them directly under `runs/<id>/...` (LFC3D/, etc.). Any
  downstream parser must glob for `**/LFC3D/*NonAggr_LFC3D.tsv` rather than assume a `<GENE>/` level.
- [FRICTION] Column-prefix is `<GENE>_<screen-stem>` (e.g. `MAP2K1_Tram`), derived from the staged TSV
  filename, NOT documented. Robust parse = strip the `_LFC3D_pos` suffix off the first matching column
  (as `robust.py` does). A machine-readable `columns.json` in the output would remove the guesswork.
- [DISCRIMINATION / high base rate] On these resistance screens BE3D calls a LARGE fraction significant in
  the positive direction: MEK1 13–20%, **MEK2 24–32%** of all residues (p<0.001→p<0.05). Per the A0 spec
  this makes p<0.05 overlaps weak evidence. Two concrete asks: (a) report the base rate in the run summary
  so users see it, and (b) the null appears not to control the positive-tail FDR well when a broad swath of
  the protein is mildly enriched (paralog MEK2 far worse than MEK1 on the SAME drug/pipeline) — worth a
  look at whether the per-residue randomization is too permissive for enrichment (vs dropout) screens.
- [FALSE-POSITIVE geometry] The dominant FP mode was low-pLDDT regions: MEK1 284–287 (αF–αG loop, pLDDT
  40–46) and MEK2 N-terminal 1–32 (disordered) were robustly "significant" across the sweep. BE3D already
  writes pLDDT in characterization but does NOT down-weight or flag hotspots in low-pLDDT geometry.
  Recommend an automatic `low_pLDDT` flag (e.g. mean neighborhood pLDDT<60) on each called hotspot.
- [QA direction] Unlike the KRAS resistance arms (null QA), the MEK trametinib arms gave SIGNIFICANT KS
  (MEK1 D=0.21 p=1.9e-5; MEK2 D=0.19 p=6.2e-6) because MEK is essential downstream of BRAF-V600E in HT29,
  so knockouts DO move vs neutrals. Good reminder that the standard dropout-style QA gate is sometimes
  informative and sometimes null on resistance screens depending on target essentiality — the tool should
  not treat a null KS as a hard fail.
- [FEATURE REQUEST] The most valuable cross-target analysis here (paralog spatial conservation of the
  escape surface) required a hand-rolled Biopython superposition + ligand transfer. A built-in
  "map hotspots of gene A onto structure of paralog B (via alignment+superpose)" utility would be a natural
  BE-MetaClust3D companion for paralog/ortholog studies.

## BCL2 / MYC run (Queen "BCL2/MYC" round 2; Coelho essentiality / Control arms)
- [BLOCKER-for-purpose][data] **The shipped "Coelho = venetoclax/BCL2" assumption is wrong.** The MOESM4
  `ST2 BE z-scores` sheet has **no venetoclax/navitoclax arm** (drug arms: HT29 DebCet/Tram/Pict, H23
  Adag/Sotor, PC9 Osim/Gefit, MHHES1 Olap/Nirap). BCL2 (5,562 guides) and MYC (2,498) are tiled but the only
  target-relevant readout is **essentiality via the per-line Control (no-drug) arms**. **Lesson/fix:** a
  pre-flight "does an assay arm relevant to this target exist, and is the target essential/altered in the
  screened line?" check would flag un-runnable-for-purpose targets (BCL2: not essential in HT29/H23/PC9 →
  nonsense guides *enrich*, no functional window) before a full sweep is spent.
- [SCIENCE][calibration][HIGH] **Base rate / over-call on weak- or null-signal inputs.** On BCL2 (a genuine
  null: no venetoclax, non-essential) BE3D calls **~50% of scored residues** significant at p<0.05 and ~95% of
  "robust" (≥3/4-run) hotspots are false positives; groove enrichment OR<1. Even on MYC (real essentiality)
  the base rate is ~32%. The raw p<0.05 randomization threshold is **not a usable prioritizer** — it needs an
  **effect-size floor + FDR/BH correction + an explicit base-rate report** in `RUN_COMPLETED.txt` so a user
  sees "K sig / N scored = X%" and is warned when X is high (per the A0 spec). Magnitude/rank + the
  independent enrichment test are what actually discriminate.
- [SCIENCE][resolution] **Domain-level, not residue-level, on compact folds.** On MYC's bHLH-LZ, BE3D cleanly
  separates the folded functional domain from the IDR (OR 5.8–9.9, p→3e-7, 4× gap) BUT **within** the domain
  does not resolve DNA/MAX contacts from neighbours (within-domain enrichment ≈ chance) — a 6 Å neighborhood
  smears signal across a small helical bundle. Worth documenting as an expected limit; an adaptive/smaller
  radius or per-residue (not neighborhood-summed) channel might sharpen compact domains.
- [SCIENCE][value-add, positive] **BE3D recovers low-pLDDT functional boxes a folding/AlphaMissense baseline
  misses** — MYC MBI (60-65) and all of MBII (145-158) are robustly significant despite pLDDT 51-59. This is
  the genuine edge over "is it folded" (folded R/E=3.45 beats BE3D's ~2 on the *whole-protein* contact set, but
  cannot see the disordered boxes). Keep/showcase this behaviour.
- [FRICTION][numbering] **MYC 454-vs-439 frame trap.** `AF-P01106-F1-model_v6` and the Coelho screen are the
  **454-aa CUG-initiated MYC1** frame; literature/1NKP/AlphaMissense are **439-aa** (canonical-454 = 439 + 15).
  The screen matched AF at offset 0 (so BE3D I/O is self-consistent in 454), but every external ground-truth
  residue needed +15. A `structureid`↔`screen`↔`UniProt-canonical` frame reconciliation report would prevent
  silent 15-residue mis-annotation.
- [CONFIRM] `--qa-controls "No Mutation" Silent` parses the two-word token correctly (KS file header shows
  `Nonsense_Splice-donor_vs_No Mutation_Silent`); the sibling-`be3d_configs/` per-run-parent-dir workaround
  again avoided the shared-config race across an 8-run sweep.

## Round-2 PI3K/AKT run (PIK3CA P42336, AKT1 P31749; Coelho 2024 screen) — issues & findings
- [SCIENCE][base-editing reachability — the dominant limiter, again] The screen simply cannot install the
  #1 oncogenic hotspots. **PIK3CA H1047R/H1047L, M1043, G1049, H1065, C420, E418 have NO guide** (kinase +
  some C2 hotspots unreachable by CBE/ABE). **AKT1 E17K is not installed** — codon 17 only yields the
  non-activating **E17G** (ABE A>G), never oncogenic E17K (needs G>A); T308 also unreachable. BE3D is
  structurally blind to these by construction. A **reachability report** (per validated-hotspot: is the
  required codon change achievable by the editors in the library? which substitution IS installed?) would
  stop users mistaking an editor gap for a biological negative, and flag "wrong-substitution" traps like
  E17K→E17G. This is now the recurring #1 caveat across KRAS/MYC/PIK3CA/AKT1.
- [SCIENCE][base-rate blowup on compact folds, confirmed on AKT1] AKT1 (480 aa) p<0.01 base rate = 24%
  (r6), rising to 31% at r8; union enrichment < 1 (below chance) and the **discrimination gap inverts**
  (tolerant hit-rate 0.24 > functional 0.11). When the input screen is near-null (no on-target drug),
  6–8 Å neighborhood aggregation over a small fold paints ~25% of residues as "hotspots." The
  `RUN_COMPLETED.txt` should print base rate + a loud warning above ~15–20%, and ideally refuse to emit a
  "hotspot list" when union enrichment vs an independent functional set is not > 1.
- [SCIENCE][value-add is drug-DIRECTIONALITY, not recall] On PIK3CA, BE3D does NOT beat AlphaMissense on
  functional-residue recall (AM R/E 2.29 p 1e-5 vs BE3D union R/E 1.7 p 0.15, n.s.). Its genuine edge is
  that the **screen signal is drug-directional**: under trametinib (MEK bypass) the resistance hotspots
  concentrate on the **activating helical E542/E545 nSH2 interface (ACT R/E 2.9–3.6, p 0.008–0.02)** and
  the ATP/drug pocket is silent (POCKET R/E = 0), whereas under the on-target pictilisib there is NO signal
  anywhere. AlphaMissense scores all of these uniformly pathogenic and cannot make that distinction.
  Recommend BE3D lead its characterization with the **functional-subset-resolved** enrichment (activating
  vs pocket vs catalytic), not a single whole-protein number — that split is where it actually informs.
- [SCIENCE][mean >> sum on these screens] `function_for_lfc3d=sum` inflated base rate and DESTROYED the
  activating-hotspot enrichment on PIK3CA (r6_sum ACT R/E 1.8 vs mean 3.6; union dropped below 1). Mean is
  the safer default for aggregating multi-guide neighborhoods; sum over-weights dense-coverage regions.
- [DATA-SCHEMA] Reconfirmed: Coelho ST2 category lives in `most_severe_consequence` (not
  `variant_classification`); drug arms are shared across all tiled genes, so a target gene (PIK3CA/AKT1)
  is scored under drugs that don't target it — the informative arm must be chosen by pathway logic
  (PIK3CA: pictilisib on-target = null, trametinib bypass = signal). The brief's assumption that each gene
  has a matched inhibitor (alpelisib/capivasertib) is WRONG for this file — those drugs are absent.
- [CONFIRM] 16-run sweep (2 genes × 2 arms × 4 configs) all EXITCODE 0 / RUN_COMPLETED SUCCESS; the
  per-run-parent `be3d_configs/` workaround again prevented the shared-YAML race; DSSP placeholder fine.

## PARP1 / PARP2 run (round 2 — olaparib/niraparib resistance, Coelho 2024)
- [SCIENCE][direction matters: essentiality is recovered, resistance is not] On PARP1, BE3D's **negative
  (depletion)** direction cleanly re-finds the catalytic NAD+/inhibitor pocket + HD autoinhibitory helix
  (hypergeometric p ≈ 1e-10, all 12 pocket residues, beats burial AND AlphaMissense: R/E 6.1 vs 2.9 vs 2.0).
  But the **positive (resistance)** direction — the biologically interesting one — has base rate 10–28% and
  does NOT beat chance/AlphaMissense (olaparib R/E 0.7–1.1 n.s.; AlphaMissense R/E 2.1 p 5e-6). BE3D should
  report enrichment **per direction AND per functional-subset** (catalytic-depletion vs resistance-enrichment);
  a single whole-protein number hides that it wins on essentiality and loses on resistance specificity.
- [DATA-SCHEMA][PARP inhibitors present ≠ brief] The Coelho ST2 file contains ONLY **olaparib (Olap) and
  niraparib (Nirap)** for PARP1/PARP2 (cell line MHH-ES-1, **ABE only**, 2104/809 guides). There is **no
  talazoparib or rucaparib** arm — the brief's four-inhibitor comparison is not possible from this dataset.
  Any turnkey recipe should print available `L2FC_*zscore` arms per gene, not assume a fixed inhibitor panel.
- [DATA-SCHEMA][SILENT numbering offset — high-risk] PARP2 `Amino_Acid_Position` is numbered in a
  **canonical−13 short-isoform frame**; a naive canonical mapping gave only **14% refAA match** (looks like
  noise, not an offset — the per-offset histogram was flat at ±5). Only a wider offset scan revealed +13 →
  **86% match**. BE3D silently drops refAA-mismatched edits, so this would have quietly gutted PARP2 coverage
  with no error. **Recommend BE3D emit a refAA-match-rate report and auto-suggest an integer offset** when the
  match rate is low but a single shift recovers it (isoform-frame detection).
- [SCIENCE][coverage/power gates discrimination] PARP2 (ABE-only, 809 guides, ~200/583 residues scored) FAILS
  discrimination entirely (catalytic R/E ≈ 1, p 0.3; precision@10 = 0; loses to AlphaMissense), while PARP1
  (2104 guides) succeeds on the same pocket. Same tool, same pocket, different power. A **coverage/power
  warning** (e.g. "<X% of the domain tiled → discrimination unreliable") in RUN_COMPLETED would prevent
  over-reading under-powered targets.
- [SCIENCE][per-inhibitor 3D maps are breadth-different, not cleanly drug-specific] Olaparib vs niraparib
  positive maps: Jaccard 0.22 (PARP1) / 0.33 (PARP2). The contrast is niraparib producing a ~2× broader
  DNA-binding-centric hotspot, not an orthosteric-vs-distal split like KRAS sotorasib/adagrasib. High base
  rate blurs per-drug geometry; lower/adaptive thresholds may be needed to resolve drug-specific signatures.
- [CONFIRM] 16-run sweep (2 genes × 2 arms × 4 configs) all EXITCODE 0 / RUN_COMPLETED SUCCESS; PARP2 BE-QA
  KS is significant (D=0.28, p=1.3e-4) — a genuine pass; per-run-parent `be3d_configs/` workaround again
  avoided the shared-YAML race; DSSP placeholder fine (characterization-only).
