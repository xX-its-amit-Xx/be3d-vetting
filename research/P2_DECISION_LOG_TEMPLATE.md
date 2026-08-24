# BE3D Run — Decision Log TEMPLATE  (copy per target → P2_run_<GENE>.md)

> Every runner-worker MUST fill this out. The point is a transparent, reproducible chain from
> raw public data → BE3D parameters → results → NOVEL, actionable proposals. No unexplained choices.

## 0. Target & the open gap
- Gene / protein: … (UniProt: …, isoform: …)
- Structure used: AlphaFold (AF-…-model_vN) | PDB … (chain …, liganded? complex/PPI?)
- Source dataset: paper (author, year, DOI) + exact data URL(s) fetched
- **The OPEN question we attack** (NOT the paper's headline): …
- **What NEW actionable output we intend to produce**: …

## 1. Data provenance & reformatting  (decision log)
| Decision | Choice | Why | Alternative rejected |
|---|---|---|---|
| Source file(s) fetched | URL + sha/size | verified downloadable | … |
| Rows kept / filtered | e.g. dropped controls w/o LFC | … | … |
| Mutation_list column | mapped from `<col>` | … | … |
| Mutation_type column | mapped from `<col>` (+ mapping of categories) | … | … |
| sgRNA_score (LFC) column | mapped from `<col>` | which replicate/condition & why | … |
| Multi-edit guides | how bystander/multi-AA edits handled | … | … |
| Screens combined | which screens → meta, and why | … | … |
- Reformatting script: path + one-line description. Row counts in/out. Sanity checks run.

## 2. QA (BE-QA)  (decision log)
- cases = […]  controls = […]  (why these categories)
- H1/H2 result (MWU + KS p-values); **screen accept/reject decision** and rationale.

## 3. Parameter sweep grid  (EXHAUSTIVE — this is the mandate)
Record every run. Vary and justify each:
| run_id | structure_radius | clustering_radius | function_for_lfc | function_for_lfc3d | function_for_meta | nRandom | p-thr | atom_level | scope (monomer/PPI) | notes |
|---|---|---|---|---|---|---|---|---|---|---|
| r01 | 6 | 6 | mean | mean | SUM | 1000 | .05/.01/.001 | F | monomer | baseline |
| r02 | 8 | 6 | … | … | … | … | … | … | … | radius sensitivity |
| … | | | | | | | | | | |
- Radius sweep (e.g. 4/6/8/10 Å): how do hotspots move? (robustness)
- Aggregation (mean vs sum vs median): effect on ranking.
- nRandom (e.g. 500 vs 2000): tail-p stability at p<0.001.
- atom-level vs residue-level neighbors: interface sensitivity.
- meta strategy if >1 screen (SUM vs MEAN): consensus vs additivity.

## 4. Results per run
- # residues scored; # significant (neg/pos) at each p-thr; # 3D clusters (active) at chosen radius.
- Top hotspot residues (position, AA, LFC3D, direction, significance).
- Cluster membership for the leading hotspots (residue ranges).
- Robustness: which hotspots persist across the parameter sweep? (report a stability table)
- Figures saved (paths).

## 5. Interpretation & NOVELTY check
- Map top hotspots to known biology (domains, active/binding sites, PPI interfaces, ligand pockets, known disease mutations via UniProt/ClinVar/cBioPortal).
- **Which hotspots are NOT explained by the source paper's stated findings?** (the novel part)
- Confidence caveats: pLDDT of hotspot residues; disordered regions; structure choice effects.

## 6. NEW actionable proposals for other groups (the deliverable)
- Concrete, testable next steps (e.g. "mutate residues X/Y/Z — predicted allosteric patch not previously highlighted"; "residue cluster at the A–B interface suggests a PPI-disruption handle"; "resistance hotspot spatially distinct from the active site → candidate allosteric drug site").
- Who should act (which community/lab type) and what experiment/analysis.

## 7. BE3D issues encountered → append to BE3D_IMPROVEMENTS.md
- List anything that broke, was awkward, or affected correctness, with the fix/workaround.
