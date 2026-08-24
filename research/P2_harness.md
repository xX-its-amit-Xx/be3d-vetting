# P2 — Plug-and-play BE3D harness for arbitrary targets

Two small scripts let any worker run BE3D on a **new** target (arbitrary UniProt +
a reformatted screen TSV) without rediscovering setup issues. Verified end-to-end
on **LYZ / UniProt P61626** (human lysozyme C, 148 aa) — a target unrelated to the
KBTBD4 / HDAC1 / MORC2 examples.

Scripts live in the scratchpad next to the BE3D clone and the venv:

```
scratchpad/
  prep_target.py       # download AF model + generate DSSP
  run_new_target.py    # write monomer YAML + invoke be3d_local.py
  gen_dssp.py          # (pre-existing) classic-DSSP generator, reused by prep_target
  BE3D/                # editable-installed clone (examples/be3d_local.py)
  be3dvenv/            # venv (py3.14); python = be3dvenv/Scripts/python.exe
```

Always call them with the venv interpreter:
`C:\...\scratchpad\be3dvenv\Scripts\python.exe`.

---

## 1. `prep_target.py` — fetch structure + DSSP

Downloads the AlphaFold model (`AF-{uniprot}-F1-model_{v6|v4|v2}.pdb`, trying v6→v4→v2)
and generates a matching classic-DSSP (so BE3D never needs the missing `mkdssp` binary).

```
python prep_target.py --uniprot P61626 --gene LYZ --outdir ./targets/LYZ
```

- `--uniprot` : **BASE** accession only (e.g. `P61626`). AlphaFold DB has no isoform
  files, so an isoform id like `Q9NVX7-2` will 404 — for isoform numbering supply your own PDB.
- `--outdir`  : where the `.pdb` and `.dssp` are written.
- `--gene`    : optional, cosmetic (used in the summary line).

Prints machine-parseable paths:
```
PDB_PATH=...\targets\LYZ\AF-P61626-F1-model_v6.pdb
DSSP_PATH=...\targets\LYZ\AF-P61626-F1-model_v6.dssp
AF_VERSION=v6
```

The DSSP is a placeholder (coil SS, filler ASA/torsions); it only feeds the final
*characterization* step and never affects LFC3D / significance / clustering.

---

## 2. `run_new_target.py` — write YAML + run BE3D (monomer)

Writes a self-contained monomer YAML (conservation OFF, PPI dict null,
`user_pdb`/`user_dssp` wired in) and invokes `examples/be3d_local.py` on it.

```
python run_new_target.py \
    --screen ./plumbing_test_LYZ_P61626.tsv \
    --gene LYZ --uniprot P61626 \
    --pdb  ./targets/LYZ/AF-P61626-F1-model_v6.pdb \
    --dssp ./targets/LYZ/AF-P61626-F1-model_v6.dssp \
    --outdir ./output/LYZ_smoke \
    --mut-list-col Mutation_list --mut-col Mutation_type \
    --val-col sgRNA_score --gene-col Gene \
    --nrandom 100 --qa-cases Nonsense --qa-controls "No Mutation"
```

Key flags:
- `--screen`  : one TSV, or comma-separated list for **meta-aggregation** across screens.
- column map  : `--mut-list-col` (edits, e.g. `M1V;K2G;`), `--mut-col` (category,
  e.g. `Missense;`), `--val-col` (numeric score), `--gene-col` (rows filtered to `--gene`).
- tunables    : `--nrandom` (default 500), `--structure-radius`/`--clustering-radius`
  (default 6.0), `--single-screen-pthr`/`--multi-screen-pthr`, `--chain` (default `A`),
  `--qa-cases`/`--qa-controls` (accept multiple values).
- `--dry-run` : write the YAML but don't run.

The generated YAML is written to a sibling `be3d_configs/` dir (NOT inside `--outdir`),
because `be3d_local.py` copies the config into `output_dir` and would otherwise hit a
Windows `WinError 32` (src==dst). Screens are staged into `<outdir>/screens/`.

### Screen TSV format
Tab-separated, with the four named columns. Edits are `refAA + pos + altAA`, `;`-joined
(`M1V;`), `*` for a stop (`R21*;`); Silent = `refAA==altAA` (`A6A;`); `No Mutation`
controls use non-coding tokens (`utr;`). Positions/refAA must match the structure sequence.

---

## Verified smoke test (LYZ / P61626)

Fabricated plumbing-test screen `plumbing_test_LYZ_P61626.tsv` (refAA drawn from the real
P61626 sequence): 40 Missense + 20 Silent + 12 Nonsense + 25 No-Mutation guides.

- Single screen → **EXITCODE 0**, `RUN_COMPLETED.txt: status: SUCCESS`.
- Two screens (`a.tsv,b.tsv`) → **EXITCODE 0**, `meta-aggregate/` with
  `LYZ_MetaAggr_LFC.tsv` + `LYZ_MetaAggr_LFC3D.tsv` produced.

Output tree under `./output/LYZ_smoke/` (single-screen run):

```
LYZ_smoke/
├── RUN_COMPLETED.txt            # status: SUCCESS
├── LYZ_run.yaml                 # archived config (copied by be3d_local.py)
├── sequence_structure/          # processed PDB, coord + struc features, DSSP parsed (10 files)
├── screendata/  screendata_rand/
├── screendata_sequence/  screendata_sequence_rand/
├── hypothesis_qa/               # QA / KS tests (12 files)
├── LFC/                         # LFC bidirectional / NonAggr / dis_wght (6 files)
├── LFC3D/                       # LFC3D NonAggr + randomized LFC3Dr (8 files)
├── cluster_LFC/                 # LFC clusters + dendrograms (10 files)
├── cluster_LFC3D/               # LFC3D clusters + dendrograms (10 files)
├── cluster_union/               # union hits (14 files)
├── characterization/            # enrichment tests + plots (4 files)
└── g2p_visualization/           # gene-to-protein hit/cluster export (3 files)
```

Multi-screen runs additionally produce `meta-aggregate/`.

---

## Blockers found + fixed (also logged in BE3D_IMPROVEMENTS.md)

1. **[BLOCKER] Config copy src==dst.** `be3d_local.py` copies the config into `output_dir`;
   if the source YAML lives there too, Windows raises `WinError 32`. Harness fix: write the
   YAML to a sibling `be3d_configs/` dir.
2. **[FRICTION] Isoform accessions have no AF file.** Pass the base UniProt to `prep_target.py`;
   use a user PDB for isoform numbering.
3. **[FRICTION] DSSP from raw vs processed PDB.** Fine for single-chain AF F1 models (residue
   numbers == unipos); DSSP only affects characterization regardless.
