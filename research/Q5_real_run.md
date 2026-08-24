# Q5 — BE3D real run on KBTBD4 (Yeo 2025 ABE + CBE screens)

**Date:** 2026-08-23
**Result:** Full pipeline ran end-to-end (`EXITCODE=0`) producing real LFC3D per-residue
scores, per-threshold significance, meta-aggregation, and 3D hotspot clusters at 6 Å.

- Config used: `examples/yaml/KBTBD4_chain_B_run.yaml` (copy of `KBTBD4_chain_B.yaml`
  with `output_dir: ./output/KBTBD4_run` and `user_dssp` set).
- Output dir: `.../scratchpad/BE3D/examples/output/KBTBD4_run` (332 files).
- Key outputs copied to: `D:\Users\ashenoy00000\.cursor\BE3D Vetter\real_output\`.

---

## Which fallback worked: **Fallback 1 — user-supplied DSSP** (plus 3 small code fixes)

The confirmed blocker was `run_dssp()` shelling out to `mkdssp` (not installed →
`FileNotFoundError [WinError 2]`). Providing a `user_dssp` file makes
`sequence_structural_features` **skip** the `mkdssp` call entirely.

`DSSPparser.parseDSSP` reads fixed-width columns; `parse_dssp` only needs, per chain-B
residue: `inscode` (= uniprot/PDB residue number), `chain`, `aa`, `struct` (SS; may be
coil placeholder — only used in the final characterization step), and float-parseable
`acc`, `phi`, `psi`. I generated a valid classic-DSSP file from the processed PDB with
Biopython (`gen_dssp.py`), writing 1374 residues (chain B: 504, residues 26–534) with
placeholder SS (coil), `acc=100`, `phi=-60`, `psi=-45`. It parsed cleanly and chain-B
`inscode`s matched the UniProt positions in the coord table.

Generator: `scratchpad/gen_dssp.py`; output DSSP:
`examples/pdb/PDB-Q9NVX7-2_processed.user.dssp`.

### Three additional blockers fixed (all pandas-3.x / Windows incompatibilities, minimal edits)

**1. `to_csv` positional `sep` no longer allowed** (pandas removed positional args).
`beclust3d/aggregate/nonaggregate.py:303` and `beclust3d/aggregate/metaaggregate.py:355`:

```python
# before
df_z.to_csv(filename, "\t", index=False)
df_meta_Z.to_csv(filename, "\t", index=False)
# after
df_z.to_csv(filename, sep="\t", index=False)
df_meta_Z.to_csv(filename, sep="\t", index=False)
```

**2. Windows-illegal `<` in output filenames.** Significance labels `p<0.05` etc. are
used BOTH as DataFrame values (compared downstream) AND baked into filenames. I left the
data labels intact and sanitized only the filename usages in
`beclust3d/lfc3d/clustering_plot.py` (`<` → `lt`, `>=` → `ge`):

```python
# in the `for gid, gcont in cluster_dist_pd.groupby('pthr'):` loop
gid_safe = str(gid).replace('<', 'lt').replace('>=', 'ge')
clust_filename = ... f"..._{gid_safe}_Aggr_Hits_List.tsv"
plot_filename  = ... f"..._{gid_safe}_cluster_distance.{save_type}"

# in the `for name, pthr, colname in zip(...):` loop
pthr_safe = str(pthr).replace('<', 'lt').replace('>=', 'ge')
# dend_filename / analysis_color_file / txt_filename now use {pthr_safe}
```
(Cluster filenames therefore read e.g. `..._Dendrogram_plt0.05_6A.txt`.)

**3. `Series[0]` positional indexing** (now label-based → `KeyError: 0`) in the runner
`examples/be3d_local.py` `find_union` (line ~583):

```python
def find_union(input, pthr_str):
    input = list(input)  # positional access (pandas Series[0] is now label-based)
    if input[0] == f'p<{pthr_str}' or input[1] == f'p<{pthr_str}':
        ...
```

None of these edits change scientific logic; they are API/portability fixes.

---

## Results — actual numbers

### LFC3D per-residue (NonAggr, `LFC3D/KBTBD4_NonAggr_LFC3D.tsv`, 534 residue rows)

Residues with non-zero LFC3D: **ABE = 504**, **CBE = 501** (chain B, residues 26–534).

Significant residues (per screen, LFC3D z-based p):

| screen | dir | p<0.05 | p<0.01 | p<0.001 |
|---|---|---|---|---|
| ABE | neg | 109 | 86 | 78 |
| ABE | pos | 114 | 90 | 69 |
| CBE | neg | 121 | 82 | 59 |
| CBE | pos | 76 | 73 | 61 |

### Meta-aggregate (SUM across both screens, `meta-aggregate/KBTBD4_MetaAggr_LFC3D.tsv`)

| dir | p<0.05 | p<0.01 | p<0.001 |
|---|---|---|---|
| neg | 145 | 131 | 117 |
| pos | 160 | 146 | 122 |

**Top 12 residues by |meta SUM LFC3D| (all positive-direction, all significant at p<0.001):**

| unipos | aa | direction | SUM LFC3D |
|---|---|---|---|
| 93 | N | pos | 3.773 |
| 94 | L | pos | 3.670 |
| 76 | R | pos | 3.184 |
| 88 | S | pos | 2.996 |
| 92 | S | pos | 2.956 |
| 89 | M | pos | 2.900 |
| 90 | F | pos | 2.900 |
| 77 | L | pos | 2.866 |
| 91 | T | pos | 2.746 |
| 78 | V | pos | 2.521 |
| 80 | S | pos | 2.510 |
| 95 | K | pos | 2.439 |

The dominant per-residue signal is a tight positive-LFC (gain-of-function /
enrichment-direction) hotspot spanning **residues ~76–95** on chain B.

### 3D hotspot clusters (meta LFC3D, 6 Å radius; total / active i.e. size≥2)

| file | total | active (size≥2) |
|---|---|---|
| Meta Positive p<0.05 | 17 | 15 |
| Meta Positive p<0.01 | 22 | 17 |
| Meta Positive p<0.001 | 18 | 14 |
| Meta Negative p<0.05 | 19 | 11 |
| Meta Negative p<0.01 | 18 | 10 |
| Meta Negative p<0.001 | 15 | 9 |

**Largest positive cluster (p<0.05, 6 Å):** Cluster 3.0, **58 residues** —
B-52…62, B-74…99, B-134…154 — i.e. the 76–95 peak fused with neighboring surface loops.
Other notable positive clusters: Cluster 1.0 (20 res: B-388,389,390,397,398,399,
416,417,438–441,450,451,458–460,467–469) and Cluster 0.0 (13 res incl. B-296…305,
B-334…339, B-353, B-385,386).

**Largest negative clusters (p<0.05, 6 Å):** Cluster 1.0 (**46 res**, B-262…266 +
B-453…455 + B-481…530) and Cluster 0.0 (**38 res**, B-180…247 core).

### KBTBD4 biology note
KBTBD4 is a CUL3 BTB-Kelch substrate receptor. The recurrent medulloblastoma
gain-of-function in-frame indels sit in the Kelch β-propeller around **R379–R390**.
The positive-direction meta cluster **1.0 includes B-388/389/390** (and 397–399),
overlapping that functional indel hotspot; cluster 0.0 picks up B-385/386 as well. The
strongest raw per-residue signal (76–95) lies N-terminal of the Kelch repeats.

---

## Caveats
- **Secondary structure is placeholder (coil).** The user_dssp has SS9 blank and
  `acc/phi/psi` set to constants, so only the *characterization* step's SS/burial
  annotations are non-physical. LFC3D, significance, meta-aggregation, and 3D clustering
  (which use only Cα coordinates + LFC) are fully real and unaffected.
- `nRandom=500` (as configured); p-values are from the z-normalized empirical null.
- Isoform Q9NVX7-2, chain B of local PDB 8voj; monomer mode. Directions: `pos` =
  positive LFC (enrichment), `neg` = negative LFC (depletion/essential).

## Files copied to `real_output\`
`LFC3D/KBTBD4_NonAggr_LFC3D.tsv`, `LFC3D/KBTBD4_LFC3D_bidirectional.tsv`,
`meta-aggregate/KBTBD4_MetaAggr_LFC3D.tsv`, `meta-aggregate/KBTBD4_MetaAggr_LFC.tsv`,
all `cluster_LFC3D/*_6A.txt` cluster-membership files + `*_Aggr_Hits.tsv`,
`g2p_visualization/*.tsv` (G2P-formatted), and 5 representative `.svg` plots.
