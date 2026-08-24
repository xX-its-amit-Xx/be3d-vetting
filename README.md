# BE3D Vetting & Real-Use Benchmark

An **independent** capability assessment of **[BE3D](https://github.com/broadinstitute/BE3D)** (`beclust3d`),
the Broad Institute's structure-function tool for base-editor tiling-mutagenesis screens — plus a set of
**real, benchmarked applications** of BE3D to public screens it had never seen.

**▶ Live dashboard:** https://be3d-vetting.vercel.app

> Not affiliated with or endorsed by the BE3D authors (Iqbal Lab, Broad · Liau Lab, Harvard).
> This is an outside analysis. Underlying screen data belongs to the original papers (cited); only the
> derived analysis, code, and text here are ours (MIT).

---

## What's here

1. **A competitive/《design》vetting** of BE3D (`research/Q1`–`Q5`): the field landscape and where BE3D sits
   (its one direct rival is **ProTiler-Mut**, Cell Systems 2026), 20 user personas + market, cross-domain
   critique, a code-level design deep-dive, and a real reproduction of BE3D's shipped KBTBD4 example.
2. **BE3D applied to real problems** (`research/P2_run_*`, `research/P2_user_*`): we installed BE3D and ran it
   on **five public base-editing screens** — KRAS, MCL1, DNMT3A, TP53, BRCA1 — each on a genuinely open 3D
   question, with an exhaustive parameter sweep and a full decision log.
3. **Honest discrimination benchmarks** (`research/P2_user_*`): for each target we test whether BE3D
   *discriminates* or merely *over-calls* — base rate, enrichment vs chance (hypergeometric/Fisher), the
   functional-vs-tolerant discrimination gap, precision@K, and whether BE3D beats a trivial baseline
   (burial / AlphaMissense). False positives and false negatives are flagged explicitly.
4. **Shareable outreach packages** (`outreach/<GENE>/`): a per-user email + 1-page brief + a hotspot table
   (with TP/novel/false-positive flags) + a [Genomics 2 Proteins Portal](https://g2p.broadinstitute.org)
   TSV for interactive 3D viewing.
5. **A reusable harness** (`harness/`) that runs BE3D on any UniProt + screen, and a running
   **[BE3D_IMPROVEMENTS.md](BE3D_IMPROVEMENTS.md)** log of real-usage bugs/gaps we hit and fixed.

## Headline honest findings

- BE3D's **base rate varies a lot** (≈13–16% of residues flagged on KRAS, but **36–45% on TP53/MCL1/DNMT3A**),
  so its raw p<0.05 flag can't be trusted uniformly — **you must rank**. This is the concrete face of its
  "no FDR/q-values" gap.
- On **ranked/enrichment** terms it genuinely discriminates and beats trivial baselines:

  | Target | Base rate | Enrichment vs chance | Discrimination gap (func vs tol) | Verdict |
  |---|---|---|---|---|
  | KRAS (adagrasib) | 15.9% | R/E 5.6×, OR 57, p=1.4e-6 | 45% vs 1.9% (~24×) | strong, validated (8/9 known resistance res) |
  | KRAS (sotorasib) | 12.7% | pocket NS (p=0.51) | 1.4× | **unproven** — distal α5 finding is a hypothesis |
  | TP53 | 36% | R/E 2.63×, p=4.2e-7 | 94% vs 14% | discriminates; **AlphaMissense out-ranks it on the conserved core** |
  | MCL1 (drug arms) | ~40–45% | OR ∞, p≈6e-4 | 100% vs 12.5% (8×) | strong; **only 3D map available (MCL1 isn't point-mutated)** |
  | DNMT3A | 37% | OR ~12, p≈1e-11 | 83% vs 24% | discriminates on rank; **R882 is base-edit-invisible (false negative)** |

- **Where BE3D wins:** functional-screen phenotype + effect direction + editor-accessibility, and
  **allosteric/resistance/interface** sites that recurrence tools can't see. Its clearest edge is on
  **druggable-but-not-recurrently-mutated** proteins (MCL1: somatic 3D-hotspot tools return *nothing*).
- **Where it loses:** raw over-calling (needs FDR); **AlphaMissense out-discriminates on conserved cores**
  (TP53); and it is blind to functionally critical residues base editing simply cannot reach.

See the dashboard and each `research/P2_user_<GENE>.md` for the full, cited analysis.

## Reproduce

Environment: Python ≥3.10 (tested on 3.14), a conda/venv with BE3D's deps. BE3D upstream:
`pip install git+https://github.com/broadinstitute/BE3D`. We applied small portability patches (see
[BE3D_PATCHES.md](BE3D_PATCHES.md)) — all logged in [BE3D_IMPROVEMENTS.md](BE3D_IMPROVEMENTS.md).

```bash
# 1. fetch a structure + generate a DSSP for any UniProt
python harness/prep_target.py <UNIPROT> <out_dir>
# 2. run BE3D on a reformatted screen TSV (cols: Gene, Mutation_list, Mutation_type, sgRNA_score)
python harness/run_new_target.py --screen <tsv> --gene <G> --uniprot <U> --out <dir> \
       --mut_list_col Mutation_list --mut_col Mutation_type --val_col sgRNA_score --gene_col Gene \
       [--screens a.tsv,b.tsv   # meta-aggregate]  [--function-for-lfc3d mean|sum] [--atom-level]
```

**Data** is *not redistributed* (it belongs to the source papers). Each `research/P2_user_<GENE>.md` and
`research/P2_scout_datasets.md` gives the exact download URL and the column mapping used. Figures under
`figures/` are BE3D's own SVG outputs from these runs; the dashboard (`site/index.html`) is self-contained.

## Layout

```
site/index.html            the dashboard (self-contained)
figures/                   real BE3D output figures (SVG)
research/Q1..Q5            competitive / design vetting dossiers
research/P2_scout_*        the public-dataset scouting + recipes
research/P2_run_<GENE>     per-target run decision logs (exhaustive sweeps)
research/P2_user_<GENE>    per-target discrimination benchmark + competitor comparison
research/P2_DECISION_LOG_TEMPLATE.md , P2_harness.md , OUTREACH_GUIDE.md
outreach/<GENE>/           copyable email + brief + hotspot TSV (flagged) + G2P TSV
harness/                   prep_target.py, run_new_target.py, gen_dssp.py
BE3D_IMPROVEMENTS.md       real-usage bug/gap log with concrete fixes
```

## Credits & disclaimer

BE3D is by the **Iqbal Lab (Broad)** and **Liau Lab (Harvard)** (Hu, Myung, Mani, Iqbal); it pairs with the
**Genomics 2 Proteins Portal** (Kwon et al., *Nat Methods* 2024). Screen data: Coelho/Dincer 2024
(*Nat Genet*), Hanna 2021 (*Cell*), Lue 2023 (*Nat Chem Biol*), Sánchez-Rivera 2022 / MaveDB (TP53),
Cuella-Martin 2021 (*Cell*) — see each decision log for DOIs. This repository is an **independent** analysis;
novel hotspots are **hypotheses**, not validated claims. Code & prose: MIT (`LICENSE`).
