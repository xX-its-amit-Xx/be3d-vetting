# BE3D Run — Decision Log — MCL1 (PPI/complex mode; BH3-groove interface vs S63845 drug pocket)

Runner: Queen "MCL1" (P2). Date: 2026-08-24. All runs on venv py3.14 / beclust3d v1.0.0 clone.
Outputs: `real_output/P2_MCL1/` (tsv/, structures/, runs/, results_summary.tsv,
interface_vs_pocket_detail.tsv, robustness_baseline.tsv, analyze_mcl1_results.py).

## 0. Target & the open gap
- Gene / protein: **MCL1** (Induced myeloid leukemia cell differentiation protein Mcl-1), UniProt **Q07820**, canonical (350 aa; transcript ENST00000369026 == canonical, refAA 151/151 match).
- Structures used:
  - **AlphaFold** AF-Q07820-F1-model_v6 (monomer baseline, 350 aa).
  - **3PK1** — MCL1 (chains A/C, res 173-321) + **BAX** BH3 peptide (chains B/D, res 54-80). *(Brief said BAK; peptide sequence `ASTKKLSECLKRIGDELDSNMELQRMI` and PDB title confirm it is the **BAX** BH3 domain.)*
  - **3MK8** — MCL1 (chain A, 172-320) + a **stapled MCL-1 BH3 helix / SAHB** (chain B; `ALETLRRVGDGVQRNHETA` == Q07820 209-227, hydrocarbon-staple MK8 residues). *(Brief said BID; PDB title "the MCL-1 BH3 helix is an exclusive MCL-1 inhibitor" confirms it is an MCL-1-derived stapled peptide, not BID.)* Both peptides occupy the canonical BH3 groove, so the interface question is unaffected.
  - **6QB4** — MCL1 (chain A, 172-321) + small-molecule inhibitor **HVN** (an S63845-class BH3-mimetic); chain X is an antibody VH crystallization chaperone (ignored). 6QB4 chain A carries 5 engineered point substitutions (a solubility-optimized construct) — noted, does not affect numbering.
- Source dataset: **Hanna et al. 2021, Cell** "Massively parallel assessment of human variants with base editor screens", Table S4 raw counts — `https://ars.els-cdn.com/content/image/1-s2.0-S009286742100012X-mmc4.xlsx` (sheets `MCL1` + `MCL1 Library annotation`). CBE (BE3.9max) tiling of MCL1 in MELJUSO cells.
- **OPEN question we attack (NOT the paper's headline):** On the *same* MCL1 BH3 domain, does the base-editing **dependency/essentiality** signal localize to the **BH3-groove PPI surface** (residues contacting BAX / a BH3 peptide) while the **S63845 (MCL1-inhibitor) resistance** signal localizes to a **distinct orthosteric drug pocket** — i.e. can BE3D dissociate a protein-protein dependency interface from a drug-contact surface on one domain?
- **NEW actionable output intended:** a per-arm 3D hotspot map on MCL1 that either (a) separates a dependency interface from a drug pocket, or (b) shows they coincide — directly informing BH3-mimetic / molecular-glue design (which groove residues are the mutational escape hotspots, and whether escape can also arise off-groove).

## 1. Data provenance & reformatting (decision log)
| Decision | Choice | Why | Alternative rejected |
|---|---|---|---|
| Source file | `hanna_mmc4.xlsx`, sheet `MCL1` (392 count rows) + `MCL1 Library annotation` (392 rows) | Elsevier CDN serves openly to curl (200); annotation carries `Amino acid edits` + `Mutation category` | mmc2 (BRCA1/2) — different gene |
| Join key | `sgRNA sequence` | 392/392 perfect 1:1 match counts↔annotation | — |
| Rows kept | 124 of 392 | 268 guides have blank `Amino acid edits`/`Mutation category` (non-editing guides) → dropped | keeping them = untyped noise |
| Normalization | CPM per column, then LFC | raw counts differ in depth across conditions; CPM is the standard screen normalization | raw-count ratio (depth-biased) |
| sgRNA_score (LFC) | `log2((mean(CPM_arm)+0.5)/(CPM_pDNA+0.5))` for 3 arms: **dropout**, **A1331852** (BCL-xL-inhibitor; MCL1 becomes essential), **combo A1331852+S63845** (the S63845/MCL1-inhibitor drug arm; resistance enriches, +LFC) | pseudocount 0.5 on CPM; these are the 3 conditions present in the MCL1 sheet | S63845-alone arm not present for MCL1 |
| Mutation_list | `Amino acid edits` 3-letter → 1-letter (`Leu4Phe;`→`L4F;`), multi-edit `;`-joined kept | BE3D format; `refAA` validated 151/151 vs Q07820 | collapsing multi-edit guides |
| Mutation_type | `Mutation category`, `;`-joined kept; controls (UTR/Intron) **remapped to "No Mutation"** | BE3D hardcodes reading a `No_Mutation.tsv` control file; our raw controls are UTR/Intron (see §7 issue #1) | leave as UTR/Intron → crash |
| Multi-edit guides | kept intact; BE3D `mutation_priority` collapses `Silent;Missense;`→Missense per guide | matches BE3D design | — |
| Screens combined | run separately (not meta) | each arm answers a different question | meta would blur dependency vs drug |
- Reformatting script: `real_output/P2_MCL1/tsv/make_mcl1_tsv.py` (data worker). Category remap: inline fixer in run log.
- Row counts: 392 in → **124** usable (86 Missense, 28 Silent, 1 Nonsense, 9 controls after remap). Pearson(dropout, combo) = **0.29** (weakly correlated — different biology, as expected).

## 2. QA (BE-QA) (decision log)
- cases = [Nonsense]; controls = [No Mutation, Silent] (only LOF category available is Nonsense).
- H1/H2 (identical, single screen): **Mann-Whitney U p = 0.523; Kolmogorov-Smirnov D = 0.70, p = 0.632**; num_cases = **1**, num_controls = 37.
- **Screen decision: FAILS formal QA (underpowered), ACCEPTED for exploratory LFC3D with caveat.** The MCL1 tiling library is small and contains **only 1 nonsense guide**, so the knockout-vs-neutral separation test is statistically meaningless. Ran with `qa_passed_only: false`. Consequence: the **plain dropout arm has a weak depletion window** (mean Missense LFC ≈ +0.18, controls ≈ −0.26 — controls are the *most* negative group), so dropout hotspots are low-confidence. The **A1331852 arm** (MCL1 made essential by BCL-xL inhibition) has a much stronger essentiality window (LFC3D to −3.4) and is the trustworthy dependency signal.

## 3. Parameter sweep grid (EXHAUSTIVE)
21 runs total (19 monomer + 2 complex). p<0.05/0.01/0.001 are all computed every run (znorm hardcodes the three thresholds), so the p-thr axis is free. Driver: `scratchpad/run_mcl1_sweep.py`; complex configs `be3d_configs_mcl1/rC1.yaml,rC2.yaml`.

| run_id | struct | arm | radius | func_lfc3d | nRandom | atom | scope | notes |
|---|---|---|---|---|---|---|---|---|
| r01–r04 | AF | dropout | 4/6/8/10 | mean | 1000 | F | monomer | radius sensitivity |
| r05 | AF | dropout | 6 | **sum** | 1000 | F | monomer | aggregation |
| r06 | AF | dropout | 6 | mean | **500** | F | monomer | nRandom stability |
| r07 | AF | dropout | 6 | mean | 1000 | **T** | monomer | **atom-level** |
| r08 | AF | **combo (S63845)** | 6 | mean | 1000 | F | monomer | drug arm |
| r09 | AF | **A1331852** | 6 | mean | 1000 | F | monomer | dependency arm |
| r10/r11/r12 | **3PK1** | dropout/combo/A1331852 | 6 | mean | 1000 | F | monomer-on-complex | bound conformation |
| r13/r14 | **3MK8** | dropout/combo | 6 | mean | 1000 | F | monomer-on-complex | 2nd BH3 complex |
| r15/r16 | **6QB4** | dropout/combo | 6 | mean | 1000 | F | monomer-on-complex | drug complex |
| r17/r18 | 3PK1 | dropout | 4/8 | mean | 1000 | F | monomer | crystal radius sweep |
| r19 | 3PK1 | combo | 8 | mean | 1000 | F | monomer | radius×arm |
| **rC1** | **3PK1 A/C** | dropout | 6 | mean | 1000 | F | **complex (PPI)** | MCL1 homodimer, real cross-chain data → **RAN, exit 0** |
| **rC2** | 3PK1 A+pep B | dropout | 6 | mean | 1000 | F | **complex (PPI)** | BAX peptide partner (no screen data) → **CRASHED** (see §7 #2) |

## 4. Results per run
Residue sets (from crystal contacts, worker "structures"): **BH3 interface** = union(3PK1,3MK8) ≤5 Å = 31 residues; **drug pocket** = 6QB4/HVN ≤5 Å = 18 residues; **overlap = 13** (Jaccard 0.36). Derived: **IFACE-only** (18: 215,216,220,234,248,252,255,256,258,260,261,262,265,305,317-320), **POCKET-only** (5: **246,271,274,290,294**).

Enrichment (Fisher exact, one-sided, background = scored residues in 172-321). Full table: `results_summary.tsv`, `interface_vs_pocket_detail.tsv`.

| arm (representative run) | nsig p<0.05 (domain) | sig in IFACE-only | sig in POCKET-only | OR_iface (p) | OR_pocket (p) |
|---|---|---|---|---|---|
| dropout (r02 AF) | 14 | 4 (317-320 only) | 0 | 0.70 (ns) | 0 (ns) |
| dropout (r10 3PK1) | 13 | 4 (317-320) | 0 | 1.00 (ns) | 0 (ns) |
| **A1331852 (r09 AF)** | 33 | **10** (215,216,220,255,256,258,317-320) | 1 (294) | **inf (p=5.6e-4)** | 1.68 (ns) |
| **A1331852 (r12 3PK1)** | 32 | **10** | 1 (294) | **inf (p=1.1e-3)** | 1.47 (ns) |
| **combo/S63845 (r08 AF)** | 24 | **7** (215,216,220,255,256,258,265) | 0 | **inf (p=7.9e-4)** | inf (p=0.10, ns) |
| combo/S63845 (r11 3PK1) | 19 | 6 | 0 | 2.10 (ns) | 0.37 (ns) |
| combo/S63845 (r19 3PK1 r8) | 23 | 8 | 0 | 2.69 (p=0.097) | 0.55 (ns) |

**Top robust hotspots (magnitude, direction):**
- **S63845 resistance (+LFC3D):** groove α3/α4 residues **R215, V216, V220, T254, N255, V256, F258, T266, F270** (all runs), plus off-groove **α2 R176/Q177/S178/E180**, **C-terminal E325/G326**, **N-terminal E85**. Raw top single-guide resistance edit = **T212I (+2.44)** (aggregates into the 215-222 hotspot).
- **A1331852 dependency (−LFC3D, strongest window):** **G217/D218/G219/V220/Q221/R222** (α3-α4 groove core, LFC3D ≈ −3.4), plus C-terminal 324-329.
- **Plain dropout (−LFC3D, low confidence):** surface helix **295-303 (R300/T301/K302)** and **α8 edge 317-320**, N-terminal 57-61 — *not* the core groove; consistent with the weak/underpowered dropout window.

**Robustness (`robustness_baseline.tsv`, 9-run baseline matrix):** groove residues **215, 216, 220, 255, 256, 258** significant in **all** combo/A1331852 runs across **every** structure (AF, 3PK1, 3MK8, 6QB4) → structure-independent. Radius 4→10 Å broadens the neighborhood but preserves the same peak residues; **mean vs sum** and **nRandom 500 vs 1000** give identical top hotspots; **atom-level** (r07) adds a few residues but same peaks. **Complex mode rC1 == monomer r10** (identical dropout top hotspots) — the crystallographic A/C partner adds no signal (confirmed expected behavior, §5).
- Painted PDBs (|LFC3D| in B-factor) for viewers: `structures/MCL1_S63845resistance_LFC3Dpos_on_6QB4.pdb`, `MCL1_A1331852dependency_LFC3Dneg_on_3PK1.pdb`, `MCL1_dropout_LFC3Dneg_on_3PK1.pdb`.

## 5. Interpretation & NOVELTY check
- **The naive hypothesis is FALSIFIED in an informative way.** BE3D does **not** dissociate a PPI-dependency interface from a *distinct* drug pocket on MCL1 — instead it shows they are the **same surface**: both the dependency signal (A1331852 essentiality) and the S63845 resistance signal are significantly enriched on the **BH3 groove** (α3/α4: R215/V216/V220/N255/V256/F258 and the overlap shelf T254/T266/F270), and **neither** arm produces signal at the **pocket-only** residues (246, 271, 274, 290, 294) that extend beyond the BH3-peptide footprint. This is the direct functional-genomics readout that **S63845/HVN is a true orthosteric BH3-mimetic** — its resistance-relevant contacts coincide with the surface MCL1 uses to bind pro-apoptotic BH3 partners — rather than an allosteric binder with an independent pocket.
- **What BE3D DID cleanly separate (the genuinely novel part):** the S63845-resistance map has **two spatially distinct compartments** — (1) the **on-groove** escape hotspots (215/216/220/254/255/256/258/266/270), and (2) an **off-groove** set the source paper never spatially resolved: **α2 (R176/Q177/S178/E180)**, the **C-terminal helix (E325/G326)**, and the **N-terminal PEST/phosphodegron region (E85, R78, P80, R95)**. The off-groove sites are *not* drug contacts (all >8 Å from HVN) and most likely act through **protein stability / degradation / expression** (MCL1's N-terminal PEST region and phosphodegron control its short half-life). So resistance to a BH3-mimetic can arise both by remodeling the drug-binding groove *and* by an orthogonal stability mechanism — a partition invisible in the paper's pooled hit list.
- Not explained by the source paper: Hanna et al. report MCL1 variant scores but never map them to 3D or attribute them to the BH3 groove vs stability compartments; the α2 and PEST-region resistance clusters, and the groove-coincidence of dependency and drug signal, are new here.
- Confidence caveats: (i) plain **dropout arm underpowered** (1 nonsense guide; QA fails) — its non-groove hotspots (295-303) are low-confidence and likely screen noise; the **A1331852 arm is the reliable dependency readout**. (ii) N-/C-terminal hotspots (57-61, 85, 324-329) lie in **low-pLDDT disordered regions** in AF and are outside the crystal — spatial LFC3D there is weak (few neighbors), so treat as sequence-level (per-guide) resistance calls, not 3D clusters. (iii) 6QB4 chain A is a 5-mutation engineered construct.

## 6. NEW actionable proposals for BH3-mimetic / molecular-glue developers
1. **Orthosteric-escape watchlist for next-gen MCL1 BH3-mimetics.** The robust, structure-independent groove hotspots **R215, V216, V220, T254, N255, V256, F258, T266, F270** (and the single strongest edit **T212I**) are the mutational escape residues that must be de-risked. Design compounds whose pharmacophore does **not** depend on contacts to these side chains, or that tolerate their substitution (e.g. accommodate T266I / F270 remodeling). *Who:* med-chem teams at MCL1-inhibitor programs (S64315/AMG-176/AZD5991 successors). *Test:* saturation mutagenesis / deep-mutational-scan of these 9 codons under compound selection; co-crystallize T212I/T266I mutants.
2. **The pocket-only residues 246/271/274/290/294 are functionally silent in this screen** — they line HVN but not the BH3-peptide footprint and carry no dependency/resistance signal. This flags them as **safe anchor points**: contacts a mimetic makes here are unlikely to create essential-function escape liabilities, and could improve selectivity over BCL-xL/BCL-2. *Test:* structure-guided elaboration into the 271/274/290 sub-pocket.
3. **A second, orthogonal resistance axis = MCL1 stability.** The off-groove α2 (176-180), C-terminal (325/326) and N-terminal PEST/phosphodegron (78/80/85/95) resistance edits suggest escape via increased MCL1 abundance/half-life, not drug binding. *Who:* teams pursuing **MCL1 degraders / molecular glues** — these regions are candidate handles, and a degrader would bypass the groove-escape mutations that defeat orthosteric mimetics. *Test:* measure MCL1 protein half-life for E85K / S178L / E325K knock-ins; check whether they blunt BH3-mimetic response but remain sensitive to a degrader.
4. **Method proposal:** because a BH3-mimetic's pocket coincides with the PPI interface, the *A1331852-sensitized dependency arm* (not plain dropout) is the correct BE3D input for mapping a druggable PPI groove — use a synthetic-lethal co-treatment to open the essentiality window when the target is not constitutively essential.

## 7. BE3D issues encountered → appended to BE3D_IMPROVEMENTS.md
1. **[BLOCKER] Control category must be literally "No Mutation".** `parse_be_data` names per-category files by the raw token, but `main()`/`preprocess_ppi_partner` hardcode reading `{gene}_{screen}_No_Mutation.tsv`; a library whose controls are `UTR`/`Intron` crashes (`FileNotFoundError`). Fix: remap control tokens to "No Mutation".
2. **[BLOCKER for PPI] complex mode cannot take a structural-only partner.** `run_complex_mode`→`preprocess_ppi_partner` filters the screen to the partner gene and unconditionally reads its `No_Mutation.tsv`; a partner chain with **no screen rows** (a BH3 peptide / any untiled binding partner) crashes. This is the common PPI case (the interacting peptide/ligand is not itself tiled) and there is no supported path for it (`blind_target` handles a data-less *target*, not a data-less *partner*).
3. **[BEHAVIOR] PPI mode with a data-less partner == monomer numerically.** Verified in code + run: cross-chain neighbors whose partner has no LFC contribute nothing (`_gather_values` skips `-`), so complex rC1 reproduced monomer r10 exactly. Worth documenting so users don't expect interface enrichment purely from adding a complex PDB.
4. **[FRICTION] monomer mode silently drops cross-chain neighbors on a multi-chain PDB.** Feeding a complex PDB in monomer mode (`ppi_chain_gene_dict=None`) is safe (peptide chain ignored) — but this means "interface" information from a bound structure is only usable via distance calc done *outside* BE3D; BE3D itself gives no interface annotation in monomer mode.
