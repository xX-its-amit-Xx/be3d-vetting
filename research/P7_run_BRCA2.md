# BE3D Run — BRCA2 (huge, mostly-disordered tumor suppressor; where BE3D STRUGGLES)

Runner: Queen "BRCA2" (P7). Date: 2026-09-03. Patient-impact framing: BRCA2 = hereditary breast/ovarian
cancer; thousands of missense VUS. Question probed: **on a 3418-aa, ~4/5-disordered protein, does BE3D
concentrate signal on the folded functional core (DNA-binding domain / BRC repeats / tower) or OVER-CALL
the disordered bulk — and can it even be run at all?**
Outputs: `real_output/P7_BRCA2/`, `outreach/BRCA2_HANNA/`. Clone: patched fork
`xX-its-amit-Xx/BE3D@fix/portability-and-dssp-fallback` at `C:/Temp/p7_brca2`; runs under `C:/Temp/p7out/`.

## 0. Target, data, structure
- **Gene/protein:** BRCA2, UniProt **P51587, 3418 aa** (canonical). Architecture: a long N-terminal +
  central **intrinsically-disordered / low-complexity region (res ~1–2478)** carrying the 8 **BRC repeats**
  (RAD51-binding motifs: 1002–1036, 1212–1246, 1421–1455, 1517–1551, 1664–1698, 1837–1871, 1971–2005,
  2051–2085; UniProt), then the single globular module — the **DNA-binding domain (DBD, ~2470–3190)**:
  helical domain (2479–2667, binds DSS1/SEM1), OB1 (2670–2803), OB2 (2809–3048) with the **tower domain**
  (2831–2872, binds dsDNA), OB3 (3052–3190) — then a disordered C-terminal tail incl. TR2/RAD51 site
  (3265–3330). **Only the DBD (~721 aa = 20.8% of the protein) is a large folded module**; the rest is
  disordered/low-complexity (MobiDB flags disorder; structurally, ~79% has no solved or predictable fold).
- **Screen:** Hanna et al. 2021 *Cell* "Massively parallel assessment of human variants with base editor
  screens" (PMID 33461210), Table **S2** (`mmc2.xlsx`, curl 200). CBE tiling of BRCA1/2; **BRCA2 sheet =
  raw counts** (pDNA + HAP1 Dropout ×3, MELJUSO Dropout/Cisplatin/Talazoparib ×3) — **not** a ready
  Z-score (the ready-Z "BRCA1 & BRCA2 missense" sheet has only 13 validation rows), so **LFC computed from
  counts**. 769 guides join 1:1 to the BRCA2 Library annotation sheet (`Amino acid edits`, `Mutation
  category`, `Clinical significance`=ClinVar). refAA match vs P51587 = **348/348 (100%)** → screen
  numbering == canonical.
- **Assay direction (documented):** dropout — **LOF = negative LFC** (no negation; LOF lands in BE3D's
  `neg` channel). Verified internal controls (HAP1): nonsense mean **−1.46** vs No-Mutation **−0.01**,
  silent −0.18, missense −0.30. Primary arm chosen by LOF window strength (see §1).
- **Structure (THE FIRST STRUGGLE — see §4.1):** AlphaFold DB **does not host BRCA2** (>2700-aa cap), so
  BE3D's auto-fetch is impossible; the PDB has only tiny BRCA2 fragments (2271–2343, 3260–3308). Used a
  **SWISS-MODEL homology model of the DBD, res 2470–3190** (template 1IYJ/1MIU, human P51587 numbering,
  GMQE 0.051 — low; via 3D-Beacons). DSSP = pure-Python placeholder (characterization-only).

## 1. Data reformatting + arm choice (decision log)
Script `convert_brca2.py` → BE3D 4-col TSVs (`brca2_<arm>.tsv`). CPM-normalize each count column, then
`LFC = log2((mean(CPM_arm)+0.5)/(CPM_pDNA+0.5))`. 3-letter→1-letter AA edits (`Ile21Ile;`→`I21I;`),
multi-edit `;`-joined; non-editing guides (462, blank category) → **No Mutation** controls; Intron kept.
Category counts (per-guide, pre-collapse): 197 Missense / 51 Silent / 33 Nonsense / 26 Intron / 6 Splice /
462 No-Mutation. AA-edit positions span **2–3418** (whole protein, sparse: 348 missense edits — 219 pre-DBD
incl. 86 in BRC region, 84 in DBD, 45 C-terminal).

**Arm QA (nonsense vs No-Mutation, KS/MWU), 4 arms computed:**

| arm | nonsense mean | No-Mut mean | KS D | KS p |
|---|---|---|---|---|
| **HAP1 dropout (PRIMARY)** | **−1.46** | −0.01 | **0.77** | **2.5e-19** |
| MELJUSO dropout | −0.60 | −0.03 | 0.50 | 1.6e-7 |
| MELJUSO cisplatin | −0.71 | +0.02 | 0.61 | 2.4e-11 |
| MELJUSO talazoparib (2nd) | −0.69 | +0.01 | 0.53 | 1.3e-8 |

**HAP1 dropout** (haploid line — BRCA2 directly essential) has by far the cleanest LOF window →
**PRIMARY**. **Talazoparib** (PARPi; synthetic-lethal with BRCA2-LOF, like BRCA1's olaparib arm) run as a
**robustness** second arm. QA PASSES decisively on both.

## 2. A0 DISCRIMINATION BENCHMARK — the primary result (lead with this)
Run: monomer, SWISS-MODEL DBD structure, radius 6, nRandom 500, mean, neg (LOF) channel. Because the
structure covers only the DBD, BE3D emits `LFC3D_neg` for **233 residues** — 134 inside the DBD (real 3D
neighborhoods) and **99 OUTSIDE the structure** (see §2.1 — a bug, not a feature).

| Metric (HAP1 dropout, p<0.05 neg) | Value | Honest reading |
|---|---|---|
| **Base rate overall** (sig/scored, 233) | **38.6%** | WAY above the 15–20% "everything-is-a-hotspot" bar |
| **Base rate in DBD** (folded, 134) | **41.8%** (56) | compact-fold over-call (as BCL2/AKT1/MEK/MYC) |
| **Base rate in IDR / out-of-structure** (99) | **34.3%** (34) | flags the disordered bulk almost as often |
| **Discrimination gap DBD − IDR** | **only +7.4 pts** | NOT concentrating on the folded core |
| **Functional enrichment** (DBD ClinVar-P missense; M=2 scored of 9) | R/E=**1.2**, OR=**1.4**, **p=0.66** | **n.s.** — no residue-level discrimination |
| **Disc. gap functional vs tolerant** (DBD) | +0.07 (0.50 vs 0.43) | essentially none (bystander smearing) |
| **Precision@10 / @20** (ClinVar-P) | **1/10, 1/20** | top ranks = DBD folding-core, ClinVar-unvalidated |
| **Beats "is-it-folded (pLDDT>70)" baseline?** | **NO** | same verdict as BRCA1: domain signal only, no fine resolution — and here pLDDT doesn't even exist (§4.1) |

**Talazoparib arm (robustness):** DBD base rate **33.6%**, **IDR base rate 37.0% (> DBD!)**; 39 DBD
hotspots shared with HAP (Jaccard **0.60**). The over-call of the disordered bulk is reproduced — in this
arm the IDR is flagged *more* than the folded DBD.

**Verdict:** On BRCA2, BE3D **fails the discrimination test**. It over-calls (base rate 34–42%), and the
IDR-vs-DBD gap is ~0 (+7 pts HAP, −3 pts TALA) — it does **not** concentrate signal on the folded
functional core, and its only real signal (the DBD folding core) is exactly what a trivial burial/pLDDT
predictor already gives. This is a stronger negative than BRCA1 (which at least had a clean +0.13 domain gap).

### 2.1 The self-leak that drives the over-call (verified mechanism)
For all 99 out-of-structure residues, **`LFC3D_neg` == the residue's own 1-D `LFC_neg` exactly**
(max|Δ| = 0.0). BE3D's `calculate_lfc3d` includes the focal residue's own value and **does not require
≥1 structural neighbor**, so any edited residue with no coordinates (the entire disordered bulk) is emitted
as a "3D hotspot" indistinguishable in the output from a genuine spatially-aggregated call. Result: **34
(HAP)/37 (TALA) "significant 3D hotspots" scattered across res 69–2395 and 3290–3405 with zero spatial
support** — pure single-guide LFC. Restricting the structure to the folded DBD (the safe move) does **not**
protect against this. This is the single sharpest new finding of the run.

## 3. Hotspot calls (FP/FN flagging) — `outreach/BRCA2_HANNA/BRCA2_HANNA_hotspots.tsv` (113 residues)
- **TP / genuine recovery (coarse):** the **DBD folding core** is the strongest, robust (2-arm) cluster —
  **D2723** (ClinVar-**Pathogenic**, a classic folding-disruptive missense) plus the buried hydrophobic
  OB-fold core **2582/2588/2595, 2657, 2686–2725 (incl. W2725), 2792/2793, 2972–2999, 3054–3080,
  3092–3125**. These are structurally coherent mutation-intolerant positions (37 robust in both arms). So
  BE3D **does** localize the DBD's intolerant core — but only at domain-block resolution, and mostly
  ClinVar-unvalidated because base-editing coverage of the *validated* pathogenic residues is sparse (FN).
- **LIKELY FALSE POSITIVES:** (a) **48 self-leak IDR "hotspots"** (§2.1) — flagged explicitly, no
  structure, single-guide LFC (e.g. 69, 297, 561, 942, 1263, 1687, 1932, 3290, 3389). (b) **DBD-benign
  3079** (ClinVar Benign) flagged in both arms = aggregation bystander (shares its neighborhood's value,
  like BRCA1's L63/K65). (c) The identical clustered LFC3D values (−2.795, −2.695, −2.585 across whole
  neighborhoods) show the hard-radius smearing that inflates the DBD base rate.
- **FALSE NEGATIVES:** (a) **Reachability** — only **2/9 DBD ClinVar-pathogenic residues** are covered by a
  base-editable missense guide; the other 7 (and most literature DBD pathogens) are unreachable by CBE →
  missed by construction, not by the algorithm. (b) **The BRC repeats / RAD51 interface (1002–2085) and the
  entire disordered functional bulk have NO structure** → cannot be spatially analyzed at all; TR2 (3265–
  3330) likewise. A huge fraction of BRCA2's known biology is structurally invisible to BE3D.

## 4. WHERE BE3D STRUGGLED (explicit, with numbers)
### 4.1 It cannot be run out-of-the-box at all (the headline)
- **AlphaFold DB does not host BRCA2** (P51587): `api/prediction/P51587` → `{}`, and `AF-P51587-F1-model_v{1..6}.pdb`
  all **404**, whereas **BRCA1 P38398 → 200** (worked in the earlier run). AFDB excludes proteins >2700 aa
  (Titin Q8WZ42 also `{}`). BE3D's flagship workflow — *"give a UniProt ID, we auto-fetch AlphaFold"* —
  **dies at the structure stage with a 404** for any large protein. No full-length experimental structure
  exists either (PDB: only 2271–2343 & 3260–3308 fragments). I had to hand-source a SWISS-MODEL DBD model
  via 3D-Beacons. **For the ~2/3 of the disordered proteome that is large and/or unmodeled, BE3D has no
  structure to stand on.**
- **No pLDDT exists** for BRCA2 (no AF model) → any "pLDDT gate / confidence shortlist" (PR #23) has no
  input; a fallback confidence signal (disorder predictor, model QMEAN, or a "no-structure" flag) is required.

### 4.2 Over-calling with no spatial support (§2, §2.1)
Base rate 34–42%; discrimination gap ~0; **48 self-leak IDR FPs**; functional enrichment n.s. (p 0.66).

### 4.3 Scaling of the neighbor search (O(N²), pure Python)
`count_aa_within_radius` is a **naive double `for` loop over all residues with per-pair `math.sqrt`, no
KD-tree, no vectorization**. On the 721-residue DBD the full run is **68 s (HAP) / 98 s (TALA)**, sub-GB
RAM — fine. But it is **O(N²)**: full-length residue-level (N=3418) ≈ **22× the pairwise ops** (~11.7M);
**atom-level** (`atom_level_naa=True`, ~27,000 atoms) ≈ **7×10⁸ pairwise sqrt in pure Python** — the real
blow-up (minutes→hours + large neighbor-list RAM). A `scipy.spatial.cKDTree` makes this O(N log N). So even
if a full BRCA2 model existed, atom-level analysis of it would be impractical as written.

## 5. PROPOSED CHANGES (cited by the analysis above; PR status noted)
1. **[NEW, highest-value] Require ≥1 structural neighbor before emitting an LFC3D call** (or split
   "spatial LFC3D" from "self-only 1-D LFC" in the output + significance). Directly kills the **48 self-leak
   IDR false positives** (§2.1). **Not covered by any open PR.** This is the cleanest, most surgical fix and
   the one most specific to large/partially-modeled proteins.
2. **[NEW] Graceful handling when AlphaFold is absent (protein >2700 aa / not in AFDB).** Detect the 404 /
   empty API and (a) fail with an actionable message, (b) auto-try 3D-Beacons/SWISS-MODEL model providers,
   and (c) fall back to a disorder-predictor confidence track when no pLDDT exists (§4.1). Feeds PR #23's
   shortlist, which currently assumes pLDDT is available. **Not covered.**
3. **[NEW] Scalable neighbor search (KD-tree).** Replace the O(N²) double loop with `cKDTree.query_ball_point`
   → O(N log N); essential for large multidomain proteins and any atom-level run (§4.3). **Not covered.**
4. **[partly covered] Base-rate + FDR reporting and a folded-only / domain-restricted analysis mode.** Print
   `#sig/#scored` (and per-region base rate) in `RUN_COMPLETED`, apply an effect-size floor + BH/FDR, and let
   the user restrict scoring/significance to a residue range or a pLDDT/disorder threshold. FDR q-values are
   **PR #20**; the confidence-weighted "what-to-validate-next" shortlist (FDR+pLDDT+reachability, A/B/C tiers)
   is **PR #23**; a GP/pLDDT spatial kernel + cluster-level spatial-FDR are in **issue #24**. BRCA2's 34–42%
   base rate and 48 self-leak FPs are direct motivation.
5. **[covered — PR #19] Base-editing reachability report.** Only **2/9 DBD ClinVar-pathogenic residues** are
   editable here; the reachability report would mark the other 7 as assay-missed, not BE3D failures.

## 6. Competitor context (brief)
| Tool | Signal on BRCA2 | Edge / limit |
|---|---|---|
| **BE3D** | DBD folding-core cluster from HAP1 dropout (D2723 + OB core), robust across arms | functional screen, germline-capable; but over-calls (34–42%), self-leak IDR FPs, needs a structure it can't fetch, coarse (domain-level) resolution |
| **AlphaMissense** | scores all 3418 residues incl. the whole IDR, no screen or structure needed | proteome-wide, folded-core strong; a predictor not a measurement — but it *works* where BE3D can't even fetch a model |
| **BRCA2 HDR / SGE functional assays** (Guidugli, Mesman, Ikegami, Richardson) | per-variant DBD function (the real gold standard) | orthogonal truth; DBD-focused, heavy engineering |
| **3dhotspots / cBioPortal** | ~nothing (BRCA2 is a germline tumor suppressor, non-recurrently somatic) | blind to germline TS → BE3D's niche in principle |

**Where BE3D could win:** germline tumor suppressor + functional dropout readout. **Where it loses here:**
it needs a structure large disordered proteins don't have, over-calls, and adds no residue-level resolution
beyond "is it in the folded DBD" — which AlphaMissense/pLDDT give for free across the *whole* protein.

## 7. Deliverables
- `outreach/BRCA2_HANNA/BRCA2_HANNA_hotspots.tsv` (113 sig residues + TP/novel/self-leak-FP/bystander flags +
  ClinVar + region), `BRCA2_HANNA_G2P.tsv` (NEG G2P view), `BRCA2_HANNA_brief.md`.
- `real_output/P7_BRCA2/`: per-arm `*_NonAggr_LFC3D.tsv`, reformatted `brca2_HAPdrop/MELtala.tsv`, audit TSV.
- Scratch: `scratchpad/p7/` (convert_brca2.py, struct/BRCA2_DBD.pdb, analysis).
```
```
