# BE3D Run — Decision Log — MYC (HONEST HARD discrimination test: folded bHLH-LZ vs disordered bulk)

Runner: Queen "BCL2/MYC" (round 2). Date: 2026-08-24. venv py3.14 / beclust3d clone.
Outputs: `real_output/P2_MYC/runs/{AF_r4_mean,AF_r6_mean,AF_r8_mean,AF_r6_sum}/MYC/`,
`real_output/P2_MYC/discrim_bm_report.json`, `outreach/MYC/`.

## 0. Target, data reality, and the honest question
- Gene / protein: **MYC** (UniProt **P01106**). **Numbering hazard (critical):** the AlphaFold model
  `AF-P01106-F1-model_v6` and the Coelho screen are the **454-aa (CUG-initiated MYC1 / SV=2)** frame;
  the classic literature / 1NKP / AlphaMissense are the **439-aa** frame. **canonical-454 = 439-number + 15.**
  The screen's `Amino_Acid_Position` matches the AF sequence at **offset 0** (328/328 single-pos missense
  refAA match) — so all BE3D I/O is in **454-numbering**; every literature residue below is stated in 454.
- Structure: **AF-P01106** (primary). Folded blocks (pLDDT>70): only **203-209** (blip) and **365-453**
  (the **bHLH-LZ**). Everything else, **1-364, is a low-pLDDT intrinsically disordered region (IDR)** —
  exactly the brief's expectation (23% of residues pLDDT>70; mean pLDDT 60). Corroborating complex: **1NKP**
  (MYC-MAX-DNA); MYC=chain D, auth→454 offset = authD−132.
- **Data reality (important):** the Coelho ST2 sheet has **no MYC-targeted drug**; the MYC-relevant readout
  is **essentiality** via the **Control (no-drug) dropout arms**. MYC is essential in these lines, so this is
  a clean, powerful readout (unlike a resistance arm). Chosen arm: **HT29 Control** (strongest KO window:
  nonsense mean z = **−2.67**; H23 −1.60; PC9 −1.39). Interpreted direction = **NEG (depletion = essential)**.
- **The honest HARD test (like BRCA1):** does BE3D **concentrate** essentiality signal on the small folded
  **bHLH-LZ** (DNA-binding basic region + MAX-dimerization + LZ) and **avoid the disordered bulk**, or does it
  **over-call across the IDR**?

## 1. Data provenance & reformatting
- Source: Coelho, Dincer et al. 2024 *Nat Genet* (10.1038/s41588-024-01948-8), MOESM4 `ST2 BE z-scores`,
  local `scratchpad/coelho_S4.xlsx`. **2,498 MYC guide rows** (pooled ABE+CBE).
- Converter `scratchpad/convert_bm.py` (same schema as `convert_kras.py`): category from
  `most_severe_consequence`; edits from `zip(Amino_Acid_Position, Edited_AA, New_AA)`; controls→`nc`.
  Output `scratchpad/bm_tsv/MYC_HT29ctrl.tsv`: **2,478 rows** = 801 Missense (337 distinct positions,
  span 1-454), 157 Silent, 34 Nonsense, 25 Splice-donor, 6 Intron, 1,455 No-Mutation.
- QA (cases=Nonsense+Splice-donor, controls=No Mutation+Silent): **KS D=0.354, p=7.1e-7** (59 vs 1,612) —
  a real, strong essentiality window. **ACCEPT.** (Coelho drug arms would be null here; the Control arm is
  the correct essentiality gate.)

## 2. Parameter sweep (4 runs, all EXITCODE 0)
`structure_radius ∈ {4,6,8}` mean + `r6 sum`; clustering_radius 6; nRandom 1000; AF monomer.
Driver `scratchpad/run_bm_sweep.sh`. p<0.05/0.01/0.001 computed every run.

## 3. A0 DISCRIMINATION BENCHMARK (the primary result) — `discrim_bm_report.json`
Ground-truth sets are **BE3D-independent** (1NKP contacts + UniProt domains):
- **FUNCTIONAL** = 1NKP DNA-contacts {370,371,374,375,378,379,381,382,386,406,407} ∪ MAX-dimerization
  {384,385,388,389,391,392,395,408,411,412,414,415,417-419,421,422,424,425,428,429,431,432,435,436,438,439,442,443,445-449} (M=45).
- **TOLERANT/IDR bulk** = `range(79,143) ∪ range(159,369)` (274 res), i.e. disordered, **excluding** the
  functional-but-disordered conserved boxes **MBI (60-78)** and **MBII (143-158)**.
- Baselines: burial (AF CN top-tertile, 191 res), **folded (pLDDT>70, 105 res)**, AlphaMissense>0.7 (143 res).

**Baseline (AF, r6, mean, NEG direction), N_neg universe = 255:**

| metric | value |
|---|---|
| base rate p<0.05 | **32.5%** (83/255 sig) — high-ish, but MYC is genuinely essential |
| **FUNCTIONAL enrichment** p<0.05 | **R/E = 2.11, OR = 5.84, p = 7.7e-6** |
| FUNCTIONAL enrichment p<0.01 | R/E = 2.31, **OR = 6.89, p = 1.2e-6** (r8: OR 6.8, **p = 2.6e-7**) |
| **discrimination gap** (hit_FUNC/hit_IDR) | **4.31×** (68.8% vs 15.9%); r4 4.29×, r8 4.2× |
| domain hit-rate | **bHLH-LZ 64.8% vs IDR 15.9%** (r8: 75% vs 17%) |
| precision@10 (FUNC) | **0.6** (r6), **0.7** (r8); top-10 = R379/R381/R382 (DNA), Y417/I418/L419/V421 (MAX-LZ), C148/M149 (MBII) |

**Sweep-stable:** enrichment OR 3.7–9.9, p ≤ 3e-3 across all four runs; gap 3.1–4.3×. Radius 4→8 broadens
the neighborhood but the FUNC concentration and the folded-vs-IDR gap **hold at every radius** (r6-sum is the
weakest, OR 3.1–3.7, but still significant).

**Beats-baseline (the decisive nuance):**
- vs IDR bulk: **YES, decisively** — BE3D flags the bHLH-LZ at ~4× the IDR rate; it does **not** diffusely
  over-call the disordered region. This is a genuine PASS of the folded-vs-disordered stress test.
- vs "is it folded": at the **whole-protein** level, "call all pLDDT>70" gives FUNC R/E=**3.45** — a *stronger*
  raw predictor than BE3D (R/E~2), because every FUNC residue lives in the one folded domain. **Within the
  folded bHLH-LZ, BE3D does NOT further resolve DNA/MAX contacts from non-contact domain residues**
  (within-domain FUNC R/E ≈ **1.0**, OR 1.3–1.5, p 0.33–0.78; contact hit 69–73% vs non-contact-in-domain
  59–78%). So BE3D's discrimination is **DOMAIN-level, not residue-level**.
- **What BE3D adds beyond folding (the honest value-add):** it recovers the **functional IDR boxes a pure
  folding/pLDDT predictor cannot see** — **all of MBII (145-158)** and part of **MBI (60-65)** are robustly
  significant (low pLDDT 51-59, so invisible to a burial/folded baseline). These are true functional sites
  (MBII's DCMWSGF core recruits TRRAP/HATs; C148/M149 are among the strongest depletions, LFC3D −2.9/−3.5).

## 4. Robust hotspots (≥3/4 runs sig p<0.05) — `outreach/MYC/MYC_hotspots.tsv`
**74 robust hotspots.** Annotated: **19 validated 1NKP contacts** (6 DNA: R379,R381,R382,K386,K370-region;
13 MAX: Y417,I418,L419,V421,L428,L435,L442… the leucine-zipper), **12 other bHLH-LZ-domain**, **16 MBI/MBII**
(functional-IDR), **1 novel-folded**, **26 likely-FP disordered-IDR**. → **48/74 (65%) functional/plausible;
26 residual IDR over-calls.** Strongest depletions (LFC3D_neg): **Y417 −4.49, A416 −4.39, L419 −4.07,
I418 −4.02, M149(MBII) −3.52, R379(DNA) −3.47, S420 −3.15, Q380 −2.94, C148(MBII) −2.90, V421 −2.71,
R382/R381 −2.67/−2.65** — i.e. the strongest signals sit **exactly on the DNA-binding basic arginines, the
MAX leucine-zipper, and MYC Box II**.

## 5. Interpretation & HONEST verdict
- **BE3D PASSES the discrimination stress test at the domain level.** From screen signal alone it concentrates
  MYC essentiality on the folded bHLH-LZ (OR up to ~10, p to 3e-7; 4× over the IDR) and correctly lights up the
  functional conserved boxes MBI/MBII inside the IDR — it does **not** over-call the disordered bulk uniformly.
  This is the opposite of a BRCA1-style "everything is a hotspot" failure.
- **Honest limits:** (1) BE3D's resolution is **domain-level, not residue-level** — within the bHLH-LZ it does
  not separate the actual DNA/MAX contacts from neighbouring domain residues (a 6 Å neighborhood on a compact
  helical domain smears signal across the whole fold). (2) It still emits **~26 scattered IDR false positives**
  (base rate ~32%), clustered around res 37-47, 92-95, 120-127, 317-319 — treat p<0.05 as a screen, not a
  verdict; use magnitude/rank + the enrichment test. (3) The screen cannot say *which* function (DNA vs MAX vs
  cofactor) each residue serves — that annotation comes from 1NKP, not BE3D.
- **Value BE3D genuinely adds over trivial predictors:** it flags the **low-pLDDT MBI/MBII functional boxes**
  that a folding/burial/pLDDT baseline misses, and it does so **from a functional dropout screen** (no prior
  structure knowledge of which region is engaged). **AlphaMissense** (R/E 2.09, OR 5.71, p=1e-5, after the
  +15 frame correction) is an **equally strong** whole-protein concentrator of the contact set — BE3D does not
  beat it on the folded contacts; BE3D's distinctive add is the **low-pLDDT MBI/MBII recovery** that AM under-weights.

## 6. False positives / false negatives
- **FALSE POSITIVES:** 26 IDR residues (37-47, 92-95, 120-127, 317-319) — low-pLDDT, no known function,
  likely bystander/aggregated-guide noise amplified by the neighborhood sum. Flagged in the hotspot TSV.
- **FALSE NEGATIVES / non-resolution:** individual DNA-contact specificity within the domain is lost
  (within-domain enrichment ≈ chance). No residues are truly "missed" (the whole functional domain is hit),
  but the method cannot rank the DNA basic arginines above their folded neighbours.

## 7. Competitor note
- **AlphaMissense** (orthogonal, needs no screen) concentrates on the contact set as well as BE3D (R/E 2.09,
  OR 5.71) but, like a folding
  predictor, under-weights the low-pLDDT MBI/MBII boxes — BE3D's screen-driven call recovers them.
- **3D somatic-recurrence tools** (3dhotspots.org / HotMAPS) find little on MYC (MYC is amplified/translocated,
  rarely point-mutated), so they would be near-silent here — BE3D's functional-screen readout sees essentiality
  the recurrence tools cannot.
- **ProTiler-Mut** (tiling→3D-RRA) is the closest rival; conceptually it would also localize to the folded
  domain but shares the same domain-vs-residue resolution ceiling on a compact fold.

## 8. BE3D issues → BE3D_IMPROVEMENTS.md
(1) 454-vs-439 MYC numbering trap; (2) high base rate on essential genes needs an effect-size/FDR gate;
(3) domain-level-only resolution on compact folds; (4) no MYC-targeted arm in Coelho — essentiality (Control
arm) is the honest MYC readout. See appended entries.
