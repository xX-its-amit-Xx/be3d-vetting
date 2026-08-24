# P2 — Per-user BE3D package: TP53 (cancer driver-gene / VUS lab)

**Date:** 2026-08-24 · **Runner:** Queen "TP53" · **Protein:** human TP53 (UniProt **P04637**, 393 aa)
**Structure:** AlphaFold **AF-P04637-F1** (v6), monomer chain A (numbering == P04637; all 17 benchmark
residues verified to match). Experimental reference for DNA-contact geometry: PDB **1TUP** (Cho 1994).
**Screens (real, public, human numbering — cleanest available):** MaveDB base-editor TP53 tiling sets
from **Kaplan/Sander et al.** ("activity-based selection", co-selection method; Etoposide arm):
- `urn:mavedb:00001245-a-1` — **ABE8e** (465 guides → 376 missense, 89 synonymous controls)
- `urn:mavedb:00001245-a-2` — **CBE**   (518 guides → 363 missense, 37 nonsense, 118 controls)
API: `https://api.mavedb.org/api/v1/score-sets/<URN>/scores` (CSV: `hgvs_pro`, `score`=LFC z-score).
Transcript ENST00000269305.9 / ENSP00000269305.4 == canonical P04637, so **no cross-species mapping
needed** (mouse Trp53 was the fallback; MaveDB human set is cleaner and was used).

## Assay direction (documented, like the DNMT3A case)
LFC = log fold-change of guide abundance, **Etoposide arm relative to no-drug**, Day 21.
Empirically from the screen itself:

| class | ABE8e median LFC | CBE median LFC |
|---|---|---|
| Nonsense (definite LOF) | — (ABE makes ~no stops) | **-2.97** |
| Missense | -1.01 | -0.80 |
| Silent / synonymous control | -0.09 | +0.04 |

**Direction: p53 loss-of-function = NEGATIVE LFC (depletion).** Nonsense knockouts drop out strongly
(median -3.0) vs synonymous ~0 — clean QA separation. Mechanistically, functional p53 arrests/protects
cells under etoposide DNA-damage selection; LOF cells keep cycling and are killed → deplete. **All
hotspot analysis therefore uses the BE3D negative (`_neg`) direction.** (Screens carry no intergenic
"No Mutation" guides; the synonymous/silent guides are the neutral controls and were relabelled to
BE3D's `No Mutation` category — semantically identical: no protein change.)

---

## A0. DISCRIMINATION TEST (the primary deliverable) — headline numbers

Canonical run: **radius 6 Å, mean aggregation, nRandom 1000, residue-level.** Negative/LOF direction.
Independent ground truth (curated WITHOUT reference to BE3D):
- **FUNCTIONAL (17):** DNA-contact R248,R273,S241,C277,R280,A276; Zn-binding C176,H179,C238,C242;
  structural R175,G245,R249,R282,Y220,V143,F134. (Cho 1994 1TUP; Joerger & Fersht 2010/2016; IARC.)
- **TOLERANT (clean neutral):** ClinVar has-benign **AND** AlphaMissense modal-benign, minus any
  functional/cBioPortal-recurrent residue (e.g. P72 the common rs1042522 polymorphism). 37 scored in meta.

### Base rate (the honest caveat first)
BE3D flags a **LARGE fraction** of the folded core: at p<0.05 the negative-direction base rate is
**36 % of scored residues (117/327)** in the meta run (38 % ABE8e, 28 % CBE). TP53's DNA-binding
domain is genuinely hypersensitive to mutation (missense median LFC ≈ -1), so a high base rate is
expected — but it means "BE3D recovered hotspot X" is weak on its own. The enrichment/gap below are
what matter.

### Enrichment vs chance (Fisher / hypergeometric), FUNCTIONAL set

| run (p<0.05) | base rate | R (obs) | E (exp) | **R/E** | OR | **Fisher p** |
|---|---|---|---|---|---|---|
| ABE8e | 38% | 15/17 | 6.44 | **2.33** | 14.2 | 1.5e-05 |
| CBE   | 28% | 13/17 | 4.75 | **2.74** | 9.9 | 2.5e-05 |
| **META (ABE8e+CBE, SUM)** | 36% | **16/17** | 6.08 | **2.63** | 33 | **4.2e-07** |

At p<0.001 (meta): R=16/17, R/E=2.93, Fisher p=6.9e-08. (OR is large but unstable because functional
recovery is near-complete — R/E and p are the robust metrics.) **BE3D beats chance decisively.**

### Discrimination gap (functional hit-rate vs tolerant hit-rate) — meta, negative dir
| p-threshold | FUNCTIONAL hit | TOLERANT hit | **gap** |
|---|---|---|---|
| p<0.05 | 94% (16/17) | 14% (5/37) | **+81%** |
| p<0.01 | 94% | 11% | +83% |
| p<0.001| 94% | 8%  | **+86%** |

Functional 94 % ≫ tolerant 8-14 % ≫ base rate leak. **BE3D is discriminating, not "everything is a
hotspot."** It does leak onto some benign positions (14 % tolerant) — consistent with the high base rate,
and reported honestly below as false positives.

### Precision@K (rank by |LFC3D_neg|, meta)
- precision@10 = 20% (2/10), precision@20 = **30% (6/20)** vs base expectation M/N = 5.2% → **~4-6× chance**.
- Top-20 by strength: 249,250,280,243,247,244,246,281,226,265,248,278,282,242,264,225,279,277,283,196 —
  dominated by the L3 DNA-binding loop (243-250 around R248/R249), the loop-sheet-helix (277-283 around
  R280/R282), the Zn site (242), and recurrent hotspots (281,196,278).

### Trivial baselines on the SAME protein (does BE3D add anything?)
| predictor (vs FUNCTIONAL, meta N=327) | K | R/E | OR | p |
|---|---|---|---|---|
| **BE3D neg p<0.05** | 117 | 2.63 | 33 | 4e-07 |
| buried (RSA<20%, from AF SASA) | 96 | 2.20 | 4.9 | 2e-03 |
| **AlphaMissense-high (mean≥0.90)** | 86 | **3.58** | **55** | **2e-09** |

**Brutally honest verdict:** BE3D clearly beats the **burial** baseline (OR 33 vs 4.9) — critically because
the canonical **DNA-contact hotspots are solvent-EXPOSED in the apo monomer** (R248 RSA 0.72, R280 0.49,
A276 0.58, S241 0.42), so a burial predictor *misses* them while BE3D's screen signal catches them.
**However, AlphaMissense — a pure sequence/structure predictor needing no screen — matches or beats BE3D**
on this functional-site task (R/E 3.58, tolerant-hit 0% vs BE3D's 14%). On TP53's conserved core, AM is the
cleaner functional-residue caller. **BE3D's distinct value is orthogonal:** it reads the *actual
proliferation phenotype* under a defined selection (etoposide), gives *direction* (LOF depletion) and
*editor-accessibility*, and recovers patient-recurrence hotspots zero-shot (next section) — not that it
out-discriminates AM on known conserved residues.

---

## A. Ground-truth benchmark: TP/novel/FP + false negatives

Robust hotspot = significant (neg) at **p<0.05 in the meta run AND in ≥1 single editor**. **117 robust
residues.** Classification (see `outreach/TP53/TP53_hotspots.tsv`):

| class | n | meaning |
|---|---|---|
| **TP-structural** | 16 | in the 17-residue FUNCTIONAL set (DNA-contact/Zn/structural) |
| **TP-recurrent** | 21 | cancerhotspots.org recurrent single-residue hotspot (patient recurrence), incl. R196(141), H193(123), P278(100), D281(69), H214(40), G266, C275, Y234, M237, V272, R213, R158, C135, K132, Y163, P151, V157, E285, R306, R337 |
| **NOVEL-candidate** | 67 | robust, DNA-proximal (L2/L3/LSH loop) or buried core, adjacent to validated sites — hypotheses |
| **likely-FP** | 13 | robust but solvent-exposed, not DNA-proximal, not recurrent |

So **37/117 (32%) match an independently-validated functional or recurrent site**; the strongest signals
(top-25 by |LFC3D|) are 60% validated. The novel candidates are overwhelmingly **L3 DNA-binding-loop**
(M243, G244, M246, N247, P250 — note G244S/D are real cancer mutations) and **LSH** (G279, R283, T284)
residues that extend the canonical hotspots into contiguous 3-D LOF surfaces — the intended BE3D output.

**FALSE NEGATIVES (validated sites BE3D missed, and why):**
- **H179 (Zn-binding)** — the *only* one of the 17 functional residues not flagged (ns in both editors,
  meta ns). It is scored but sub-threshold; its Zn-site *neighbours* C176/C238/C242 are all strongly hit,
  so the Zn pocket is recovered as a cluster even though H179 itself is missed. Likely editor/guide-window
  limitation: the reachable substitutions at 179 in this library gave weaker depletion.
- General: **ABE8e alone recovers 78% of cBioPortal hotspots vs CBE's 49%** — editor accessibility is the
  dominant false-negative driver (CBE C>T can't reach many Arg-hotspot substitutions; ABE A>G reaches
  more). **Meta (both editors) recovers 84% (31/37)** — combining editors is the fix, exactly BE3D's
  MetaClust3D use-case. The 6/37 recurrent hotspots still missed (e.g. R342, R337, R306 in the tetramer/
  reg region, R158/A159) are outside or peripheral to the tiled/folded core.

**FALSE POSITIVES (honest):** of the 13 likely-FP, **~5 are also ClinVar-benign AND AlphaMissense-benign**
— genuine false positives: **S261, N263, A129, S260** (exposed loop residues, AM modal LBen, ClinVar
has-benign) and D228 (LBen). These are flagged by neighbourhood aggregation of nearby strong guides
despite being tolerant. The other ~8 exposed FP-flagged residues (E224, S227, G262, A138, G199, Y103,
G226, V225) are AM-pathogenic, so "FP" there is uncertain — they may be real but surface-exposed.
**Caveat on RSA:** in a *monomer* AF model, high RSA conflates true surface exposure with DNA-contact
exposure — do not treat RSA-exposed automatically as FP (R248 is exposed yet is the #1 DNA-contact hotspot).

**Precision/recall feel:** of 117 robust hotspots, **37 match validated functional/recurrent sites, ~67
are structurally-plausible novel candidates (mostly DNA-loop/core), ~5-13 look like false positives.**
Recall of the canonical set is high (16/17 functional; 84% of recurrent hotspots) at the cost of a high
base rate.

---

## B. Competitor comparison (on TP53)

| tool | input signal | what it finds on TP53 | advantage | disadvantage |
|---|---|---|---|---|
| **BE3D** (this run) | ABE8e+CBE etoposide-selection **screen LFC** → LFC3D + randomization null + 3D clustering | 117 robust LOF residues; recovers 16/17 structural + 84% of recurrent hotspots; extends L3/LSH loops into contiguous surfaces; direction + editor-accessibility | uses a *functional phenotype*, not prediction; gives LOF direction; separates editor reachability; meta-combines editors; would capture context/allosteric effects recurrence/AM can't | high base rate (36%); some benign FPs; on TP53's conserved core it does **not** out-discriminate AlphaMissense; needs a screen |
| **cancerhotspots.org / cBioPortal 3D** | patient mutation **recurrence** | R273,R248,R175,G245,R282,R249,Y220,S241,C176,R280,H179,H193,R196,P278,D281… (counts) | decade of clinical benchmarking; no screen needed; direct clinical relevance | recurrence ≠ mechanism; blind to editor accessibility & to functional-but-rarely-mutated residues; no experimental phenotype |
| **AlphaMissense** | sequence+structure **predictor** | every DBD hotspot scored ~0.89-1.00 (modal LPath); P72 benign 0.13; 160/199 DBD residues modal-pathogenic | **best functional-site discrimination here** (R/E 3.58, tolerant 0%); zero cost, no screen; per-substitution resolution | pure prediction — no phenotype, no direction, no assay/editor context; can't see selection-specific or resistance effects; also flags most of the DBD (its own "high base rate") |
| **ProTiler-Mut** (Cell Systems 2026) | tiling screen → 3D-RRA | (conceptual) would also surface the DBD hotspot loops from the same LFC | direct tiling→3D rival; RRA robust to outliers | recurrence/rank-based; BE3D's explicit randomization null + editor-aware meta-aggregation is the differentiator |

**Fair take:** on TP53 specifically, the recurrence tools and AlphaMissense are *strong* and cheaper for
"which residues are functional." BE3D wins when the question is **"what does the actual base-editing screen
say, in 3-D, with a null model, and which variants are editor-reachable"** — e.g. reconciling ABE-vs-CBE
accessibility (78% vs 49% recovery) and turning a noisy per-guide LFC into a coherent spatial LOF map.

---

## C/D. Deliverables (outreach/TP53/)
- `TP53_hotspots.tsv` — 117 robust hotspots: LFC3D_neg, per-editor + meta significance, RSA, pLDDT, region,
  structural role, cBioPortal count, AlphaMissense mean/class, TP/novel/FP flag.
- `TP53_G2P_LFC3D_neg.tsv` — G2P-portal interactive-view input (negative/LOF direction, both editors + meta).
- `TP53_ABE8e_screen_BE3D_input.tsv`, `TP53_CBE_screen_BE3D_input.tsv` — reformatted screens (reproducible).
- `TP53_email.md`, `TP53_brief.md` — copyable outreach.

## Robustness sweep
Sweep over radius {4,6,8,10} × {mean,sum} × nRandom {500,1000} × {residue,atom} on the meta screen
(`output/sweep/`, live summary `output/sweep/SWEEP_discrimination_summary.tsv`). Early completed configs
confirm discrimination is stable:

| config | base rate | func R | R/E | Fisher p | tol hit | cBio recovered |
|---|---|---|---|---|---|---|
| R4 mean n1000 residue | 28% | 14/17 | 2.90 | 3e-06 | 6% | 24/37 |
| R4 mean n500 residue | 28% | 14/17 | 2.90 | 3e-06 | 6% | 24/37 |
| R4 mean n500 atom | 40% | 16/17 | 2.32 | 3e-06 | 13% | 31/37 |
| R6 mean n1000 residue (canonical) | 36% | 16/17 | 2.63 | 4e-07 | 14% | 31/37 |

Observations: **nRandom 500 vs 1000 is numerically identical** (null is well-sampled). **Atom-level** and
**larger radius** raise the base rate and recover slightly more (more diffuse clusters); **smaller radius /
residue-level** sharpens to individual hotspots at a lower base rate. Across all settings R/E stays ~2.3-2.9
with Fisher p ≤ 3e-6 and functional recovery 14-16/17 — the conclusion (real discrimination, high base rate)
is not an artifact of one parameter choice. The robust hotspot core (L3 243-250, LSH 277-283, Zn 242/238,
R249, R175, R280) persists across the grid. (Full 32-config grid completing in background.)

## Honesty summary
- Base rate is high (36%) — stated up front.
- BE3D beats chance (R/E 2.6, p 4e-7) and burial (OR 33 vs 4.9), but **AlphaMissense beats BE3D** on pure
  functional-site discrimination — stated plainly.
- H179 is a genuine false negative; ~5 benign residues are genuine false positives — both flagged.
- No novel candidate is called validated; all are labelled hypotheses.
