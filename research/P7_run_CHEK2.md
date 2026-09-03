# BE3D Run — CHEK2 (clean kinase; VUS / ClinVar discrimination stress-test)

> Runner: Queen "ATM/CHEK2". Real public DDR base-editing screen (Cuella-Martin 2021 *Cell*) →
> BE3D LFC3D + 3D clustering → does BE3D discriminate ClinVar/AlphaMissense-pathogenic kinase
> residues from benign, from screen signal alone? A direct **VUS-relevance** test on a compact,
> well-folded FHA+kinase protein. Central deliverable = the **A0 DISCRIMINATION BENCHMARK**.
> Patient-impact framing: CHEK2 is a moderate-penetrance breast/prostate cancer gene swamped with VUS.

## 0. Target, data, direction
- **Gene/protein:** CHEK2 (CHK2), UniProt **O96017**, 543 aa. Domains (UniProt): **SCD/SQ-TQ 19–69**
  (contains **T68**, the ATM phosphosite), **FHA 113–175** (pThr-peptide binding), **kinase 220–486**.
- **Structure:** **AlphaFold AF-O96017-F1-model_v6** (543 aa, whole-protein numbering = unipos; supplied
  as `user_pdb` so its B-factor column carries the real per-residue **pLDDT**). Monomer mode, chain A.
  DSSP via the fork's pure-Python fallback (SS/RSA are placeholders → do NOT use BE3D's RSA/burial
  characterization; burial baseline below is computed from `Naa_count` instead).
- **Screen:** Cuella-Martin et al. 2021 *Cell* (PMID 33592168), DDR CBE (BE3-FNLS) tiling, Table **S5**
  (`cm_mmc5.xlsx`, curl 200), sheet **`MCF10A treatments Subl1`** (non-transformed line; same file/sheet
  as the BRCA1 run). **205 CHEK2 guides.** `AAChg` per-token parsed (recipe = P2_run_BRCA1 §1):
  `pX{n}Y`→Missense, `pX{n}=`→Silent, `pX{n}*`→Nonsense; delins/complex skipped. `ClinVar` retained.
- **Assay direction (documented — DIFFERS from BRCA1):** internal controls show **nonsense guides are
  ENRICHED, not depleted** (OLAP nonsense mean **+1.17** vs synonymous **+0.23**; KS D=0.60 p=1.3e-3;
  CISP/CPT/DOX similar). Biology: CHK2 loss ablates the DNA-damage checkpoint/apoptosis → cells **survive
  and out-grow** under a DNA-damaging drug → **LOF = POSITIVE LFC (a survival/resistance screen)**, the
  **opposite** sign to BRCA1's synthetic-lethal dropout. We therefore **NEGATE** the LFC so LOF lands in
  BE3D's `neg` channel (matching the project-wide "neg = LOF" convention; cf. DNMT3A). **Primary arm =
  olaparib (OLAP)** — the single best control separation for CHEK2 (KS D=0.60); CISP/CPT are comparable.
- **BE-QA (in-run, negated):** cases=[Nonsense] vs controls=[No Mutation]: **KS D=0.49, p=5.5e-3;
  MWU p=…** → screen **PASSED**. Reformat: `p7data/convert_cm.py`; input `CHEK2_OLAP.tsv`
  (99 Missense / 30 Silent / 14 Nonsense / 2+2 Splice / 58 No-Mutation guides).

## 1. Sweep
Monomer AF, r4 / **r6 (PRIMARY)** / r8, nRandom **1000**, mean/mean, meta SUM (single screen). All
EXITCODE 0. Outputs `C:/Temp/p7_out/CHEK2_r{4,6,8}/`. Analysis: `benchmark.py`, `make_hotspots.py`.

## 2. A0 DISCRIMINATION BENCHMARK (the primary result) — neg (LOF) channel, r6
Universe = 543 modeled residues; "scored" = 260 residues a missense guide actually reaches.
Independent gold standard = **AlphaMissense** (full-19 per-residue, authoritative AF-DB table for O96017:
frac mean-AM>0.564 = **50.1%**) + a small **cited catalytic/FHA set** (UniProt O96017 + Li 2002 *Mol Cell*
9:1045 FHA-pThr; Cai 2009; UniProt active/binding sites).

| Test | Number | Honest reading |
|---|---|---|
| **Base rate** (sig / scored) | **38.8%** (101/260) p<.05; 38.1% p<.001 | **HIGH — over-calls.** Of modeled residues 18.6%. The no-FDR gap, loud. |
| base rate vs radius (of scored) | r4 34.5% / r6 38.8% / r8 40.5% | radius-robust over-call |
| **Enrichment** for AlphaMissense-pathogenic (M=272) | R=78, E=50.6, **R/E=1.54, OR=4.34, p=6.7e-10** | **highly significant** — hotspots concentrate on predicted-functional residues |
| **Discrimination gap** (hit-rate AM-high vs AM-low) | **28.7% vs 5.1% (+0.236)** | real: functional ≫ tolerant |
| **Precision@10 / @20** (frac AM-pathogenic) | **1.00 / 1.00** (chance 0.50, **lift 2.0×**) | **top signals are 100% functional** — ATP pocket 303–357, catalytic-loop 350–352, FHA 113–161 |
| domain hit-rate kinase / FHA / SCD | **26.2% / 19.0% / 19.6%** | concentrates on the kinase domain |
| confidence: hit-rate pLDDT<70 vs ≥70 | **10.8% vs 22.3%** | **does NOT over-call the disordered N-term** (folds-first) |
| **Beat burial baseline** (top-K by Naa_count) | burial enr **1.78 > BE3D 1.54** | a naive "most-buried" caller matches/edges BE3D on AM-enrichment |
| **Beat AlphaMissense** (curated catalytic set, n=7 modeled) | **BE3D 2/7 vs AM 7/7** | **AM out-discriminates** on the exact catalytic residues |

**Honest verdict.** CHEK2 is the case BE3D handles *well*: it **passes the discrimination bar decisively**
(OR 4.3, p 7e-10; gap +0.24; **precision@10 = 100% functional**) and **concentrates on the folded kinase
ATP pocket, catalytic loop, and FHA** — not the disordered SCD. This is a genuine, VUS-relevant zero-shot
recovery from screen signal alone. **BUT** two honesty checks temper it: (1) the **base rate is high
(~39% of scored residues flagged)** — without an FDR/rank filter BE3D over-calls; (2) it **does not beat
AlphaMissense** on the curated catalytic set (2/7 vs 7/7) and only ties a trivial burial baseline on
AM-enrichment. → **BE3D on CHEK2 = complementary functional evidence (a real, orthogonal PS3-style signal),
not a standalone VUS classifier.**

## 3. FP / FN / the recovery (VUS-relevant detail)
- **THE RECOVERY (spatial-aggregation win):** the top robust cluster is the **catalytic loop 350–352**
  (P350/E351/**N352 mean-AM 0.998**) — the immediate 3D neighbourhood of the **catalytic aspartate D347**.
  **D347 itself is base-editing-UNREACHABLE** (a CBE cannot install the required codon change; no missense
  guide → unscored), yet BE3D flags its catalytic loop **via base-editable neighbours**. Likewise the **ATP
  cleft 303–308** (around K249) and **FHA R145/I146**. **TP-functional:** R145 (FHA-pThr), K249 (ATP-Lys)
  hit directly; the catalytic loop hit by proximity. This "recover the catalytic *neighbourhood* when the
  catalytic *residue* is unreachable" is BE3D's distinctive value over a per-variant predictor.
- **FALSE NEGATIVES (reachability-driven):** of 7 curated catalytic residues, **4 are base-editing-unreachable**
  (R117, D347, D368, T383 — no missense guide) and **1 reachable-but-missed** (S140). BE3D directly flags
  only 2/7. The misses are the **assay's** (CBE codon coverage), not the algorithm's — but they are real
  recall holes for VUS at the catalytic core.
- **LIKELY FALSE POSITIVES:** of 101 robust hotspots, **19 flagged likely-FP** — **13 in the low-pLDDT SCD
  N-terminus** (pLDDT<50, no fold) and **6 AM-benign**. These are avoidable with a pLDDT/AM gate.
- **Robust hotspot table:** 101 residues (all persist across r4/r6/r8): **2 TP-functional, 66
  likely-functional (AM-high, folded kinase/FHA), 14 novel-candidate, 19 likely-FP.** →
  `outreach/CHEK2/CHEK2_hotspots.tsv` (+ flags, mean-AM, ClinVar, pLDDT), `CHEK2_G2P.tsv`.

## 4. WHERE BE3D STRUGGLED (CHEK2)
1. **Over-calling / no FDR.** 38.8% of *scored* residues flagged at p<0.05 — and the count barely moves from
   p<0.05→p<0.001 (101→99): the randomization z-scores sit far in the tail, so the three nominal thresholds
   are **not** an effective stringency dial. A residue-level BH-FDR + a pLDDT gate would cut the 19 FP.
2. **Does not beat cheap predictors.** BE3D's AM-enrichment (1.54×) is **below** a naive most-buried caller
   (1.78×), and **AlphaMissense flags 7/7 catalytic residues vs BE3D's 2/7**. On a compact folded domain,
   BE3D adds *orthogonal phenotype evidence*, not superior discrimination.
3. **Reachability ceiling.** 4/7 catalytic residues (incl. the catalytic D347 and DFG D368) are simply not
   base-editable with a CBE → structurally invisible to the screen. BE3D should say so up front (PR #19).
4. **Placeholder-DSSP invalidates BE3D's own burial characterization** (RSA all >0.25) — the `pLDDT_RSA_scatter`
   and burial-enrichment plots are meaningless without a real DSSP; users may not notice.

## 5. Competitor comparison
| Tool | Signal | On CHEK2 | Advantage | Disadvantage |
|---|---|---|---|---|
| **BE3D** | CBE tiling → 3D LFC3D | recovers ATP pocket + catalytic loop + FHA from a survival screen; P@10=100% functional | orthogonal *functional* readout w/ direction; recovers catalytic *neighbourhood* of an unreachable residue | high base rate; ≤AM/burial on discrimination; CBE-reachability-limited |
| **AlphaMissense** | seq+structure DL | flags 7/7 catalytic, 50% of protein >0.564 | instant, proteome-wide, no experiment, near-saturated on catalytic core | predictor not measurement; no direction/mechanism; over-calls half the protein too |
| **ClinVar (this file)** | clinical | **benign = all synonymous**, P/LP mostly nonsense/splice | clinical truth | **no benign-missense at residue level here** → can't build a residue tolerant set from it |
| **3dhotspots/cBioPortal** | somatic recurrence | ~nothing (CHEK2 is a germline moderate-penetrance gene, not recurrently somatically mutated) | gold-standard for recurrent oncogenes | blind to CHEK2's use case |
→ **BE3D wins** = orthogonal functional/directional evidence on a germline VUS gene the recurrence tools can't
see; **AlphaMissense wins** = finer, cheaper residue discrimination on the folded catalytic core.

## 6. Cited change proposals (patient-impact-ordered)
1. **[highest value] Residue-level FDR (BH q-values) + a pLDDT confidence gate on hotspot calls.** Base rate
   39% of scored, and 13/19 CHEK2 FPs sit in the low-pLDDT SCD. Cited by §2 base-rate + §3 FP list.
   *Already partly covered by open PRs #20 (FDR q-values) and #23 (validation shortlist with pLDDT).*
2. **Emit a base-editing reachability report** (which functional residues are even installable by the screen's
   editor). Cited by §3 FN: 4/7 catalytic unreachable. *Covered by PR #19.*
3. **Warn when DSSP is a placeholder** so RSA/burial characterization isn't silently trusted (§4.4).
