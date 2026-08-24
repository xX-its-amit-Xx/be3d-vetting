# BE3D → per-user package — MCL1 (BH3-mimetic / apoptosis drug-development lab)

Synthesizer "MCL1-outreach". Date 2026-08-24. Built on the completed P2 run
(`real_output/P2_MCL1/`, decision log `research/P2_run_MCL1.md`). No BE3D re-run.
Analysis scripts: `real_output/P2_MCL1/a0_discrimination.py`, `build_hotspot_table.py`
(system py3.14, pandas 2.3.3 / scipy 1.17.1). Deliverables in `outreach/MCL1/`.

Target: MCL1 (Q07820, 350 aa). Screen: Hanna et al. 2021 Cell CBE tiling in MELJUSO,
3 arms — **dropout** (plain essentiality, QA-failed/underpowered), **A1331852**
(BCL-xL inhibition makes MCL1 essential → dependency, NEG LFC3D), **combo A1331852+S63845**
(MCL1-inhibitor resistance, POS LFC3D). Structures: AF Q07820; 3PK1/3MK8 (BH3 groove);
6QB4 (HVN/S63845-class inhibitor pocket).

---

## 1. A0 DISCRIMINATION BENCHMARK (primary; the credibility gate)

### Base rate — high among tested residues, so raw overlap is weak evidence (stated honestly)
BE3D only computes LFC3D for residues in an edited neighbourhood, so the honest denominators
differ. `n_scored` = residues with a directional LFC3D value.

| arm | n_scored | p<0.05 nsig | /n_scored | /350 aa | p<0.01 /n_scored | p<0.001 /n_scored |
|---|---|---|---|---|---|---|
| dropout | 47 | 18 | **38.3%** | 5.1% | 36.2% | 36.2% |
| A1331852 dependency | 113 | 51 | **45.1%** | 14.6% | 39.8% | 37.2% |
| combo / S63845 resistance | 95 | 41 | **43.2%** | 11.7% | 30.5% | 14.7% |

**Reading:** among residues BE3D can actually score, it flags ~40–45% at p<0.05 — a HIGH base
rate. "BE3D recovered the groove" is therefore weak on its own; the case rests on the
enrichment / discrimination-gap tests below, not on raw overlap. (Over the whole 350-aa protein
the rate is a more modest 5–15%, because ~2/3 of residues are never scored — see FN §2.)

### Discrimination gap — hit-rate on FUNCTIONAL vs TOLERANT (p<0.05, scored-domain universe)
FUNCTIONAL = BH3 groove (3PK1/3MK8 union ≤5 Å, n=31) ∪ S63845/HVN pocket (6QB4 ≤5 Å, n=18) = 36.
TOLERANT = scored domain residues, solvent-exposed, distal (>8 Å from any groove/pocket set).

| arm | FUNCTIONAL hit-rate | TOLERANT hit-rate | gap | ratio |
|---|---|---|---|---|
| **combo / S63845 resistance** | **11/11 = 100%** | 1/8 = 12.5% | **+87.5%** | **8.0×** |
| **A1331852 dependency** | **12/13 = 92.3%** | 6/20 = 30.0% | **+62.3%** | **3.1×** |
| dropout | 4/8 = 50% | 5/8 = 62.5% | **−12.5%** | 0.8× |

### Enrichment vs chance (Fisher exact, one-sided; obs R vs E=K·M/N)
| arm | target | obs R | E | R/E | OR | p |
|---|---|---|---|---|---|---|
| combo/S63845 | groove | 11 | 6.44 | 1.71 | ∞ | **7.9e-4** |
| combo/S63845 | pocket-only | 0 | — | — | 0 | ns |
| A1331852 | groove | 11 | 6.05 | 1.82 | ∞ | **5.6e-4** |
| A1331852 | groove+pocket | 12 | 7.15 | 1.68 | 14.9 | **2.0e-3** |
| A1331852 | pocket-only | 1 | 1.1 | 0.91 | 0.81 | ns |
| dropout | groove | 4 | 4.48 | 0.89 | 0.7 | ns |

### Precision@K (top-K by |LFC3D|, domain-ranked; chance = M/N)
| arm | P@10 | P@20 | chance |
|---|---|---|---|
| combo/S63845 | 30% | **45%** | 27% |
| A1331852 | 30% | 25% | 22% |
| dropout | 40% | 30% | 32% |
- **Honest nuance:** P@10 is only modest because the LARGEST-magnitude hits are OFF-groove
  (whole-protein top-10 for combo = 324,325,216,215,180,178,177,219,220,217 → only 215/216/220
  are groove; for A1331852 the top-10 is 217–221 groove-core + 324–328 C-term). The strongest raw
  signals sit in the off-groove **stability compartment**, not the groove — this is the two-compartment
  finding, not a failure, but it does mean magnitude-ranking alone under-weights the groove.

### Beat a trivial baseline
- **Burial:** the pipeline's RSA/ACC columns are saturated (min RSA 0.35, ACC all 100 — unusable);
  used BE3D's categorical `exposure`. "Call all buried (medburied)" gives FUNCTIONAL enrichment
  OR=0.75 (combo, ns) / 2.0 (A1331852, ns) — i.e. **burial does NOT predict the groove** (the BH3
  groove is a solvent-exposed surface feature). BE3D's functional enrichment **holds among
  exposed-only residues** (combo OR=∞ p=7.8e-3; A1331852 OR=6.4 p=0.08) → the signal is not a
  burial proxy.
- **AlphaMissense (competitor as baseline):** AM calls **86% of the folded BCL domain (209–320)
  pathogenic** (mean 0.805) — "call all AM-high" would flag nearly the whole domain (huge base
  rate, no phenotype specificity). BE3D isolates the ~9 escape residues *within* that uniformly-high
  region. AM cannot separate a resistance/escape residue from a generically-essential one (§3).

### VERDICT (honest)
**BE3D DISCRIMINATES on the two drug-relevant arms and OVER-CALLS on the underpowered dropout arm.**
- Combo/S63845 resistance and A1331852 dependency: functional hit-rate 92–100% vs tolerant 12–30%
  (3–8× gap), groove enrichment OR=∞, p<0.001, robust across 5 structures, and it beats both the
  burial and (in specificity) the AlphaMissense baselines. This is genuine discrimination, not
  "everything is a hotspot" — despite the high raw base rate.
- Dropout arm: NEGATIVE discrimination gap (−12.5%), no enrichment; consistent with the QA failure
  (1 nonsense guide). Its top hotspots (surface helix 295–303) are treated as noise/FP.

---

## 2. FALSE POSITIVES / FALSE NEGATIVES

### Robust hotspot classes (full table: `outreach/MCL1/MCL1_hotspots.tsv`)
- **TP (validated groove/pocket), structure-independent (sig in all 5 combo or all 5 dep runs):**
  - S63845 resistance (+): **R215, V216, G217, D218, G219, V220, Q221** (α3), **S255, D256, G257,
    V258, T259, S269** (α4/shelf), F254 (n=4). Raw strongest single edit T212I (+2.44).
  - A1331852 dependency (−): **G217–R222 core** (LFC3D ≈ −3.4), α8 edge **E317/F318/F319/H320/V321**
    (in 3PK1 groove-shell 5 Å).
  - These found residues carry AlphaMissense-pathogenic support (R215 0.94, V216 0.96, V220 0.93,
    F254 0.99, D256 0.99, V258 0.81). The buried pocket residues BE3D MISSED are also AM-high
    (T266 0.98, F270 0.99, anchor R263 0.997) — see false negatives in §2.
- **NOVEL-CANDIDATE — the off-groove stability compartment (flag honestly; real biology OR passenger):**
  - **α2: R176/Q177/S178/E180** (resistance, +, n=2–4); **C-terminal: E322/L324/E325/G326/G327/I328/R329**
    (+resistance and −dependency, n=2–5); **N-terminal PEST: R78/E85/R95** (+, n=3–4) and scattered
    123–138.
  - **Plausible mechanism (cited):** MCL1's disordered N-terminal PEST region (~1–170) carries the
    FBW7 phosphodegron (GSK3-primed S159/T163; Wertz et al. Nature 2011, PMID 21358673), plus
    MULE/HUWE1 and βTrCP sites — mutations that stabilise MCL1 / raise its abundance could confer
    BH3-mimetic resistance WITHOUT touching the drug pocket (all off-groove sites are >8 Å from HVN).
  - **Why it could be a passenger/artifact (honest caveat):** these lie in low-pLDDT disordered AF
    regions with few 3D neighbours, so LFC3D there is essentially a per-guide (sequence-level) call
    with weak spatial support; several are single-run (n=1). Treat as **hypotheses**, not clusters.
- **LIKELY FALSE POSITIVE:** the **α6 surface patch 295–303 (R300/T301/K302)** — no known function,
  solvent-exposed helix distal to the groove, called mainly by the QA-failed dropout arm (and
  moderately by A1331852). Best explained as screen noise / a surface bystander cluster, not
  biology. Also low-confidence: dropout N-terminal 58–61, 70, 122.

### FALSE NEGATIVES — functional groove residues BE3D missed (and WHY)
Of the 36 FUNCTIONAL (groove+pocket 5 Å) residues:
- **19 were NEVER scored** (no CBE-reachable disruptive edit and/or no guide coverage):
  224, 227, 228, **231 (M231, P2 wall)**, 234, 235, 246, 248, 249, **250 (M250)**, 252, 253, 260,
  261, 262, **270 (F270, P2 floor)**, **271 (G271)**, **274 (V274)**, 305. This is the intrinsic
  **base-editing blind spot** — CBE only makes C→T (mostly missense to a limited set), so the codon
  changes needed to disrupt these hydrophobic-pocket side chains aren't reachable.
- **1 scored but never significant:** L290.
- **The conserved anchor R263** (the residue whose salt bridge to the BH3 aspartate / drug
  carboxylate defines the groove; AM 0.997) was only weakly caught (1 run) — effectively a near-miss.
- Net: BE3D recovers the α3/α4 **rim** of the groove (the CBE-accessible face) strongly, but misses
  the buried P2-pocket floor — a coverage limitation, not a scoring failure. **AlphaMissense catches
  these (it needs no guide); a functional screen cannot edit them.** Complementary, not redundant.

**Precision/recall FEEL:** of ~30 robust hotspots — ~13 are validated groove/pocket TPs, ~12 are
off-groove novel-candidate (stability) hypotheses, ~5 look like false positives (295–303 patch +
dropout N-term). Recall of the groove rim is high; recall of the buried pocket is low (edit access).

---

## 3. COMPETITOR COMPARISON (B)

| tool | input signal | what it finds on MCL1 | advantage | disadvantage |
|---|---|---|---|---|
| **BE3D** (this run) | base-editing screen LFC → 3D randomization null | escape/dependency rim of BH3 groove (R215/V216/V220/D256/V258…) + a novel off-groove **stability** compartment; separates the two | functional-phenotype signal; sees drug-resistance & stability escape that no static predictor or recurrence tool can; editor-specific escape geometry | high base rate; blind to CBE-inaccessible buried pocket (F270/M231); underpowered dropout arm over-calls |
| **AlphaMissense** (Cheng 2023, PMID 37733863) | sequence+structure pathogenicity predictor | **86% of the folded domain 209–320 pathogenic** (mean 0.805); N-term PEST mostly benign (0.33) | zero-shot, no screen needed; covers the buried pocket residues BE3D misses (R263 0.997, F270 0.99) | cannot separate a **resistance/escape** residue from a generically-essential one — phenotype-agnostic; flags the whole groove uniformly |
| **3dhotspots.org / cancerhotspots** (Gao 2017, PMID 28115009) | somatic mutation **recurrence** | **NOTHING** — live API `?hugoSymbols=MCL1` → `[]`; cancerhotspots MCL1 = 0 residues | decade-long benchmarking where recurrence exists | blind to MCL1 entirely |
| **cBioPortal 3D hotspots / HotMAPS / OncodriveCLUST** | recurrence clustering | nothing (MCL1 not a driver — Bailey 2018, PMID 29625053) | robust for point-mutation drivers | same recurrence blind spot |
| **ProTiler-Mut** (Cell Systems 2026, DOI 10.1016/j.cels.2026.101651) | base-editing tiling → residue/substructure/SoF | (conceptual) would also map the groove from screen signal; adds separation-of-function class | closest peer; SoF classification; multi-condition | 3D-RRA/segmentation vs BE3D's randomization-null + arm-specific contrast; not run here |

**The MCL1-specific BE3D advantage (state explicitly):** MCL1 is one of the most frequently
**AMPLIFIED** genes in cancer (1q21; Beroukhim et al. Nature 2010, PMID 20164920) but is **NOT a
recurrent point-mutation gene** — 0 somatic 3D hotspots, 0 cancerhotspots residues, not a Bailey
2018 driver. So every recurrence-based 3D-hotspot tool returns nothing for MCL1. This is the clean
context where a **functional-screen** tool (BE3D / ProTiler-Mut) is the ONLY 3D-hotspot method that
produces a map at all — and where BE3D adds what AlphaMissense structurally cannot: which groove
residues are drug-**escape** residues, and a second **stability-escape** axis.

---

## 4. OUTREACH (C, D) — files in `outreach/MCL1/`
- `MCL1_hotspots.tsv` — robust hotspot table, residue/AA/direction/LFC3D/n_robust_runs/groove·pocket
  membership/**TP·NOVEL·FP flag**.
- `MCL1_S63845resistance_G2P.tsv`, `MCL1_A1331852dependency_G2P.tsv` — G2P interactive-3D inputs.
- `MCL1_email.md` (copyable), `MCL1_brief.md` (1-page forwardable), `_false_negatives.txt`.
- Frame: dependency AND drug both localise to the groove → S63845 is a true orthosteric BH3-mimetic;
  the off-groove stability route → implications for degraders; the 9 groove escape residues to de-risk.

**Email readiness:** READY. Positive controls recovered zero-shot (groove escape rim, structure-
independent across AF/3PK1/3MK8/6QB4) + honest FP/FN + the amplification-not-mutation competitor
gap give an honest, non-overclaiming pitch to a BH3-mimetic / degrader program.
