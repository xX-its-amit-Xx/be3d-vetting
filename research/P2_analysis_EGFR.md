# BE3D Run — EGFR TKI resistance (gefitinib 1st-gen vs osimertinib 3rd-gen), 3D

> Runner: Queen "EGFR/BRAF" (round 2). Real public base-editing screen → BE3D LFC3D + 3D clustering.
> **A0 DISCRIMINATION BENCHMARK leads.** Honest verdict up front: on the *resistance* axis BE3D is
> throttled by base-editing chemistry, not by BE3D itself.

## 0. Target, data, decision log
- **Gene/protein:** EGFR, UniProt **P00533** (1210 aa; kinase domain ≈712–979). Canonical UniProt
  numbering == clinical numbering (verified: T790, L858, C797, G719, G724, L718, L792, L861, E746–A750
  all land on the expected residues).
- **Structure:** AlphaFold `AF-P00533-F1-model_v6` (WT, full-length) — **primary**. Osimertinib-bound
  crystals **6LUD** and **4ZAU** (chain A) run as overlays for pocket-distance annotation.
- **Source:** Coelho, Dincer et al. 2024, *Nature Genetics*, "Base editing screens define the genetic
  landscape of cancer drug resistance." DOI 10.1038/s41588-024-01948-8, **PMID 39424923**. Sheet
  `ST2 BE z-scores`; EGFR = **7,030 guide rows**.
- **Drug arms (schema check, per instruction):** the z-score columns are NOT gene-specific — every gene
  carries all cell-line×drug arms. The clinically relevant line for EGFR is **PC9** (EGFR-mutant NSCLC):
  `L2FC_PC9_Gefit_plasmid_average_zscore` (**gefitinib, 1st-gen**) and `L2FC_PC9_Osim_…` (**osimertinib,
  3rd-gen**) — run as two separate screens. **RESISTANCE = positive z** (enrichment under drug).
- **Reformatting:** `real_output/P2_EGFR/input/convert_egfr.py` (mirrors the KRAS recipe: category from
  `most_severe_consequence`; edits reconstructed from `zip(Amino_Acid_Position, Edited_AA, New_AA)`;
  ABE+CBE pooled). 627→ (EGFR) distinct missense positions. Converter output → `input/EGFR_{Gefit,Osim}.tsv`.
- **Sweep (18 runs, all EXITCODE 0 / RUN_COMPLETED SUCCESS):** AF r4/r6/r8/r10 mean, r6 sum, r6 n500,
  r6 atom, + 6LUD & 4ZAU r6, × {Gefit, Osim}. Driver `scratchpad/run_grid_egfr.sh`.
- **QA (BE-QA):** KS H1 (nonsense/splice vs no-mut/silent) **D=0.37, p=4.2e-52** (Osim) — knockouts
  DEPLETE, i.e. EGFR is essential in PC9. Unlike the KRAS resistance-only screen, the standard QA gate
  **passes**; the screen is quality-accepted.

## 1. ★ A0 DISCRIMINATION BENCHMARK (the primary deliverable) ★
Ground-truth sets curated INDEPENDENTLY of BE3D (worker-sourced, PMIDs in §4). FUNCTIONAL = 25 validated
resistance/driver/catalytic residues; TOLERANT = 16 ClinVar-benign/germline positions (ectodomain + C-tail).
Baseline structure = AF r6 mean n1000, **positive (resistance) direction**.

| metric | **Gefitinib (1st-gen)** | **Osimertinib (3rd-gen)** |
|---|---|---|
| scored-positive residues N | 272 | 304 |
| **base rate** p<0.05 / 0.01 / 0.001 | 32% / 26% / 19% | **40% / 31% / 24%** |
| FUNCTIONAL residues *scored positive* | **7 / 25** | **5 / 25** |
| enrichment R/E @p<0.01 (R,E) | **3.28** (6, 1.83) | 1.29 (2, 1.55) |
| odds ratio; Fisher p @p<0.01 | **18.5; p=1.5e-3** ✅ | 1.5; p=0.49 ✗ |
| discrimination gap FUNC vs TOL @p<0.01 | **86% vs 37%** | 40% vs 37% (none) |
| precision@10 | 1/10 | 0/10 |
| burial baseline R/E (p) | 2.19 (0.078) | 1.94 (0.17) |
| pLDDT of hits: median / frac<50 | 63 / 32% | 57 / 45% |

**Reading it honestly:**
1. **Base rate is high (~30–40%).** Any bare "recovered residue X" overlap is weak evidence; only formal
   enrichment counts.
2. **Base-editing inaccessibility dominates.** Only **5–7 of 25** validated functional residues receive
   ANY positive signal. **~18/25 are unreachable by CBE/ABE**: T790M (ACG→ATG needs the right transition
   but the screen installed only T790A), **C797S** (needs a transversion → only C797R/G/Y made, all
   strongly *depleting* LOF, LFC −3…−5), **L858R** (CTG→CGG transversion → only L858P made), **exon19del
   E746–A750** (not a point edit), G719, K745, D855. This is the DNMT3A-R882H failure mode at scale.
3. **Gefitinib PASSES the A0 test; osimertinib FAILS it.** Gefit shows real enrichment (OR 18.5,
   p=1.5e-3) with a genuine discrimination gap (86% vs 37%) and beats the burial baseline. Osim shows no
   enrichment and no discrimination (tolerant residues hit as often as functional). See §2 for why.
4. **Precision@10 ≈ 0.** The strongest |LFC3D| signals sit in the **disordered C-terminal tail**
   (res 1004–1019, 1096–1098, 1184–1186; pLDDT<50) — the autophosphorylation tail, **likely false
   positives** (many guides, low structural confidence), not pocket resistance.
5. **Negative (essential) direction is the real structural signal but it is NOT resistance.** In the neg
   direction 20/25 functional residues ARE scored and the ATP-pocket core (L718, G724, T790, L792, C797)
   lights up, with a clean discrimination gap (FUNC 30–45% vs TOL 12.5%). This recovers **kinase
   essentiality** (reachable pocket edits kill the kinase), orthogonal to the resistance question.

**Verdict:** On the resistance axis, BE3D is **weakly positive for gefitinib, negative for osimertinib**,
and its top signals are tail false-positives — but the binding constraint is base-editing chemistry, not
the algorithm. BE3D's negative direction is a solid essential-pocket detector.

## 2. What BE3D DOES recover (the reachable signal) & the 1st-vs-3rd-gen contrast
The Gefit enrichment is driven almost entirely by the **αC-helix / exon20 region**: robust positive hits
at **S768, V769, D770, N771, H773** (+ modifier E709). These are real: **exon20 insertions and S768I are
validated INTRINSIC resistance loci to 1st-gen TKIs** (Yasuda 2013 PMID 24353160) — exactly the arm where
they enrich. So BE3D re-finds a bona-fide gefitinib-resistance compartment zero-shot (TP).

**NEW question (do 1st- vs 3rd-gen hotspots occupy different 3D compartments?) — honest answer:**
partially unanswerable. The classic compartments the question targets — the **C797/covalent site** and the
**ATP-pocket rim** (T790/L718/G724/L792) — are **base-editing-inaccessible**, so BE3D sees neither for
either drug. The one reachable resistance compartment, the **αC/exon20 face (768–773)**, is **shared**
(Gefit robustly, Osim weakly: 769/770 only). Osimertinib simply has fewer reachable escape residues in
this screen, which is why its A0 enrichment collapses. A contiguous C-lobe patch **806–812** (and 929–933
for Osim) recurs robustly but is **NOVEL-candidate with FP risk** (a run of adjacent residues → possible
overlapping-guide/bystander artifact; activation-segment-adjacent but unvalidated).

## 3. Top hotspots (robust ≥5/7 AF runs, positive direction) — TP / novel / FP
Full tables: `outreach/EGFR/EGFR_{Gefit,Osim}_hotspots.tsv`.
- **TP (validated):** V769, D770 (exon20 anchors) — Gefit & Osim; S768, N771, H773, E709 — Gefit.
- **NOVEL-candidate (kinase, structurally plausible, unvalidated):** 806–812 (C-lobe/activation-segment
  adjacent), 929–933 (αG/substrate lobe, Osim), 735/738 (β2–β3 near G-loop, Osim). Flag: contiguity ⇒
  possible bystander artifact.
- **LIKELY FALSE POSITIVE:** the bulk of robust hits — C-terminal tail >979 (1004–1019, 1096–1098,
  1184–1186), 41–51 of ~50–65 robust hits, pLDDT<50, disordered.
- **FALSE NEGATIVES (validated, missed):** T790M, C797S, L718Q, G724S, L792H, G796R, L858R, exon19del,
  G719, D761Y, T854A — **base-editing-inaccessible** (require transversions/indels CBE/ABE cannot install),
  the dominant miss category; not a BE3D fault.

## 4. Competitor note
- **3dhotspots.org / cBioPortal:** report EGFR's recurrence codons — L858R, exon19del, T790M, G719 —
  i.e. **activation/driver** sites, orthogonal to per-drug escape and (crucially) the very residues BE
  can't reach. They would "win" on the canonical hotspots by construction.
- **AlphaMissense:** drug-agnostic; scores the whole kinase pocket high regardless of drug and needs no
  screen. It cannot separate gefitinib- from osimertinib-specific escape.
- **BE3D's unique edge:** per-drug escape geometry — demonstrated here by recovering the **exon20/S768
  gefitinib-resistance face** from screen signal alone. But that edge is **gated by base-editing
  reachability**, a limitation the competitors don't share.

## 5. Ground-truth citations (worker-sourced)
Activating: exon19del/L858R Lynch 2004 PMID 15118073, Paez 2004 PMID 15118125; L861Q/L792/G796 PMID
29857056; S768I/E709 doi 10.3389/fphar.2022.976731. Catalytic/pocket: K745/T790/D855 Yun 2007 PMID
17349580, Stamos 2002 PMID 12196540; M793 PMID 12196540. Resistance: T790M Kobayashi 2005 PMID 15728811 /
Pao 2005 PMID 15737014; C797S Thress 2015 PMID 25939061; L718Q/G724S Oztan 2017 PMID 28838405; T854A Bean
2008 PMID 19010870; D761Y PMID 17085664; exon20ins Yasuda 2013 PMID 24353160. Tolerant: R521K PMID
21896992 + ClinVar benign ectodomain/C-tail positions.

## 6. BE3D issues → appended to BE3D_IMPROVEMENTS.md
(a) No base-editing REACHABILITY report: BE3D silently scores only the substitutions a BE can install, so
a user comparing to clinical hotspots sees unexplained false negatives (here 18/25). A "target allele not
reachable by CBE/ABE" annotation would prevent misreading. (b) High positive base rate (~40% of scored)
with no built-in base-rate warning. (c) Disordered/low-pLDDT residues (C-terminal tail) dominate the top
positive hits with no pLDDT gate on the hotspot call. (d) Contiguous single-neighborhood runs (806–812)
not flagged as possible bystander artifacts.
