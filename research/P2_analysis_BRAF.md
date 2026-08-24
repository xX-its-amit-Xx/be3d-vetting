# BE3D Run — BRAF RAF-inhibitor resistance (dabrafenib+cetuximab), 3D: pocket vs dimer interface

> Runner: Queen "EGFR/BRAF" (round 2). Real public base-editing screen → BE3D LFC3D + 3D clustering.
> **A0 DISCRIMINATION BENCHMARK leads.** The novel finding: the *reachable* BRAF resistance signal BE3D
> recovers sits at the **dimer interface (L505)**, not the ATP pocket — matching RAF-inhibitor biology.

## 0. Target, data, decision log
- **Gene/protein:** BRAF, UniProt **P15056** (766 aa; kinase domain ≈457–717). Canonical numbering ==
  clinical (verified: V600, T529 gatekeeper, R509 dimer, K483, G-loop G464/466/469, DFG D594/F595/G596).
- **Structure:** AlphaFold `AF-P15056-F1-model_v6` (WT, full-length) — **primary** (carries the αC/dimer
  interface and V600 region). Vemurafenib-bound 3OG7 noted as optional overlay (kinase-domain only; not
  required — AF has the dimer face).
- **Source:** Coelho, Dincer et al. 2024, *Nature Genetics*, DOI 10.1038/s41588-024-01948-8, **PMID
  39424923**. Sheet `ST2 BE z-scores`; BRAF = **7,176 guide rows**.
- **Drug arms (schema check + worker confirmation from paper full text):** the BRAF-V600E line is
  **HT29** (colorectal), with three arms — **DebCet = dabrafenib + cetuximab** (BRAF inhibitor +
  anti-EGFR; **the RAF-inhibitor resistance arm**), **Tram = trametinib** (MEK1/2i, downstream), **Pict =
  pictilisib** (PI3Ki, parallel). *There is NO vemurafenib/dabrafenib-monotherapy column in this dataset*
  — DebCet is the closest RAF-inhibitor arm and is treated as **primary**; Tram/Pict are contrast arms.
  **RESISTANCE = positive z.**
- **Reformatting:** `real_output/P2_BRAF/input/convert_braf.py` (KRAS/EGFR recipe, ABE+CBE pooled).
  **627 distinct missense positions, span 1–762; refAA-vs-P15056 check: 47/2509 mismatches (1.9%)** →
  numbering confirmed. Output → `input/BRAF_{DebCet,Tram,Pict}.tsv`.
- **Sweep (6 runs, all EXITCODE 0 / SUCCESS):** DebCet AF r4/r6/r8 mean + r6 sum (focused sweep); Tram r6
  and Pict r6 baselines. Driver `scratchpad/run_grid_braf.sh`.
- **QA (BE-QA):** KS H1 **D=0.26, p=5.6e-14** (DebCet) — knockouts deplete → BRAF essential in HT29 →
  screen **quality-accepted**.

## 1. ★ A0 DISCRIMINATION BENCHMARK (primary deliverable) ★
FUNCTIONAL = 15 validated residues split into **ATP-POCKET** {464,466,469,483,529,581,594,595,596,597,
600,601} and **DIMER-INTERFACE** {505,506,509}. TOLERANT = 20 ClinVar-benign/near-tolerant (N-term/linker/
surface). Baseline structure = DebCet AF r6 mean n1000, **positive (resistance) direction**.

| metric | **DebCet (dabrafenib+cetuximab)** |
|---|---|
| scored-positive residues N | 306 |
| **base rate** p<0.05 / 0.01 / 0.001 | **39% / 31% / 25%** |
| FUNCTIONAL residues *scored positive* | **3 / 15** |
| enrichment R/E @p<0.01 (R,E) | 2.12 (2, 0.94) |
| odds ratio; Fisher p @p<0.01 | 4.45; **p=0.23** ✗ |
| discrimination gap FUNC vs TOL @p<0.01 | 67% vs 43% (weak) |
| precision@10 / @20 | 0/10 ; 2/20 |
| burial baseline R/E (p) | 2.12 (0.23) — **BE3D does not beat it** |
| pLDDT of hits: median / frac<50 | 66 / 46% |

**Reading it honestly:**
1. **Base rate high (~39%)** → overlap alone is weak evidence.
2. **Base-editing inaccessibility is severe.** Only **3/15** functional residues get positive signal.
   **V600 (the V600E driver) has ZERO guide coverage** — GTG→GAG is a transversion no BE installs; and
   **R509 (dimer interface) also has zero coverage.** Gatekeeper **T529** gets only LOF-depleting T529A;
   DFG D594/F595/G596/L597 unreachable. So the canonical pocket-resistance and driver residues are simply
   absent from the reachable variant space.
3. **The formal A0 test does NOT clear chance** (R/E=2.1, p=0.23; base rate 39%; precision@10=0; does not
   beat the burial baseline). The tiny reachable functional set (3) makes formal significance unattainable.
4. **BUT the two functional positive hits are both DIMER-INTERFACE and biologically real:** **L505 and
   R506.** L505 is reachable (CTC→TTC, C>T) → the screen installs **L505F (LFC +2.02)**, and **L505H is a
   validated ACQUIRED vemurafenib-resistance mutation** (PMID 25515853). **These 505/506 hits are ROBUST
   across all three arms** (DebCet, Tram, Pict). Precision@10 is 0 because the top |LFC3D| ranks are
   N-terminal CR1/CR2 residues (364/365/403/404) — passenger/FP-risk.
5. **Negative (essential) direction** recovers the pocket that positive can't: gatekeeper **T529**, G469,
   catalytic-loop **N581** are neg-significant (reachable pocket edits kill the kinase) — essentiality,
   not resistance.

## 2. ★ NEW question — ATP pocket vs dimerization interface ★
**Answer: the only *reachable* resistance geometry BE3D recovers is the DIMER-INTERFACE / αC-helix face,
not the ATP pocket** — consistent with the established mechanism that RAF-inhibitor resistance is
dimerization-driven (Poulikakos 2011 PMID 22113612; Rajakulendran 2009 PMID 19727074). Specifics:
- **Dimer/αC face lit up, robustly:** L505, R506, K507, T508 (the R-spine αC residues at the dimer face)
  + an extended contiguous αC/β3 patch **492–508, 521–533, 543–550** across all arms.
- **ATP pocket dark in the positive direction:** V600E and gatekeeper/DFG resistance alleles are
  base-editing-inaccessible; the pocket only appears in the *negative/essential* direction.
- **Caveat (do not overclaim):** this rests on **2 formally-validated residues (505/506)** plus a
  contiguous novel patch; enrichment is NOT statistically significant (p=0.23). It is a **biologically
  coherent, hypothesis-generating** result, not a proven discrimination win. The αC patch 521–550 is
  NOVEL-candidate (structurally in the dimer/αC region, but contiguity ⇒ possible bystander artifact).

## 3. Top hotspots (robust ≥3/4 DebCet AF runs, positive) — TP / novel / FP
Tables: `outreach/BRAF/BRAF_{DebCet,Tram,Pict}_hotspots.tsv`.
- **TP (validated, dimer interface):** **L505, R506** — all 3 arms. L505H = acquired vemurafenib R.
- **NOVEL-candidate (αC-helix / dimer region, plausible, unvalidated):** 492/493, 507/508, 521–524,
  532/533, 543–550 (β3–αC, gatekeeper-adjacent), 688–695 (αI). Flag: contiguous runs ⇒ bystander risk.
- **LIKELY FALSE POSITIVE:** N-terminal CR1/CR2 and inter-domain residues (top |LFC3D| 364/365/403/404),
  the 607–609 patch (pLDDT<50), and A689 (ClinVar near-tolerant) — ~55/80 robust hits sit outside the
  kinase domain.
- **FALSE NEGATIVES (validated, missed):** **V600E (driver — zero coverage), R509 (dimer — zero
  coverage), T529N/M (gatekeeper R — only LOF T529A reachable), DFG D594/F595, G-loop G466/G469, K601E**
  — base-editing-inaccessible or LOF-only. Dominant miss category; not a BE3D fault.

## 4. Competitor note
- **3dhotspots.org / cBioPortal:** report BRAF **V600E** (and G466/G469/L597/K601) by recurrence — the
  activation driver, orthogonal to per-drug escape, and the residue BE cannot reach. They win on the
  canonical driver by construction.
- **AlphaMissense:** drug-agnostic; flags the pocket/DFG high regardless of drug; cannot see the
  dimerization-driven resistance logic.
- **BE3D's unique edge:** per-drug escape geometry — here it re-finds the **dimer-interface L505**
  acquired-resistance residue from screen signal alone, which a recurrence- or pathogenicity-based tool
  (both V600-centric) would not surface. Edge remains **gated by base-editing reachability**.

## 5. Ground-truth citations (worker-sourced)
V600E Wan 2004 PMID 15035987; gatekeeper T529 Whittaker 2010 PMID 20538618; **L505H acquired vemurafenib-R
PMID 25515853**; R509 dimer interface Rajakulendran 2009 PMID 19727074; intact-interface requirement
Röring 2012 PMID 22735454; paradox/p61 splice resistance Poulikakos 2011 PMID 22113612; class1/2/3
taxonomy Yao 2015 PMID 26343582, Yao 2017 PMID 28783719; DFG/paradox Heidorn 2010 PMID 20141835. Tolerant:
ClinVar likely-benign M117T/S136A/A145S/T313A/I208V (I208V PMID 30355600) + N-term/linker near-tolerant.

## 6. BE3D issues → BE3D_IMPROVEMENTS.md
Same as EGFR: (a) no base-editing reachability report (V600/R509 have ZERO guides — a coverage/reachability
summary per curated residue would make the false negatives explicit); (b) ~39% positive base rate with no
warning; (c) low-pLDDT and N-terminal disordered residues dominate top hits with no structural gate; (d)
contiguous-neighborhood runs (521–550) unflagged for bystander risk.
