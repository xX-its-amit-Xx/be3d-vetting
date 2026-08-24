# BE3D Run — Analysis: PIK3CA (PI3Kα) drug-resistance base-editing screen (3D)

> Runner: Queen "PI3K/AKT" (round 2). Real public base-editing screen (Coelho/Dincer 2024) →
> BE3D LFC3D + 3D clustering. **A0 discrimination benchmark leads.** Honest verdict up front.

## 0. Target, data, and the drug-arm reality (read this first)
- **Gene / protein:** PIK3CA (UniProt **P42336**, p110α, 1068 aa). Domains: ABD 1–108, RBD 187–289,
  C2 330–487, helical 517–694, kinase 797–1068.
- **Structure:** AlphaFold `AF-P42336-F1-model_v6` (WT, canonical numbering). **Primary.** (Optional
  alpelisib-bound 8TS9/6PYS not needed — the on-target arm carries no signal; see below.)
- **Source:** Coelho, Dincer et al. 2024 *Nat Genet* (DOI 10.1038/s41588-024-01948-8), suppl. table S4
  sheet `ST2 BE z-scores`; **4,202 PIK3CA guide rows** (1,357 missense guides → 836 distinct missense
  positions, span 1–1066). Score = per-drug `L2FC_<cell>_<drug>_plasmid_average_zscore`, **positive =
  resistance/enrichment**. Converter: `scratchpad/convert_pi3k.py`; TSVs in `scratchpad/pi3k_tsv/PIK3CA/`.
- **HONEST DRUG CAVEAT (important):** the requested **alpelisib is NOT in this dataset.** The only
  on-target PI3K inhibitor present is **pictilisib** (GDC-0941, pan-PI3K ATP-competitive, HT29 arm).
  Available arms (all tile PIK3CA): HT29 (dab+cetux, trametinib/MEK, pictilisib/PI3K), H23 (sotor/adag,
  KRAS), PC9 (osim/gefit, EGFR), MHHES1 (olap/nirap, PARP).
- **Signal scan (missense guides with Z>3):** pictilisib **0**, trametinib **8**, gefitinib **17**,
  osimertinib 7, dab+cetux 5, KRAS-i 2–4, PARP-i 0. **Pictilisib (the on-target PI3K inhibitor) produces
  essentially no resistance signal** (2 guides Z>2, max Z 2.1). The strong signal is under **bypass
  drugs** — MEK/EGFR inhibitors, escaped by PI3K-pathway activation.
- **Arms run here:** **Trametinib (primary, MEK-bypass — strong signal, same HT29 cell line as
  pictilisib)** and **Pictilisib (on-target contrast — null)**. This pair *is* the answer to the NEW
  question: does resistance sit in the ATP/drug pocket (on-target escape) or on the helical/regulatory
  activation hotspots (bypass)? See §4.
- **Numbering:** screen refAA vs P42336 sequence match 7904 / 8048 (98.2%) — clean.

## 1. BE-QA (screen quality)
Cases = Nonsense+Splice-donor (267), controls = No Mutation+Silent (2504). **KS D=0.239, p=1.4e-12;
Mann-Whitney p=2.2e-12** (Trametinib arm) — **strongly significant**, unlike KRAS (which was null).
Knockout guides deplete while activating missense enrich → the screen carries real, directional signal.
**ACCEPT.** All interpretation uses the positive (resistance) direction.

## A0. DISCRIMINATION BENCHMARK (primary, non-negotiable)
Ground-truth sets (curated INDEPENDENTLY of BE3D, cited — `scratchpad/groundtruth_pi3k.json`):
- **FUNCTIONAL-ACTIVATING** (M=12): 88, 345, 418, 420, 453, 542, 545, 726, 1043, 1047, 1049, 1065
  (helical/C2/regulatory + kinase-domain oncogenic hotspots; Samuels 2004 PMID 15016963, Gymnopoulos
  2007 PMID 17376864, Huang 2007 PMID 18079394, Mandelker 2009 PMID 19805105).
- **FUNCTIONAL-POCKET/CATALYTIC** (M=13): ATP/drug cleft + catalytic — K802, S774, I800, D810, Y836,
  V851, Q859, M922, I932, D933 (DFG), DRH 915–917 (Walker 1999 PMID 10580505).
- **FUNCTIONAL-UNION** M=25 (base rate M/N = **0.023**).
- **TOLERANT** (n=203): AlphaMissense-benign residues minus functional (well-powered proxy).

**Results (AlphaFold, r6, mean, n1000, positive direction):**

| arm | thr | n_sig | **base rate** | UNION R/E (p) | **ACT R/E (p)** | POCKET R/E | gap hitF / hitT |
|---|---|---|---|---|---|---|---|
| **Tram** | p<0.05 | 152 | 0.142 | 1.69 (0.13) | **2.93 (0.019)** | 0.54 | 0.24 / 0.11 |
| **Tram** | p<0.01 | 124 | 0.116 | 1.72 (0.15) | **3.59 (0.008)** | 0.0 | 0.20 / 0.10 |
| **Tram** | p<0.001 | 80 | 0.075 | 1.60 (0.29) | 3.34 (0.054) | 0.0 | 0.12 / 0.06 |
| Pict | p<0.01 | 138 | 0.129 | 0.93 (0.65) | 1.93 (0.19) | 0.0 | 0.12 / 0.13 |
| Pict | p<0.001 | 111 | 0.104 | 1.15 (0.49) | 2.41 (0.12) | 0.0 | 0.12 / 0.11 |

**Precision@10 (Tram, r6):** top-10 by |LFC3D| = `80, 468–470, 544–548, 969`. **5 of the top 10 are the
E542/E545 helical-hotspot 3D cluster (544–548)** — the single strongest 3D peak sits exactly on the
validated helical activating hotspot. prec@10 on exact functional members = 0.1 (only 545 counts as an
exact set member; 542 is #11–20), but the strongest neighborhood is unambiguously the helical hotspot.

**Beats a trivial baseline?** **NO (for generic functional recovery).** AlphaMissense-high (mean AM≥0.90)
gives **R/E 2.29, OR 6.66, p 1.2e-5** for the same FUNCTIONAL-UNION — highly significant and stronger
than BE3D's union enrichment (R/E 1.7, p=0.15, n.s.). Burial (top-tertile CN) is useless here (R/E 0.65).

### A0 verdict (PIK3CA) — HONEST
**QUALIFIED PARTIAL SIGNAL, not a clean win.**
1. **High base rate** (7–14% of residues flagged) → "everything is a hotspot" risk; any single overlap
   is weak evidence. Say so plainly.
2. **Union enrichment is NOT significant** (R/E~1.7, p≈0.15). BE3D does **not** beat AlphaMissense at
   recovering functional residues generically.
3. **BUT the signal is real and specific where it counts:** the resistance hotspots concentrate
   **selectively on the ACTIVATING helical/regulatory hotspots (ACT R/E 2.9–3.6, p=0.008–0.02)** and are
   **completely absent from the ATP/drug pocket (POCKET R/E = 0.0)**. The strongest 3D peak IS E542/E545.
4. This is information AlphaMissense **cannot** provide: AM scores H1047R, E545K and the ATP pocket as
   uniformly pathogenic; BE3D's screen signal is **drug-directional** — it identifies which functional
   residues drive *this drug's* resistance (helical activation → MEK/EGFR bypass), separating them from
   equally-conserved pocket/catalytic residues that are irrelevant to escape.

## 2–3. Parameter sweep (focused, as specified)
AF primary; structure_radius **4 / 6 / 8** + **r6 with function_for_lfc3d=sum**; nRandom=1000; clustering 6 Å.
8 runs (2 arms × 4 configs), **all RUN_COMPLETED SUCCESS**. Outputs: `real_output/P2_PIK3CA/runs/<arm>_<cfg>/PIK3CA/`.
- **Radius:** union hit count grows modestly with radius (Tram p<0.01: r4=122, r6=124, r8=123); the
  ACT-specific enrichment is stable (R/E 2.9–3.6 across r4–r8). r6 is representative.
- **mean vs sum:** sum inflates the base rate and **destroys** the ACT enrichment (Tram r6_sum ACT
  R/E 1.8, union 0.87) — mean is the better aggregator here; sum over-weights multi-guide neighborhoods.
- **Robustness:** 87 residues persist at pos p<0.01 in ≥3/4 configs (`outreach/PIK3CA/PIK3CA_hotspots.tsv`).

## 4. The NEW question — drug-escape vs activation geometry (answered)
**Under the on-target PI3K inhibitor (pictilisib): NO resistance hotspots anywhere — not in the ATP
pocket, not in the helical domain** (union R/E ~1, ACT n.s., POCKET 0). PI3Kα does not escape
ATP-competitive pan-PI3K inhibition by pocket mutation in this screen (a genuine, informative negative).
**Under the bypass drug (trametinib/MEK): resistance localizes entirely to the helical/regulatory
ACTIVATION hotspots (E542/E545 nSH2-interface cluster + E453/E726), ATP pocket silent.** So the two
axes are cleanly separated: **the resistance geometry here is "activation geometry" (p85/nSH2-interface
disinhibition driving pathway bypass), not "drug-escape geometry" (ATP-pocket).** This is the paper-level
point the pooled hit list does not make: PIK3CA is a *bypass-resistance* node, and BE3D localizes its
escape signal to the disinhibition interface, not the catalytic pocket.

## 5. False positives / false negatives (honest)
- **TP (validated, recovered):** **E542, E545, E726** (helical/regulatory activating hotspots; E453
  present at p<0.05). The dominant 3D peak = E542/E545 nSH2-interface cluster.
- **FALSE NEGATIVES — base-editing-inaccessible (the dominant limiter):**
  - **H1047R / H1047L — the single most common PIK3CA oncogenic mutation — is NOT reachable** (no CBE/ABE
    guide installs it; kinase activation-loop hotspot). Also NOT reachable: **M1043, G1049, H1065, C420,
    E418**. BE3D is structurally blind to the entire kinase-domain hotspot cluster in this library.
  - Consequence: the H1047R kinase-domain arm of the NEW question cannot be tested from this screen —
    the pocket being "silent" is partly true biology (no on-target escape) and partly a reachability gap.
- **LIKELY FALSE POSITIVES:** of 87 robust hotspots, only 3 are validated (TP). Given the 8% base rate,
  the 73 "NOVEL-CAND" residues are **dominated by base-rate noise**, not credible discoveries — flagged
  honestly. 8 are low-pLDDT (disordered ABD/linker), 3 AlphaMissense-benign. Treat the novel set as
  low-confidence; the credible claim is the *concentration* at the helical hotspots, not the long tail.
- **Precision/recall FEEL:** of 87 robust hotspots, ~3 match validated functional sites, the strongest
  3D peak is the correct one (E542/E545), the rest are largely noise; the #1 oncogenic residue (H1047R)
  is a hard false negative for library reasons.

## 6. Competitor comparison
| tool | input | on PIK3CA | advantage | disadvantage |
|---|---|---|---|---|
| **BE3D (this)** | trametinib screen LFC | helical E542/E545 as **MEK-bypass** resistance; pocket silent | drug-directional, separates bypass from pocket; sees screen phenotype | high base rate; misses H1047R (unreachable); union enrichment n.s. |
| **AlphaMissense** | seq/struct | E542K,E545K,H1047R,pocket all uniformly pathogenic | no screen needed; R/E 2.29 p 1e-5 (beats BE3D on recall) | drug-agnostic; can't say *which* residue drives *this* drug's escape |
| **3dhotspots.org / cBioPortal** | patient mutation recurrence | E542/E545/H1047 activation codons | decade-long benchmarking | recurrence-only; blind to allosteric/direction/drug; H1047R would dominate |
| **ProTiler-Mut** | tiling→3D (RRA) | conceptual rival | native tiling model | no drug-directional null; same reachability limit |

## 7. Bottom line
BE3D on the PIK3CA×trametinib screen **recovers the validated helical activating hotspots (E542/E545) as
its strongest 3D signal and correctly localizes resistance to activation geometry rather than the drug
pocket** — a genuinely useful, drug-directional result. But it **does not beat AlphaMissense on functional
recall, has a high (8–14%) base rate, and is blind to the #1 oncogenic hotspot H1047R for base-editing
reachability reasons.** Report as a qualified partial success, not a clean discrimination win.
