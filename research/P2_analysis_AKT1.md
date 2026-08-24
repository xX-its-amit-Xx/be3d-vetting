# BE3D Run — Analysis: AKT1 drug-resistance base-editing screen (3D)

> Runner: Queen "PI3K/AKT" (round 2). Real public base-editing screen (Coelho/Dincer 2024) →
> BE3D LFC3D + 3D clustering. **A0 discrimination benchmark leads. This is an HONEST NEGATIVE result.**

## 0. Target, data, and the drug-arm reality
- **Gene / protein:** AKT1 (UniProt **P31749**, 480 aa). Domains: PH 5–108, kinase 150–408,
  hydrophobic-motif/regulatory C-tail 409–480.
- **Structure:** AlphaFold `AF-P31749-F1-model_v6` (WT, canonical numbering). **Primary.**
- **Source:** Coelho/Dincer 2024 `ST2 BE z-scores`; **4,036 AKT1 guide rows** (1,071 missense guides →
  438 distinct missense positions, span 1–480). Score = per-drug average z-score, positive = resistance.
  Converter `scratchpad/convert_pi3k.py`; TSVs `scratchpad/pi3k_tsv/AKT1/`. refAA match 6602/6688 (98.7%).
- **HONEST DRUG CAVEAT:** **no AKT inhibitor (capivasertib / MK2206 / ipatasertib) is in this dataset.**
  The "ATP-competitive vs allosteric AKT-inhibitor escape" framing of the NEW question **cannot be tested
  from this screen** — there is no AKT-directed drug. The closest pathway-relevant arms are **pictilisib**
  (upstream PI3K inhibitor; AKT1 activation could bypass) and **gefitinib** (EGFR bypass). Both run here.
- **Signal scan (missense guides Z>3):** pictilisib 3, gefitinib 5, others 0–5; **max Z ≈ 4.5** — much
  weaker than PIK3CA/KRAS. **Even the canonical activating hotspot E17 has max Z 0.38** (see FN below).

## 1. BE-QA
Cases = Nonsense+Splice-donor, controls = No Mutation+Silent (pictilisib arm). Screen accepted for
processing, but note the missense signal is weak throughout (few guides Z>3). Positive direction only.

## A0. DISCRIMINATION BENCHMARK (primary) — **FAIL**
Ground-truth (independent, cited — `groundtruth_pi3k.json`):
- **FUNCTIONAL-ACTIVATING** (M=5): E17 (PH), L52, Q79, W80 (PH–KD autoinhibitory interface), D323
  (Carpten 2007 PMID 17611497; Parikh 2012 PMID 23134728).
- **FUNCTIONAL-POCKET/CATALYTIC** (M=13): ATP cleft + catalytic — K179, E228/A230 hinge, E234, E278,
  M281, T291, D292 (DFG), D274/K276 (HRD), T308 (activation loop), S473 (HM) (Alessi 1996 PMID 8978681).
- **FUNCTIONAL-UNION** M=18 (base rate M/N = **0.037**). **TOLERANT** n=117 (AlphaMissense-benign − functional).

**Results (AlphaFold, r6, mean, n1000, positive direction):**

| arm | thr | n_sig | **base rate** | UNION R/E (p) | ACT R/E (p) | POCKET R/E | **gap hitF / hitT** |
|---|---|---|---|---|---|---|---|
| Pict | p<0.05 | 137 | **0.285** | 0.39 (0.98) | 1.40 (0.44) | 0.0 | 0.11 / **0.26** |
| Pict | p<0.01 | 115 | **0.240** | 0.46 (0.95) | 1.67 (0.35) | 0.0 | 0.11 / **0.24** |
| Pict | p<0.001 | 88 | 0.183 | 0.30 (0.98) | 1.09 (0.64) | 0.0 | 0.06 / **0.21** |
| Gefit | p<0.01 | 91 | 0.190 | 0.29 (0.98) | 1.05 (0.65) | 0.0 | 0.06 / 0.12 |

**Precision@10 (Pict, r6):** top-10 = `138, 160, 319–324, 379, 468` — 1/10 functional (323). No pocket.

**Beats a trivial baseline?** **NO — fails badly.** AlphaMissense-high: **R/E 2.16, OR 23.7, p 5.9e-6**
(highly significant). Burial: R/E 1.36 (n.s.). BE3D union R/E < 1 (anti-enriched), p ≈ 0.95–0.98.

### A0 verdict (AKT1) — HONEST: **NEGATIVE / NON-DISCRIMINATING**
1. **Base rate is enormous** (24–29% of residues flagged at p<0.01–0.05) → textbook "everything is a
   hotspot." 104 of 480 residues (22%) are robust hotspots.
2. **Union enrichment < 1** (depleted), never significant. BE3D hotspots do **not** concentrate on
   functional residues.
3. **Discrimination gap is INVERTED:** BE3D flags TOLERANT residues at a **higher** rate than FUNCTIONAL
   (hitT 0.24 > hitF 0.11 for pictilisib). This is the defining failure mode — no discrimination.
4. Fails to beat AlphaMissense (which works fine here, OR 23.7). BE3D adds nothing on AKT1.
**This is a real negative and is reported as one.** The underlying cause is the input: the AKT1 screen
under these (non-AKT) drugs carries little resistance signal, and BE3D's neighborhood aggregation over a
near-null screen just paints the compact 480-aa fold broadly.

## 2–3. Parameter sweep
AF primary; radius 4/6/8 + r6 sum; nRandom 1000. 8 runs (Pict, Gefit × 4), all RUN_COMPLETED SUCCESS
(`real_output/P2_AKT1/runs/`). Base rate rises steeply with radius (Pict p<0.01: r4=20%, r6=24%, r8=31%)
— the small fold saturates. No config rescues discrimination (union R/E stays <1). sum inflates further.

## 4. The NEW question — cannot be answered as posed
There is **no AKT inhibitor** in the dataset, so "ATP-competitive vs allosteric AKT-inhibitor escape"
is untestable here. What we *can* say: under upstream PI3K (pictilisib) / EGFR (gefitinib) inhibition,
**AKT1 shows no localized, discriminating resistance signal** — neither the PH-domain activating cluster
(E17/L52/W80) nor the ATP pocket lights up above the 22–29% background. The pocket R/E = 0.0 throughout.

## 5. False positives / false negatives
- **TP (validated, weakly recovered):** L52, D323 (PH–KD interface) appear among robust hotspots, but at
  no better rate than background — not a credible recovery given the inverted gap.
- **FALSE NEGATIVES — base-editing reachability (critical):**
  - **The canonical activating hotspot E17K is NOT installed.** The only edit at codon 17 is **E17G**
    (ABE A>G: GAG→GGG), a **non-activating** substitution; oncogenic **E17K** (G>A) is not made by the
    editors used. So AKT1's #1 hotspot is effectively absent from the screen (E17 max Z 0.38).
  - **T308** (activation-loop phosphosite) is **NOT reachable**. Several pocket residues are reachable
    (K179, A230, M281, T291) but carry no resistance signal (expected — no AKT inhibitor to escape).
- **LIKELY FALSE POSITIVES:** 104 robust hotspots, ~0 credible TP, 17 low-pLDDT, 9 AM-benign. Essentially
  the entire hotspot list is base-rate noise. Do not circulate AKT1 "novel candidates."

## 6. Competitor comparison
| tool | on AKT1 | verdict |
|---|---|---|
| **BE3D (this)** | 22% of residues flagged, gap inverted, union R/E<1 | **fails to discriminate on this input** |
| **AlphaMissense** | E17, PH–KD interface, pocket pathogenic; R/E 2.16 OR 23.7 p 6e-6 | works; the right tool absent a screen |
| **3dhotspots.org** | E17K activation codon | recurrence works; BE3D can't beat it here |

## 7. Bottom line
**AKT1 is a clean negative for BE3D in this dataset.** No AKT inhibitor is present; the pathway-bypass
arms give a near-null AKT1 screen; BE3D's neighborhood aggregation over that null yields a 22–29% base
rate, sub-chance functional enrichment, and an inverted discrimination gap — it flags tolerant residues
as readily as functional ones. Compounding this, the signature AKT1 activating mutation **E17K is not
even base-editable** in this library (only non-activating E17G). Honest recommendation: **do not use
BE3D output for AKT1 from these arms**; the tool is not at fault so much as the input carries no signal —
but the result must be reported as a non-discriminating failure, per the A0 honesty rule.
