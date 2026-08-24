# BE3D Run — Decision Log — BCL2 (venetoclax/BH3-groove cross-check vs MCL1 — a HONEST NEGATIVE)

Runner: Queen "BCL2/MYC" (round 2). Date: 2026-08-24. venv py3.14 / beclust3d clone.
Outputs: `real_output/P2_BCL2/runs/{AF_r4_mean,AF_r6_mean,AF_r8_mean,AF_r6_sum}/BCL2/`,
`real_output/P2_BCL2/discrim_bm_report.json`, `outreach/BCL2/`.

## 0. Target, data reality, and the honest question
- Gene / protein: **BCL2** (UniProt **P10415**, 239 aa canonical = clinical numbering; **offset 0** — 282/282
  single-pos missense refAA match to AF; NWGR/BH1 at 143-146, G101/F104 confirmed).
- Structure: **AF-P10415** (primary). Folded blocks (pLDDT>70): **8-32** (BH4/α1), **87-111** (α2/BH3),
  **119-203** (α3-α8 core carrying the BH1/BH2 groove). The **flexible loop domain (FLD, ~33-86)** is
  low-pLDDT/disordered (matches UniProt). Corroborating complex prepared: **6O0K** (venetoclax-BCL2; ligand
  `LBM`) — used only to derive the BE3D-independent groove residue set.
- **DATA REALITY (the headline caveat):** the Coelho ST2 sheet contains **NO venetoclax / navitoclax arm.**
  Its drug arms are HT29 (DebCet/Tram/Pict), H23 (Adag/Sotor), PC9 (Osim/Gefit), MHHES1 (Olap/Nirap) — none
  target BCL2. The only BCL2-relevant readout is **essentiality** via the **Control (no-drug) dropout arms**.
  **But BCL2 is NOT essential in these solid-tumor lines:** nonsense-guide mean z is **positive** in every line
  (HT29 +0.28, H23 +0.51, PC9 +0.39) — knocking BCL2 out causes **no depletion**. So there is **no functional
  window** on BCL2 in this dataset, for either resistance *or* essentiality.
- **Consequence for the brief's question:** we **cannot** test venetoclax resistance (no arm) and there is **no
  essentiality signal** on the BH3 groove to recover. This run is therefore a **null-input / over-call test**:
  given a screen with no real BCL2 signal, does BE3D honestly find little — or manufacture hotspots?

## 1. Data provenance & reformatting
- Converter `scratchpad/convert_bm.py`; screen `scratchpad/bm_tsv/BCL2_HT29ctrl.tsv` (HT29 Control arm, the
  representative line — all Control arms are equivalently null): **5,499 rows** = 568 Missense (233 distinct
  positions, span 1-239), 131 Silent, 24 Nonsense, 37 Splice-donor, 6 Intron, 4,733 No-Mutation. Pooled ABE+CBE.
- QA (cases=Nonsense+Splice, controls=No Mutation+Silent): near-null (KS D≈0.16-0.29 across lines, p 0.03-0.49;
  the nonsense window points the **wrong** way — KO enriches). **ACCEPT for exploratory over-call test only**,
  with the explicit caveat that there is no depletion signal to map. Interpreted direction = NEG (depletion).

## 2. Parameter sweep (4 runs, all EXITCODE 0)
`structure_radius ∈ {4,6,8}` mean + `r6 sum`; clustering_radius 6; nRandom 1000; AF monomer.
Driver `scratchpad/run_bm_sweep.sh`.

## 3. A0 DISCRIMINATION BENCHMARK — `discrim_bm_report.json`
Ground-truth **FUNCTIONAL** = venetoclax/BH3-groove (6O0K ligand contacts ≤4.5 Å ∪ clinical resistance
residues): {100,101,103,104,107,108,111,112,115,133,136,137,143,144,145,146,148,149,152,153,156,198,202,113,129}
(M=25). **TOLERANT** = FLD `range(35,92)`. Baselines: burial (86 res), folded pLDDT>70 (139), AlphaMissense>0.7 (119).

**Baseline (AF, r6, mean, NEG direction), N_neg universe = 125:**

| metric | value |
|---|---|
| base rate p<0.05 | **49.6%** (62/125 sig) — **massive over-call** ("everything is a hotspot") |
| **groove (FUNCTIONAL) enrichment** p<0.05 | **R/E = 0.50, OR = 0.31, p = 0.98** (DEPLETED, not enriched) |
| groove enrichment across sweep | R/E 0.5→0.92 (r6-sum best), **OR always < 1**, never significant |
| discrimination gap (groove/FLD) | ~0.95–1.6× — **meaningless** at a 50% base rate |
| precision@10 (groove) | **0.00** (top-10 depletions are buried-core helices, not the groove) |

**Beats-baseline? NO — it loses to the trivial predictors.** On the SAME protein: "call all buried" gives
groove R/E=**1.2** (OR 1.46), "call all folded" R/E=**1.3** (OR 5.14), "AlphaMissense>0.7" R/E=**1.51** (OR 4.58)
— all weakly *positive* for the groove, whereas **BE3D is *depleted* there (R/E 0.5)**. BE3D adds no
discrimination beyond burial/conservation here; it under-performs them.

**What BE3D actually flags:** a weak folded-vs-disordered bias only (hit folded 52-53% vs disordered 32-41%;
folded R/E 1.08-1.13, p 0.13→0.019 at r8). The strongest robust depletions are **Y180 (−1.20), L181 (−1.12),
G141 (−1.26), D140/V142** — buried α4/α6 core residues (a **generic mild-destabilization dropout** pattern),
**not** the BH3 groove. Only **L137 and E136** (groove-adjacent) reach robustness, weakly.

## 4. Robust hotspots (≥3/4 runs sig p<0.05) — `outreach/BCL2/BCL2_hotspots.tsv`
**59 robust hotspots**, of which only **2 (L137, E136) are BH3-groove contacts** (~3%, below the M/N base
expectation), **43 generic-folded-core (LIKELY-FP)**, **9 low-pLDDT (LIKELY-FP)**, **5 FLD-disordered
(LIKELY-FP)**. → precision on the venetoclax groove is **~0**; essentially the entire robust set is spurious,
driven by diffuse noise on a null screen plus a faint burial bias.

## 5. Cross-target check vs MCL1/S63845 — does BE3D reproduce the BH3-groove result? **NO (and it cannot here).**
- In `P2_run_MCL1.md`, the groove result was **NOT** obtained from plain dropout (which was underpowered and
  *did not* light up the groove) — it came from the **A1331852-sensitized DEPENDENCY arm**, a synthetic-lethal
  co-treatment that **opens an essentiality window** on MCL1 (LFC3D to −3.4 on α3/α4 groove residues
  215/216/220/255/256/258), and from the **S63845 resistance arm**. Both are *drug/dependency* readouts.
- **BCL2 has neither** in Coelho: no venetoclax (resistance) arm, and no BCL2-sensitizing co-treatment to make
  it essential. The plain no-drug Control arm on BCL2 is the direct analogue of MCL1's *underpowered plain
  dropout* — and, exactly as for MCL1's plain dropout, it does **not** find the groove.
- So the honest cross-target conclusion is **not** "BE3D fails on BCL2," but **"the MCL1-style groove result is
  a property of the DEPENDENCY/RESISTANCE arm, and that arm does not exist for BCL2 in this dataset."** The
  method's cross-target consistency claim is **untestable here**; reproducing it would require a venetoclax
  screen (e.g. a CLL/AML line) or a BCL2-synthetic-lethal co-treatment.

## 6. HONEST verdict
- **Negative / null result, reported truthfully.** With no venetoclax arm and BCL2 non-essential in these lines,
  there is **no BCL2 functional signal** in Coelho to recover. BE3D correctly fails to invent a groove
  signal — but it **badly over-calls** (base rate ~50% of scored residues; ~95% of "robust" hotspots are
  likely false positives), it is **anti-enriched** on the true groove (OR<1), and it **loses to burial/folding/
  AlphaMissense** baselines. This is a genuine **false-positive-control failure on a null input** and the main
  BE3D lesson from this target: the p<0.05 randomization threshold is **not** a usable prioritizer when the
  underlying screen carries no real signal — an effect-size / FDR / base-rate gate is needed (see §7).
- **False positives:** effectively the entire robust set (Y180/L181/G141/D140/V142 buried core; FLD; low-pLDDT).
- **False negatives:** the entire venetoclax groove (G101/D103/F104/A113/V156/R107/R129 + NWGR) — because the
  screen contains no venetoclax/dependency signal that would deplete on groove disruption.

## 7. Competitor note
- For BCL2, the **recurrence/AlphaMissense** competitors are strictly better here: AlphaMissense (needs no
  screen) at least weakly concentrates on the groove (R/E 1.51), and 3dhotspots/HotMAPS carry decade-long
  benchmarking. BE3D's advantage (functional-screen resistance/allosteric signal) **requires a relevant screen
  arm**, which this dataset lacks for BCL2 — so on this input BE3D has no edge to offer.

## 8. BE3D issues → BE3D_IMPROVEMENTS.md
(1) High base rate / over-call on null-signal inputs — need an effect-size + FDR gate, not raw p<0.05;
(2) no venetoclax arm in Coelho — a "no relevant assay arm" pre-flight check would have flagged BCL2 as
un-runnable-for-purpose before the sweep. See appended entries.
