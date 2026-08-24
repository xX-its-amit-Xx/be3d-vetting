# BCL2 — BE3D structure-function brief (1 page) — an HONEST NEGATIVE

**What was run.** The public Coelho/Dincer 2024 (*Nat Genet*) base-editing dataset was queried for **BCL2**
(5,562 guides, ABE+CBE) and run through **BE3D** (LFC3D + randomization null + 3D clustering). Structure:
AlphaFold P10415 (clinical/canonical numbering; venetoclax-resistance residues G101/D103/F104 etc. map at
offset 0). Parameter sweep: structure radius 4/6/8 + r6-sum, nRandom 1000.

**The critical data caveat (up front).** This Coelho dataset contains **no venetoclax or navitoclax arm** — its
drugs (BRAF/MEK/PI3K, KRAS-G12C, EGFR, PARP inhibitors) do not target BCL2. The only BCL2-relevant readout is
**essentiality** from the no-drug arm, and **BCL2 is not essential in these solid-tumor lines** (knockout
guides do not deplete; nonsense mean z is *positive*). **So this dataset carries no BCL2 functional signal —
neither venetoclax resistance nor essentiality — to recover.** This run is therefore a null-input / over-call
control, not a venetoclax-resistance map.

**Result — a true negative, and an over-call warning.**
- **No BH3-groove signal.** Against an independent venetoclax/BH3-groove ground-truth set (6O0K drug contacts +
  clinical resistance residues), BE3D is **not enriched — it is depleted** (observed/expected ≈ 0.5, odds
  ratio 0.31, p≈0.98). **precision@10 on the groove = 0.**
- **It over-calls.** ~**50%** of scored residues clear the p<0.05 threshold; **~95% of the 59 "robust" hotspots
  are likely false positives** (buried α4/α6 core residues Y180/L181/G141 — a generic mild-destabilization
  dropout pattern — plus disordered/low-pLDDT residues). Only 2 groove-adjacent residues (L137, E136) surface,
  weakly.
- **It loses to trivial baselines.** On the same protein, "call all buried" / "call all folded" / "AlphaMissense
  >0.7" all weakly *do* concentrate on the groove (R/E 1.2–1.5); BE3D does not.

**Does it reproduce our MCL1 BH3-groove result? No — and it cannot here.** For MCL1, the groove lit up only via
a **drug/dependency arm** (S63845 resistance + an A1331852 synthetic-lethal sensitization that makes MCL1
essential); MCL1's plain no-drug dropout did **not** find the groove. BCL2's no-drug Control arm is the exact
analogue of that underpowered plain-dropout — and, consistently, it does not find the groove. The MCL1-style
result is a property of the resistance/dependency arm, which **does not exist for BCL2 in this dataset**. The
cross-target consistency check is **untestable here**; it would need a venetoclax screen (CLL/AML line) or a
BCL2-sensitizing co-treatment.

**Bottom line.** BE3D honestly did not manufacture a groove signal that isn't in the data — but on this
null input it over-calls and under-performs burial/AlphaMissense, so **do not use these BCL2 hotspots**. To
map venetoclax resistance on BCL2, BE3D needs a venetoclax base-editing arm.

**Files.** `BCL2_hotspots.tsv` (flagged, ~all LIKELY-FP), `BCL2_G2P.tsv`, this brief. Full benchmark:
`research/P2_analysis_BCL2.md`.
