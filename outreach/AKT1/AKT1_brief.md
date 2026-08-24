# AKT1 — BE3D 3D resistance analysis (1-page brief) — HONEST NEGATIVE

**What this is.** We ran the public Coelho/Dincer 2024 base-editing screen of **AKT1** through **BE3D**
on the AlphaFold model of P31749 (480 aa). Independent analysis. Positive z = resistance.

**Headline: this is a negative result, reported honestly.** We do **not** recommend using the AKT1
hotspot list for decisions. Two reasons:

1. **No AKT inhibitor is in the dataset.** Capivasertib / MK2206 / ipatasertib are absent. The only
   pathway-relevant arms are upstream (pictilisib/PI3K, gefitinib/EGFR), under which the AKT1 screen
   carries very little resistance signal (max missense z ≈ 4.5; the canonical hotspot E17 max z 0.38).
   The requested "ATP-competitive vs allosteric AKT-inhibitor escape" question is **untestable** here.

2. **BE3D does not discriminate on this input.** Base rate is ~24–29% of residues (104 of 480 are
   "robust hotspots"). Functional-residue enrichment is **below chance** (union R/E < 1, p ≈ 0.97), and
   the **discrimination gap is inverted** — BE3D flags AlphaMissense-tolerant residues at a *higher* rate
   (0.24) than validated functional residues (0.11). AlphaMissense on the same protein works fine
   (R/E 2.16, OR 23.7, p 6e-6); BE3D adds nothing.

**Base-editing reachability (a real limiter):** the signature activating mutation **E17K is NOT
installed** — the only codon-17 edit is **E17G** (ABE A>G), a non-activating substitution; oncogenic
E17K (G>A) is not made by these editors. The activation-loop phosphosite **T308 is also unreachable.**

**What we would say to an AKT1 group:** absent an AKT-directed inhibitor screen, this dataset cannot
localize AKT1 resistance, and BE3D's neighborhood aggregation over a near-null screen just paints the
compact fold broadly. If you have a **capivasertib or MK2206 base-editing screen**, re-running BE3D on
that would be worthwhile — the method needs a real, AKT-directed resistance phenotype to work with.

**Files (provided for transparency, not for action):** `AKT1_hotspots.tsv` (all flagged low-confidence),
`AKT1_G2P.tsv`. Full benchmark: `research/P2_analysis_AKT1.md`.
