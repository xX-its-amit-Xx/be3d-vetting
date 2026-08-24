# PARP1 base-editing resistance screen → BE3D 3D map (1-page brief)

**What this is.** We ran the public PARP1 base-editing drug-resistance screen (Coelho/Dincer 2024,
*Nat Genet*; olaparib & niraparib arms, MHH-ES-1, ABE) through **BE3D** (Iqbal Lab, Broad + Liau Lab,
Harvard), which aggregates per-residue screen log-fold-change over 3D structural neighborhoods (LFC3D on
AlphaFold P09874) and tests significance against a randomization null. **This is an independent analysis,
not endorsed by the BE3D authors.** We benchmarked it honestly before interpreting.

**Headline (honest).**
- **BE3D re-discovered the PARP1 catalytic machinery zero-shot** in the *depletion* direction: the entire
  NAD+/inhibitor **catalytic pocket (H862, Y896, E988, S904, Y907, G863…) and the autoinhibitory HD helix
  (678–701, 767–779)** form one tight, robust hotspot (hypergeometric p ≈ 10⁻¹⁰, all 12 pocket residues,
  3–4/4 sweep persistence). This **beats both burial and AlphaMissense baselines** — a genuine
  structure-function win and a strong positive control that the method works.
- **The resistance (enrichment) direction is broad and low-specificity.** BE3D flags 10–28 % of the protein
  as "resistance," and this set is **not enriched above chance / above AlphaMissense** for curated PARPi-
  resistance residues (olaparib at chance; niraparib R/E≈1.8, p≈0.01). Use the **ranked** top positives as
  hypotheses, not the significance labels as a discovery list.

**Resistance residues BE3D does surface (top-ranked, 4/4 persistent) — validated:** **M1** (start-loss LOF),
the **ZnF1/2 DNA-binding cluster F44/D45, K119/S120, G161/V164** (Pettitt 2018), **WGR W589/R591**, **BRCT
L390/S391**. Mechanism: loss of DNA engagement/expression → loss of PARP-trapping → olaparib/niraparib
resistance.

**Olaparib vs niraparib (3D).** Their positive resistance maps overlap only moderately (**Jaccard ≈ 0.22–0.24**).
Both center on the **DNA-binding ZnF/WGR module**, but **niraparib's map is ~2× broader and more DNA-binding-
centric**. This is a difference in *breadth*, not a clean drug-specific pocket signature. A minority of
resistance hits fall in the **HD allosteric domain** (release-of-autoinhibition escape), but HD is dominated
by the depletion signal.

**Caveats.** ABE-only arm (CBE escape alleles unreachable); positive z-scores assigned per-guide to all
edited residues; ~70 positive residues sit in low-pLDDT disordered regions (flagged likely-FP); the catalytic
pocket appears as depletion (trapping toxicity of catalytic point-mutants), not resistance.

**Files.** `PARP1_hotspots.tsv` (flagged TP/novel/FP), `PARP1_{Olap,Nirap}_G2P.tsv` (G2P 3D viewer),
full decision log in `research/P2_analysis_PARP1.md`.
