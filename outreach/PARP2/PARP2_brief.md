# PARP2 base-editing resistance screen → BE3D 3D map (1-page brief)

**What this is.** The PARP2 arm of the same public base-editing resistance screen (Coelho/Dincer 2024;
olaparib & niraparib, MHH-ES-1, ABE) run through **BE3D** (LFC3D on AlphaFold Q9UGN5). Independent analysis,
benchmarked honestly.

**Headline (honest): this run is under-powered — treat as exploratory only.**
- Unlike PARP1, BE3D's discrimination on PARP2 **fails at the residue level**: neither direction
  significantly concentrates on the catalytic NAD+/inhibitor pocket (depletion R/E ≈ 1.4, p ≈ 0.3;
  only 3 of 9 pocket residues reach significance), precision@10 = 0, and a **no-screen AlphaMissense
  baseline beats BE3D** (R/E 2.9, p 5×10⁻⁵ for the pocket).
- Most likely cause is **power**: the PARPi arm is **ABE-only, 809 guides → ~200 scored residues** over a
  583-aa protein, so the catalytic domain is thinly and unevenly tiled.
- Base rate is high (16–22 % of scored residues called at p<0.001), so the significance labels are not a
  reliable discovery set here.

**Cross-paralog conservation (the useful biological point).** The PARP2 catalytic pocket is the structural
equivalent of PARP1's (H428/Y462/E558 ↔ H862/Y896/E988). In **both paralogs the pocket falls in the
depletion (essential) direction and is never in the resistance direction** — a consistent cross-paralog
signature. BE3D recovers this pocket **strongly in PARP1 (p≈10⁻¹⁰) but only weakly in PARP2 (3/9, NS)**: the
conserved surface is real; BE3D's power to see it is paralog-dependent and coverage-limited here.

**Data note.** The screen numbers PARP2 in a transcript offset by **−13** vs canonical Q9UGN5 (a short-isoform
convention); we corrected to canonical (+13) before mapping to structure (refAA match 14 % → 86 %).

**Olaparib vs niraparib.** Weak maps, Jaccard ≈ 0.33 (more similar than PARP1), both in WGR/CAT — but given
the failed discrimination this concordance is not strongly interpretable.

**Recommendation.** Add a **CBE arm and more guides** for PARP2 before drawing residue-level conclusions;
the current run supports the paralog-conservation statement but not specific PARP2 hotspot claims.

**Files.** `PARP2_hotspots.tsv` (flagged, with under-power caveat), `PARP2_{Olap,Nirap}_G2P.tsv`,
full log in `research/P2_analysis_PARP2.md`.
