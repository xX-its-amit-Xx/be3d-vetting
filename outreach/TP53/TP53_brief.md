# BE3D on TP53 — 1-page brief

**What.** BE3D (Iqbal/Liau labs; `broadinstitute/BE3D`) maps base-editing tiling-screen LFC onto a 3-D
structure, aggregates signal over spatial neighborhoods (LFC3D), scores hotspots against a randomization
null, and clusters them. Independent analysis — not endorsed by the BE3D authors.

**Input.** Two public human TP53 base-editing screens (MaveDB `urn:mavedb:00001245-a-1` ABE8e, `-a-2` CBE;
"activity-based selection", etoposide arm). Structure: AlphaFold P04637 monomer (numbering verified).
**Direction:** p53 LOF = **negative LFC / depletion** (nonsense median -3.0 vs synonymous ~0); analysis
uses the negative direction.

**Discrimination benchmark (meta of both editors, radius 6 Å, nRandom 1000).**
- Base rate: **36%** of scored DBD residues flagged (high — TP53 DBD is mutation-hypersensitive; stated up front).
- FUNCTIONAL enrichment (17 curated DNA-contact/Zn/structural residues): **16/17 recovered, R/E = 2.6,
  Fisher p = 4e-7.**
- Discrimination gap: **94% functional hit vs 14% tolerant** (ClinVar-/AlphaMissense-benign), +81%.
- **84% of cancerhotspots.org recurrent p53 hotspots recovered zero-shot.**
- Precision@20 = 30% (≈5× the 5.2% chance rate).

**Beats a trivial baseline? Partly.** BE3D (OR 33) beats a **burial** predictor (OR 4.9) — key point: the
DNA-contact hotspots (R248, R280, A276) are solvent-*exposed* in the apo monomer, so burial misses them.
**But AlphaMissense (OR 55, tolerant-hit 0%) out-discriminates BE3D** on p53's conserved core. BE3D's edge
is orthogonal: real phenotype, LOF direction, and **editor accessibility** (ABE8e recovers 78% of recurrent
hotspots, CBE 49%, both together 84%).

**Hotspot ledger (117 robust residues).** 16 TP-structural + 21 TP-recurrent (cBioPortal) + 67 novel
candidates + 13 likely-FP. Strongest signals: R249, P250, R280, L3 loop (M243/G244/M246/N247), D281, R248,
P278, R282, Zn-C242, R196, H193, G245.

**New hypotheses (not claims).** Contiguous LOF surfaces beyond textbook hotspots: the **L3 DNA-binding
loop** (M243/G244/M246/N247/P250) and the **loop-sheet-helix** (G279/R283/T284). G244/P250 are real but
under-catalogued cancer positions.

**Honest failure modes.** False negative: **H179** (Zn) missed (neighbors recover the site). False
positives: **S261, N263, A129, D228** (exposed, ClinVar/AM-benign). Monomer RSA conflates true surface
with DNA-contact exposure — don't over-trust it.

**Reproduce.** BE3D is public (Colab + `be3d_local.py`). Reformatted screens + G2P file + full annotated
hotspot TSV attached. Interactive 3-D: Genomics-2-Proteins portal (g2p.broadinstitute.org).
