# MCL1 base-editing screen × BE3D — 1-page brief

**What / how.** The public MCL1 tiling base-editing screen (Hanna et al., Cell 2021; CBE in MELJUSO,
3 arms: plain dropout, A1331852-sensitized dependency, and A1331852+S63845 resistance) projected onto
MCL1 structure (AF Q07820; BH3-groove crystals 3PK1/3MK8; S63845-class inhibitor pocket 6QB4) with
BE3D — per-variant screen LFC → 3D spatial clustering vs a randomization null. Independent analysis.

**Headline.** The dependency signal and the S63845-resistance signal land on the **same BH3 groove**
→ S63845/HVN is a **true orthosteric BH3-mimetic**. Separately, BE3D resolves a second, **off-groove
resistance compartment** (stability/degron) that the source paper never spatially resolved.

## Discrimination (does BE3D beat chance? — the primary test)
| arm | FUNCTIONAL hit-rate (groove+pocket) | TOLERANT hit-rate (distal surface) | gap | groove enrichment |
|---|---|---|---|---|
| combo / S63845 resistance | 100% (11/11) | 12.5% | +87.5% (8×) | OR ∞, p=7.9e-4 |
| A1331852 dependency | 92.3% (12/13) | 30.0% | +62.3% (3.1×) | OR ∞, p=5.6e-4 |
| plain dropout (QA-failed) | 50% | 62.5% | −12.5% | ns → noise |
- Base rate is high (~40% of scored residues flagged), so the case rests on this gap, not raw overlap.
- Beats trivial baselines: burial doesn't predict this surface groove; AlphaMissense calls 86% of the
  domain pathogenic and can't tell an escape residue from a generically-essential one.

## The actionable residues
- **Orthosteric escape watchlist (de-risk in next-gen mimetics):** R215, V216, V220, F254, D256, V258,
  T259, S269 + the α3 core G217–Q221 — robust across all 4 structures; strongest single edit **T212I**
  (+2.44). *Test:* DMS of these codons under compound selection; co-crystallize T212I mutants.
  (The buried P2-floor F270/T266/M231 are NOT on this list — CBE can't reach them; false negatives.)
- **Safe anchor points:** pocket-only residues 246/271/274/290/294 line the drug but carry NO
  dependency/resistance signal → contacts here are unlikely to create escape liabilities (selectivity).
- **Second resistance axis = stability (degrader opportunity):** off-groove α2 (176–180), C-term
  (325/326), N-term PEST (78/85/95); >8 Å from drug; likely act via MCL1 half-life (FBW7 phosphodegron,
  Wertz Nature 2011). A degrader would bypass groove-escape mutations. *Test:* MCL1 half-life for
  E85K/S178L/E325K; do they blunt the mimetic but stay degrader-sensitive?

## Honesty — false positives, false negatives, competitor context
- **Novel-candidate caveat:** the stability compartment sits in low-pLDDT disordered regions (weak 3D
  support) — hypotheses, not clusters.
- **Likely false positive:** distal surface patch 295–303, from the underpowered dropout arm — treat as noise.
- **False negatives (coverage, not scoring):** 19/36 groove/pocket residues were never CBE-editable —
  incl. the buried P2 floor **F270, M231, V274** and the anchor **R263** (near-missed). AlphaMissense
  covers these; a base editor cannot reach them. Complementary methods.
- **Competitor context:** MCL1 is amplification-driven (1q21; Beroukhim Nature 2010), NOT point-mutated —
  3dhotspots.org returns `[]`, cancerhotspots = 0 residues, not a Bailey-2018 driver. Recurrence-based 3D
  tools produce nothing for MCL1; a functional screen is the only route to a 3D map here.

**Files:** `MCL1_hotspots.tsv` (TP/novel/FP flags), `MCL1_S63845resistance_G2P.tsv` &
`MCL1_A1331852dependency_G2P.tsv` (interactive 3D). **Key refs:** Kotschy 2016 (PMID 27760111,
S63845/R263), Czabotar 2007 (PMID 17389404, P1–P4 groove), Wertz 2011 (PMID 21358673, FBW7 degron),
Cheng 2023 (PMID 37733863, AlphaMissense), Beroukhim 2010 (PMID 20164920), Gao 2017 (PMID 28115009).
