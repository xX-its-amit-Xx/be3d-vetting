# MEK1 (MAP2K1) trametinib-resistance 3D map — 1-page brief

**What this is.** An independent re-analysis of the public Coelho/Dincer 2024 (*Nat Genet*) base-editing
screen of MEK1 in HT29, run through **BE3D** (structure-function tool; Iqbal Lab @ Broad + Liau Lab @
Harvard) to map **trametinib** resistance onto the MEK1 structure (AlphaFold Q02750; allosteric pocket
defined from PDB 7JUR, the trametinib-class crystal). Direction analysed = **positive z = resistance**.

**Does BE3D beat chance here? (honest A0 benchmark.)**
- Base rate is high: BE3D calls **13–20%** of residues significant (p<0.001 → p<0.05). At p<0.05 that is
  weak; the credible signal is at **p<0.001**.
- Against a literature-curated, BE3D-independent set of **validated trametinib/allosteric-MEKi resistance
  residues**, BE3D is enriched **~3–4× (odds ratio 6–9, p<1e-4)**; functional-site hit-rate is **5×** the
  tolerant-residue rate; **precision@10 = 40% vs 9% expected**.
- It does **NOT** beat AlphaMissense on raw enrichment (AM flags the whole conserved core), and it does not
  specifically enrich the *bare* structural pocket beyond chance. **BE3D's added value = a resistance
  DIRECTION AM can't give, and a 3× tighter call set.** Verdict: **genuine but modest** discrimination.

**The 3 decision-relevant results.**
1. **A tight allosteric-pocket resistance shell**, distinct from the ATP site: C121, L115, L118, F129,
   M143 (gatekeeper) + αC/β3 rim 112–123. These are the pocket-disrupting residues (C121S, L115P are
   clinically/experimentally validated MEKi resistance). BE3D recovered them **zero-shot from screen signal**.
2. **The N-lobe activating cluster** (helix-A 44–61: F53, Q56, K57 + G128, Y130, E203) — the same set the
   recurrence-based 3dhotspots.org tool reports. BE3D unifies both classes in one map.
3. **A structured distal candidate patch (357/371/372, C-lobe αH/αI, ~18–28 Å from the pocket)** — a
   hypothesis for an "allosteric-of-the-allosteric" escape. **Caveat:** literature assigns this region to
   *non-trametinib* MEKi resistance, so treat as a lead to test, not a claim.

**Explicitly flagged FALSE POSITIVES:** residues 284–287 (disordered αF–αG loop, pLDDT 40–46) — likely
spurious geometry. **False negatives:** P124 and some pocket residues (V211/F209/I141) — base-editing
codon-inaccessible or below threshold, not geometry failures.

**How to use it.** Files: `MAP2K1_hotspots.tsv` (per-residue LFC3D, direction, significance, pocket
distance, pLDDT, AlphaMissense, TP/novel/FP flag), `MAP2K1_G2P.tsv` (load at g2p.broadinstitute.org for
the interactive 3D view). Suggested validation: introduce pocket (C121/L115/M143) vs distal (357/371/372)
substitutions into MEK1 and measure trametinib IC50 shift.

*Independent analysis; not endorsed by the BE3D authors. Honest about uncertainty by design.*
