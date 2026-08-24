# MEK2 (MAP2K2) trametinib-resistance 3D map + MEK1↔MEK2 conservation — 1-page brief

**What this is.** Independent BE3D re-analysis of the public Coelho/Dincer 2024 base-editing screen,
**trametinib** (HT29) arm, MEK2 (AlphaFold P36507). Direction = positive z = resistance. Paired with the
MEK1 analysis to ask whether the drug-escape surface is **structurally conserved across the paralogs**.

**Does BE3D beat chance? (honest A0 benchmark.)**
- Signal is real: enrichment vs a BE3D-independent validated-resistance set is significant (**OR 3–4,
  p<1e-3**), including the *pure structural* allosteric pocket (OR 2.6–3.1, p≤0.007), and functional-site
  hit-rate is **5–9×** the tolerant rate.
- BUT the base rate is **very high (24–32%)** — BE3D flags up to a third of MEK2 — so precision@10 is only
  **2×** base (weaker than MEK1's 4.4×), and the N-terminal 1–32 disordered tail is heavily over-called.
- Does not beat AlphaMissense on raw enrichment. **Verdict: positive but noisier than MEK1 — use p<0.001 +
  a pLDDT>60 filter; hypothesis-generating, not a precise list.**

**The headline result — cross-paralog spatial conservation.**
- All **27/27** MEK1 allosteric-pocket residues map to MEK2 with a clean **+4 offset**, every one the
  **same amino acid**, CA–CA ≤0.6 Å after structural superposition (RMSD 2.14 Å over 392 Cα). C121→**C125**
  (the paralog of the C121S resistance cysteine); M143→M147 gatekeeper; etc.
- Independently, **64% of MEK1's robust resistance hotspots fall on the same aligned position that is also
  a MEK2 hotspot** (6/7 pocket hits shared). Both paralogs light the same allosteric wall (β3–αC +
  catalytic loop + activation segment) **and** the N-lobe activating cluster.
- **Conclusion: the MEK1 and MEK2 trametinib-escape surfaces are spatially conserved** — recovered
  independently from two separate screens. Practically: a co-targeting or next-gen allosteric strategy that
  raises the resistance barrier on MEK1's pocket should transfer to MEK2's, and vice versa.

**Validated MEK2 hits (zero-shot):** C125, L119, L122, F133, M147 (pocket wall); F57/Q60/K61, Y134
(activating). **Flagged false positives:** N-terminal 1–32 (disordered, low pLDDT).

**Files:** `MAP2K2_hotspots.tsv` (per-residue LFC3D/flags), `MAP2K2_G2P.tsv` (g2p.broadinstitute.org).

*Independent analysis; not endorsed by the BE3D authors.*
