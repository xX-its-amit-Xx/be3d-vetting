Subject: MEK1 trametinib-resistance, mapped in 3D from your base-editing screen

Hi [name],

I ran the public MEK1 arm of the Coelho/Dincer 2024 base-editing resistance screen through BE3D — a
structure-function tool from the Iqbal Lab (Broad) and Liau Lab (Harvard) that aggregates screen log-fold-
changes over 3D neighborhoods and tests them against a randomization null. I focused on the **trametinib**
(HT29) arm and the resistance (enrichment) direction. Sharing in case it's useful; it's an independent
analysis, not a sales pitch, and I've tried to be candid about what does and doesn't hold up.

Honest headline first: BE3D flags a fairly high fraction of MEK1 (~13% at p<0.001), so I benchmarked it
hard. Against a literature-curated set of validated trametinib/allosteric-MEKi resistance residues (curated
without BE3D), it's enriched ~3–4× (odds ratio 6–9, p<1e-4), hits functional sites 5× more often than
tolerant ones, and precision@10 is 40% vs 9% expected. It does *not* beat AlphaMissense on raw enrichment
(AM flags the whole conserved core) — BE3D's edge is a resistance *direction* AM can't give, and a tighter
call set. So: genuine but modest.

What it found, zero-shot from screen signal:
- A tight **allosteric-pocket resistance shell**, distinct from the ATP site — C121, L115, L118, F129, M143
  (gatekeeper) + the αC/β3 rim. C121S and L115P are the classic validated MEKi-resistance mutations.
- The **N-lobe activating cluster** (helix-A F53/Q56/K57, G128, Y130, E203) — the recurrence-based
  3dhotspots set — in the same map.
- A structured **distal candidate** (357/371/372, C-lobe, ~20 Å from the pocket) as an "allosteric-of-the-
  allosteric" hypothesis. Caveat: that region is linked to *non-trametinib* MEKi resistance in the
  literature, so I'd treat it as a lead to test, not a result.
- I've flagged the likely false positives too (disordered loop 284–287, low pLDDT).

If useful: the attached `MAP2K1_hotspots.tsv` has per-residue scores/flags, and `MAP2K1_G2P.tsv` loads
directly at g2p.broadinstitute.org for an interactive 3D view. A clean test would be pocket (C121/M143) vs
distal (357/371) substitutions and a trametinib IC50 readout. If you validate (or refute) any of these,
that feedback would genuinely help us tune BE3D's thresholds and null model.

Best,
[name]

Attachments: MAP2K1_hotspots.tsv, MAP2K1_G2P.tsv, MAP2K1_brief.md
