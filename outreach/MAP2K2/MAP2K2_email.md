Subject: MEK2 trametinib resistance in 3D — and it mirrors MEK1's escape surface

Hi [name],

Following up on the MEK1 analysis: I ran the MEK2 (MAP2K2) arm of the public Coelho/Dincer 2024 base-editing
screen through BE3D (structure-function tool, Iqbal Lab @ Broad + Liau Lab @ Harvard), trametinib/HT29 arm,
resistance direction. Independent analysis, and I've tried to be honest about the noise.

Candidly: MEK2 is noisier than MEK1. BE3D flags a high fraction of the protein (24–32%), and the N-terminal
disordered tail (1–32) is over-called — I've flagged those as likely false positives. That said, the
enrichment against a validated-resistance set is real (odds ratio 3–4, p<1e-3; it even enriches the bare
structural pocket, p≤0.007), and functional sites are hit 5–9× more than tolerant ones. Use p<0.001 plus a
pLDDT filter and treat it as hypothesis-generating.

The result I think is genuinely interesting is cross-paralog: the MEK1 and MEK2 **trametinib-escape surfaces
are spatially conserved**. All 27 MEK1 allosteric-pocket residues map onto MEK2 with a clean +4 offset, same
amino acids, <0.6 Å apart after superposition (C121→C125, M143→M147, …), and — recovered independently from
the two separate screens — 64% of MEK1's resistance hotspots land on the same aligned position that MEK2
also flags (6/7 pocket hits shared). Both paralogs light the same allosteric wall and the same N-lobe
activating cluster. Practically, that argues a resistance-barrier-raising allosteric strategy on one paralog
should transfer to the other.

Validated MEK2 hits pulled out zero-shot: C125, L119, L122, F133, M147 (pocket); F57/Q60/K61, Y134
(activating).

Files attached: `MAP2K2_hotspots.tsv` (per-residue scores + TP/novel/FP flags) and `MAP2K2_G2P.tsv` for the
interactive 3D view at g2p.broadinstitute.org. If you end up validating any of these — especially whether a
MEK1-derived allosteric strategy holds on MEK2 — that feedback would help us calibrate BE3D.

Best,
[name]

Attachments: MAP2K2_hotspots.tsv, MAP2K2_G2P.tsv, MAP2K2_brief.md
