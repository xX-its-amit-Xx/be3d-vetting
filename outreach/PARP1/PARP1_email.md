Subject: PARP1 base-editing resistance screen → a 3D structure-function readout (independent analysis)

Hi [name],

I ran the public PARP1 base-editing screen from Coelho & Dincer 2024 (olaparib and niraparib arms) through
BE3D, a structure-function tool from the Iqbal Lab (Broad) and Liau Lab (Harvard) that maps per-residue
screen log-fold-change onto 3D structure (LFC3D on AlphaFold P09874) and tests it against a randomization
null. This is an independent analysis; I benchmarked it against curated ground truth before drawing anything.

The clean result: BE3D re-found the PARP1 catalytic machinery zero-shot. In the depletion direction the
whole NAD+/inhibitor pocket (H862, Y896, E988, S904, Y907…) plus the autoinhibitory HD helix light up as one
tight cluster (hypergeometric p ≈ 1e-10, all 12 pocket residues, persistent across the parameter sweep) — and
this beats both a burial baseline and AlphaMissense. It's a strong positive control that the spatial
aggregation is doing real work (it has re-found sites like KBTBD4's medulloblastoma Kelch cluster the same way).

Being honest about the resistance direction, which is what you'd actually care about: it's broad. BE3D flags
10-28% of the protein as "resistance," and that set isn't enriched above chance (or above AlphaMissense) for
known PARPi-resistance residues. What is useful is the ranked top of that list — it surfaces real ones at 4/4
sweep persistence: M1 (start-loss LOF), the ZnF DNA-binding cluster F44/D45/K119/S120 (Pettitt 2018), WGR
R591, BRCT L390/S391 — all consistent with loss-of-trapping resistance. Treat these as prioritized
hypotheses, not a validated set.

One thing that might interest you: olaparib and niraparib give only moderately overlapping 3D resistance maps
(Jaccard ≈ 0.22), with niraparib's roughly 2x broader across the ZnF/WGR DNA-binding surface — a difference in
breadth rather than a clean drug-specific pocket.

Attached: the flagged hotspot table (TP/novel/FP), the G2P files for interactive 3D viewing, and a one-page
brief. If you have PARP1 variant/validation data, I'd genuinely value it back — it would help calibrate the
resistance-direction thresholds, which are the weak spot here.

Best,
[name]
