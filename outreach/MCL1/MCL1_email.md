Subject: MCL1 base-editing screen through BE3D — the BH3-groove escape rim + a second, off-groove resistance axis

Hi [Name],

I ran the public MCL1 base-editing tiling screen (Hanna et al., Cell 2021, MELJUSO CBE) through
BE3D, a structure-function 3D-hotspot tool (Iqbal Lab, Broad; Liau Lab, Harvard) that projects
per-variant screen log-fold-changes onto a protein structure and tests for spatial clustering
against a randomization null. This is an independent analysis, not endorsed by the BE3D authors —
sharing it because it speaks directly to BH3-mimetic escape and to the degrader question.

Why trust the signal here first: BE3D has recovered known functional sites zero-shot before (e.g.
KBTBD4's medulloblastoma R379–R390 Kelch site from screen signal alone). On MCL1 it independently
concentrated the S63845-resistance and the (A1331852-sensitized) dependency signals on the SAME
BH3 groove — functional hit-rate 92–100% vs 12–30% on matched distal-surface residues (3–8× gap,
Fisher p<0.001, robust across AF, 3PK1, 3MK8 and 6QB4). Honest caveat: BE3D flags ~40% of the
residues it can score, so I lean on that enrichment/discrimination gap, not raw overlap — and on
MCL1 the gap is real.

Three decision-relevant results:

1. S63845/HVN is a true orthosteric BH3-mimetic. The resistance contacts coincide with the groove
   MCL1 uses to bind pro-apoptotic BH3 partners — not a separate allosteric pocket. The robust,
   structure-independent groove ESCAPE residues to de-risk in next-gen compounds: R215, V216, V220,
   F254, D256, V258, T259, S269, plus the α3-helix core G217–Q221 (+ the strongest single edit,
   T212I). Design pharmacophores that don't hinge on these side chains, or that tolerate their
   substitution. (BE3D did NOT flag the buried P2-floor F270/T266 — see blind spots below.)

2. A second, orthogonal resistance axis = MCL1 stability. Off-groove clusters — α2 (R176/Q177/S178/
   E180), the C-terminus (E325/G326), and the N-terminal PEST/phosphodegron (R78/E85/R95) — are all
   >8 Å from the drug and most likely act via MCL1 abundance/half-life (the FBW7 phosphodegron,
   Wertz et al. Nature 2011). Implication: a DEGRADER would bypass the groove-escape mutations that
   defeat orthosteric mimetics — these regions are candidate handles. I flag these as hypotheses,
   not validated: they sit in low-pLDDT disordered regions, so treat as sequence-level calls.

3. Honest blind spots. BE3D missed the buried P2-pocket floor (F270 side, M231, V274) — CBE simply
   can't make those codon changes — and it over-calls a distal surface patch (295–303) in the
   underpowered plain-dropout arm. AlphaMissense covers the buried residues BE3D can't edit; the two
   are complementary. Recurrence tools (3dhotspots.org, cBioPortal) find nothing on MCL1 because it's
   amplification-driven, not point-mutated — a functional screen is the only way to get a 3D map here.

Next step: the attached tables list every hotspot with TP/novel/FP flags; the G2P TSVs drop into the
interactive 3D viewer. If useful, the concrete experiment is a deep-mutational scan of the 9 groove
codons under compound selection, and MCL1 half-life measurement for E85K/S178L/E325K knock-ins
(do they blunt the mimetic but stay degrader-sensitive?). If you validate or refute any of these,
we'd value the results back — they help us tune BE3D's thresholds and null model.

Happy to run any of your in-house MCL1 screens the same way.

Best,
[Your name]

Attachments: MCL1_hotspots.tsv · MCL1_S63845resistance_G2P.tsv · MCL1_A1331852dependency_G2P.tsv · MCL1_brief.md
