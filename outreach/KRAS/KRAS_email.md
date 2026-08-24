Subject: A per-drug 3D resistance map for KRAS(G12C) — adagrasib vs sotorasib (independent BE3D analysis)

Hi [Name],

I ran the public KRAS(G12C) base-editing resistance screen (Coelho, Dincer et al. 2024, Nat Genet)
through BE3D, a structure-function tool from the Iqbal Lab (Broad) + Liau Lab (Harvard) that projects
screen log-fold-changes onto 3D structure and finds spatial clusters of resistance. This is an
independent analysis (not endorsed by the BE3D authors). I think two things in it are useful to a
G12C program.

Why trust it: BE3D has re-found known functional hotspots from screen signal alone before (e.g.
KBTBD4's medulloblastoma R379–R390 Kelch site). On KRAS, benchmarked against the validated G12Ci
resistance residues independently of the tool, the **adagrasib** arm is strong: it recovers 8 of 9
known resistance residues zero-shot (S65, R68, D69, M72, D92, H95, Y96, Q99), with a hypergeometric
enrichment of ~5.6-fold over chance (odds ratio 57, p ≈ 1×10⁻⁶) and a 24-fold higher hit-rate on
functional vs mutation-tolerant residues. Top-10 strongest signals are all pocket-lining.

The two results:
1. **Sotorasib and adagrasib have different 3D escape maps.** Adagrasib resistance is a broad
   orthosteric shell across the switch-II pocket (60–72, 89–103; His95/Tyr96/Gln99 groove) — dense,
   so single pocket substitutions are efficient escapes. Sotorasib's in-pocket footprint is sparse by
   comparison (they share only 60,65,66,67,68; Jaccard ≈ 0.12). Design implication: a sotorasib-class
   covalent binder that also reaches the His95/Q99 exit (as adagrasib does) should raise the genetic
   barrier to resistance.
2. **A sotorasib-specific distal candidate at α5/SAK residues 147–156** (~13–27 Å from the drug, well-
   structured, pLDDT 97–99), which the source paper doesn't flag. This is a **hypothesis, not a
   validated site** — honestly, our benchmark shows sotorasib's signal does NOT significantly enrich
   the canonical pocket (it sits off it), so treat the α5 cluster as a testable idea about a distal /
   conformational escape route, not a claim. Prior anchor: F156L (Feng 2022) sits in this patch.

Caveats up front: positive/resistance direction only (enrichment screen); the HVR signal at 178–180
is low-pLDDT and we flag it as a likely artifact; per-guide causality within multi-edit guides is
unresolved (mitigated by aggregation over 162 positions).

How to use it: attached is a per-residue table with TP / novel-candidate / likely-FP flags and a
G2P-formatted file for interactive 3D viewing. Concrete next step to test the α5 hypothesis: introduce
T148/R149/V152/D154-region substitutions into KRAS(G12C) and measure sotorasib (vs adagrasib) IC50
shift + GTP-loading; a Boltz co-fold of KRAS(G12C)+sotorasib with these mutants is a fast in-silico
pre-screen. You can run BE3D on your own screens (Colab/GitHub — happy to point you at it).

If any of this is useful, we'd value your validation results back — confirmation or refutation of the
α5 patch would directly help us tune BE3D's thresholds and null model.

Best,
[Your name]

Attachments: KRAS_hotspots.tsv, KRAS_sotorasib_g2p.tsv, KRAS_adagrasib_g2p.tsv, KRAS_brief.md
