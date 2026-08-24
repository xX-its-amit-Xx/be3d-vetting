# Outreach email — BRCA1 (clinical-genetics / VUS lab)

**Subject:** BRCA1 base-editing screen → 3D functional hotspots (spatial PS3-type evidence; honest caveats)

Hi [Name],

I ran the public Cuella-Martin et al. 2021 (*Cell*) BRCA1 base-editing dropout screen (CBE, ± olaparib and cisplatin, MCF10A) through **BE3D**, a structure-function tool from the Iqbal Lab (Broad) + Liau Lab (Harvard) that aggregates base-editor tiling-screen signal over 3D neighborhoods to surface functional hotspots. I benchmarked the output hard against your kind of gold standard — Findlay et al. 2018 SGE and ClinVar — and want to share both the wins and the limits, because for VUS work the limits matter.

**What BE3D got right (a genuine, orthogonal win).** From screen signal alone it flagged the **RING Zn-coordinating cluster C61/C64** as a top hotspot. Both are loss-of-function in Findlay SGE and pathogenic in ClinVar — and notably **AlphaMissense mis-calls C64G/C64R as *likely benign*** (0.09–0.11), a real predictor blind spot that the functional-screen signal rescues. BE3D also concentrates its hits on the folded, functionally-constrained **RING (2–109) and BRCT (1642–1855)** domains at ~2.4–3.8× the rate of the disordered central region — i.e., it is not "flagging everything" (genome-wide significant-residue rate ~7%).

**Where it is weak (please don't use it alone).** Its *domain-level* concentration is largely recapitulated by a trivial "is this residue in a folded domain / high-pLDDT" baseline, and its *residue-level* resolution inside the domains is modest (enrichment for Findlay-LOF residues ~1.4×, p≈0.06; per-residue rank concordance with Findlay function scores ≈ 0). Its single strongest cluster, the **BRCT-linker residues 1801–1804**, is scored **tolerant by Findlay SGE** — a likely false positive (or a base-editor-specific effect the SGE didn't capture). A central cluster at **214–222** is a reproducible multi-guide dropout signal but sits in low-confidence disordered structure — a hypothesis, not a call.

**How I'd use it for VUS.** As **supporting spatial functional evidence (toward ACMG PS3), in concert with** Findlay SGE and AlphaMissense — strongest exactly where those disagree or are silent (e.g., the C64 case). Not as a standalone classifier.

Attached: the robust hotspot table (with TP / novel / likely-FP flags), the G2P interactive-view input, and a 1-page brief. If you validate any of these (esp. 1801–1804 vs the RING Zn site), I'd value the results back — it helps us tune BE3D's null model and thresholds.

Best,
[Your name]

*Independent analysis; not an endorsement by the BE3D authors. Screen: Cuella-Martin 2021, PMID 33592168. Gold standard: Findlay 2018, PMID 30209399.*
