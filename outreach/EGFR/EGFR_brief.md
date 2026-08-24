# EGFR TKI-resistance base-editing screen through BE3D — 1-page brief

**What:** We ran the public Coelho/Dincer 2024 base-editing screen (Nat Genet, PMID 39424923; PC9
EGFR-mutant NSCLC) through **BE3D** (structure-function 3D hotspot tool, Iqbal Lab @ Broad + Liau Lab @
Harvard), separately for **gefitinib (1st-gen)** and **osimertinib (3rd-gen)**. Positive z = resistance
(enrichment). This is an independent analysis, not a BE3D-author endorsement.

**Headline (honest):**
- BE3D **re-finds a real gefitinib-resistance compartment zero-shot**: the **αC-helix / exon20 face —
  S768, V769, D770, N771, H773** (exon20 insertions + S768I are validated intrinsic 1st-gen-TKI
  resistance loci). Formal enrichment for the validated set is significant for gefitinib (odds ratio 18.5,
  Fisher p=1.5e-3; hit-rate 86% on functional vs 37% on tolerant residues) and beats a burial baseline.
- For **osimertinib the same test fails** (no enrichment, p=0.49) — because osimertinib's escape residues
  are largely unreachable by base editors (below).
- **The dominant limitation is base-editing chemistry, not BE3D.** ~18/25 canonical EGFR resistance/driver
  residues **cannot be installed by CBE/ABE**: T790M, C797S, L718Q, G724S, L858R, exon19del all need
  transversions/indels the editors don't make. The screen instead made e.g. C797R/G (which *kill* the
  kinase) — so BE3D can't see the clinical resistance alleles. Expect these as false negatives.
- **Caveat on the top signals:** BE3D's strongest positive hits sit in the **disordered C-terminal tail**
  (res >1000, pLDDT<50) — treat as likely false positives, not pocket resistance. Base rate is high
  (~30–40% of scored residues called), so weight formal enrichment, not raw overlap.

**Decision-relevant residues:** validated TP = S768/V769/D770/N771/H773 (gefitinib αC/exon20 face).
Novel-candidate (unvalidated, kinase C-lobe patch) = 806–812. Everything in the C-terminal tail: discount.

**Files:** `EGFR_hotspots.tsv` (residue, AA, LFC3D, pLDDT, TP/novel/FP flag), `EGFR_G2P.tsv` (interactive
3D view input for g2p.broadinstitute.org), per-drug variants `EGFR_{Gefit,Osim}_*`.

**Next step:** the reachable, decision-relevant hypothesis is that the αC/exon20 face (768–773) is a
shared 1st/3rd-gen escape compartment; the classic pocket alleles need saturation/prime-editing screens
(not CBE/ABE) to be assayed.
