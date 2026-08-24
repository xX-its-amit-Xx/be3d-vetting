# BRAF RAF-inhibitor-resistance base-editing screen through BE3D — 1-page brief

**What:** We ran the public Coelho/Dincer 2024 base-editing screen (Nat Genet, PMID 39424923; HT29
BRAF-V600E colorectal) through **BE3D** (structure-function 3D hotspot tool, Iqbal Lab @ Broad + Liau Lab
@ Harvard). Primary arm = **DebCet (dabrafenib + cetuximab)**, the RAF-inhibitor arm; contrast arms Tram
(trametinib/MEK) and Pict (pictilisib/PI3K). Positive z = resistance. Independent analysis.

**Headline (honest):**
- BE3D **re-finds the dimer-interface acquired-resistance residue L505 zero-shot**: L505 (with neighbor
  R506) is a robust positive hotspot **across all three drug arms**. **L505H is a validated acquired
  vemurafenib-resistance mutation** (PMID 25515853); the screen installs the reachable L505F (LFC +2.02).
  This directly supports the **new question — resistance clusters at the DIMER INTERFACE / αC-helix, not
  the ATP pocket** — consistent with dimerization-driven RAF-inhibitor resistance.
- **The formal A0 discrimination test does NOT clear chance** (enrichment R/E≈2, Fisher p=0.23; base rate
  ~39%; precision@10=0; does not beat a burial baseline). Reason: only **3/15** validated functional
  residues are base-editing-reachable, so significance is unattainable — this is a **hypothesis-generating**
  result, not a proven discrimination win.
- **Base-editing inaccessibility is severe:** **V600E (the driver) and R509 (dimer interface) have ZERO
  guide coverage** (both need transversions); gatekeeper T529 only gets loss-of-function T529A. So the
  canonical BRAF hotspots are false negatives here — an intrinsic BE-screen limit, not a BE3D fault.
- **Caveat on top ranks:** the highest |LFC3D| residues are N-terminal CR1/CR2 domain positions
  (364/365/403/404) — likely passengers/false positives. Weight the L505 finding, discount those.

**Decision-relevant residues:** validated TP = **L505, R506** (dimer/αC face, all arms). Novel-candidate
(unvalidated αC/β3 patch) = 492–508, 521–533, 543–550. Discount: N-terminal CR domains, 607–609 (pLDDT<50).

**Files:** `BRAF_hotspots.tsv` (residue, AA, LFC3D, pLDDT, TP/novel/FP flag), `BRAF_G2P.tsv` (interactive
3D view), per-arm variants `BRAF_{DebCet,Tram,Pict}_*`.

**Next step:** validate whether αC/dimer-face substitutions (L505, and the 507/508 R-spine neighbours)
shift dabrafenib IC50 in BRAF-V600E cells; the ATP-pocket V600/T529 alleles require saturation editing
(not CBE/ABE) to assay.
