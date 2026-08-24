# PIK3CA (p110α) — BE3D 3D resistance analysis (1-page brief)

**What this is.** We ran the public Coelho/Dincer 2024 base-editing drug-resistance screen (Nat Genet,
DOI 10.1038/s41588-024-01948-8) of **PIK3CA** through **BE3D** (structure-function tool, Iqbal Lab/Broad
+ Liau Lab/Harvard) on the AlphaFold model of P42336 (1068 aa). Independent analysis, not endorsed by
the BE3D or paper authors. Positive z-score = resistance (proliferation under drug).

**Key result (honest, drug-directional).** The dataset has **no alpelisib**; the on-target PI3K
inhibitor present is **pictilisib** (pan-PI3K, ATP-competitive) and it produces **essentially no
resistance signal** — PI3Kα does not escape ATP-competitive inhibition by pocket mutation in this screen.
The strong signal is under **bypass drugs** (trametinib/MEK, gefitinib/EGFR): PIK3CA **activating**
mutations confer resistance by reactivating the pathway. BE3D's strongest 3D peak under trametinib is the
**E542/E545 helical-domain hotspot cluster** (the p85-nSH2 disinhibition interface) — recovered zero-shot
from screen signal alone. Resistance localizes to **activation geometry (helical/regulatory), not the ATP
drug pocket** (pocket enrichment R/E = 0).

**Benchmark (be skeptical — we are).**
- Base rate is high: BE3D flags ~8–14% of residues → weak per-site evidence. Say so.
- Functional-set enrichment (union of activating+pocket+catalytic): R/E ~1.7, **not significant** (p=0.15).
- **But activating-hotspot enrichment is real: R/E 2.9–3.6, p=0.008–0.02**, and the #1 3D peak = E545.
- **Does NOT beat AlphaMissense** for generic functional recovery (AM R/E 2.29, p 1e-5). BE3D's value-add
  is *drug-directionality* — which functional residues drive *this* drug's resistance — not raw recall.
- **False negative:** the most common PIK3CA oncogenic mutation, **H1047R, is NOT base-editing-reachable**
  in this library (no CBE/ABE guide installs it); BE3D is blind to the kinase-domain hotspot cluster.

**Validated hotspots recovered (TP):** E542, E545 (helical activating; E545 dominant), E726.
**Do not over-read** the ~73 unvalidated "novel candidates" — at an 8% base rate they are mostly noise.

**Files:** `PIK3CA_hotspots.tsv` (residue, domain, LFC3D, persistence, pLDDT, AlphaMissense, editable
substitutions, TP/novel/FP flag), `PIK3CA_G2P.tsv` (Genomics-2-Proteins portal 3D view input,
g2p.broadinstitute.org). Full log: `research/P2_analysis_PIK3CA.md`.

**Suggested next step.** If you care about *trametinib/MEK-inhibitor bypass*: the helical E542/E545
interface is the escape locus to watch (ctDNA surveillance / co-targeting). If you care about *on-target
PI3K-inhibitor resistance*: this screen says pocket mutations are not the route — look elsewhere
(feedback/adaptive, not primary p110α coding). H1047R must be assayed by a non-base-editing method.
