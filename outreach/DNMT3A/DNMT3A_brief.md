# DNMT3A base-editor scan through BE3D — 1-page brief

**What.** Independent 3D structure-function re-analysis of the Lue/Liau DNMT3A base-editor scanning
screen (Nat Chem Biol 2023, DOI 10.1038/s41589-022-01167-4) using BE3D (Iqbal Lab, Broad; Liau Lab,
Harvard). Screen = CpG-methylation→citrine reporter FACS in MOLM13, cytosine base editor (CBE).
Structure = PDB **4U7T** (DNMT3A–DNMT3L–H3; Guo et al. 2015 Nature), catalytic domain chain A
(residues 474–912). Loss-of-function (LOF) = loss of DNMT3A methylation. UniProt Q9Y6K1 numbering.
*Independent analysis; not endorsed by the BE3D or Lue authors.*

**Did it actually discriminate? (benchmarked first, before any claim).** Yes.
- Base rate: 117/314 scored residues significant at p<0.05 (**37%** — high, so we use rank + enrichment,
  not the raw flag).
- Enrichment of LOF hotspots for an independently-curated FUNCTIONAL set (catalytic pocket + RD +
  FF interfaces + clinical R882/R736/R729, M=42): **R/E = 2.2–2.5, odds ratio ~12, p ≈ 1e-11**.
- Discrimination gap: **83%** hit-rate on functional residues vs **24%** on tolerant residues (~3.4×).
- Precision@10 = **0.60** (4.5× over chance); @20 = 0.55.
- Beats a burial baseline outright (burial has *no* functional signal, p=0.69).

**The biology (ranked surfaces).** Catalytic/SAM pocket dominates overall (odds ratio 38 — expected for
an activity screen). Among protein–protein surfaces: **RD homotetramer interface ≫ FF (DNMT3A–DNMT3L)
≫ (ADD autoinhibitory, H3-reader, DNA face — all depleted).** The **RD tetramerization interface is the
dominant editable vulnerability**, robust in 8/12 parameter runs.

**The actionable result — a base-editable R882 proxy.** A clinically **unannotated** RD-interface patch —
**T671, M674, V675, R676, H873, D876, V877, S878** — is rock-solid (T671/D876/V877/S878 p<0.001 in
12/12 runs). It abuts R882 (the R882–D876 salt-bridge network) but is distinct from it. **Hypothesis:** these
phenocopy the R882H dominant-negative tetramerization defect without touching R882 — a clean allele
series to dissect tetramerization vs catalysis, and a candidate **tetramerization-disrupting drug epitope**
(inhibit DNMT3A by blocking oligomerization, avoiding the pan-methyltransferase SAM pocket).
Validated positive controls BE3D re-found: catalytic pocket (D686/V665/E664/C666/D641…) and the FF
interface **R736** (2nd-most-frequent DNMT3A cancer residue) and **R729** (abolishes processive catalysis).

**Honest limitations.**
- **R882 = false negative.** R882 is the #1 somatic DNMT3A hotspot in AML (~50–60% of DNMT3A-mutant
  AML; 53% AML vs 10.6% CHIP, Venugopal 2022) and is dominant-negative via tetramer blockade
  (Russler-Germain 2014, PMID 24656771). **BE3D misses it** because R882H is a G→A change a CBE
  cannot install, and no R882 guide exists in the screen. For R882, **somatic recurrence tools
  (cBioPortal / 3dhotspots.org) win** — BE3D's value is the editable proxy + interface resolution.
- High base rate → trust the rank, not the binary flag. ~62 weaker p<0.05 residues are lower-confidence.
- 4U7T covers 474–912 only (real PWWP hits G293/E294 are out of interface scope). DNA face is
  literature-curated (soft). Per-residue AlphaMissense and a curated ClinVar-benign set are open items.

**Competitor snapshot.** Somatic tools (cBioPortal/3dhotspots/COSMIC): win at R882, decade-benchmarked,
but recurrence ≠ function and can't attribute the interface. AlphaMissense: no screen needed, but calls the
whole catalytic domain uniformly intolerant and can't rank by editability. ProTiler-Mut (Cell Systems 2026):
same 3D-aggregation idea, screen-agnostic; on a CBE screen it inherits the same R882 blind spot.
**BE3D adds:** functional screen signal + controls, per-surface interface resolution, and named base-editable
handles (the RD proxy patch) that recurrence-based tools cannot see.

**Suggested experiment.** Test S878F / V877 / D876 / H873 / T671 (± M674/V675/R676) by CBE/prime editing
or site-directed mutagenesis in the DNMT3A methyltransferase + tetramerization assay (analytical SEC /
cooperative DNA binding). Prediction: hypomorphic/dominant-negative tetramer defect phenocopying
R882H, orthogonal to the confounded R882 codon.

**Files.** `DNMT3A_hotspots.tsv` (ranked, flagged) · `DNMT3A_G2P_LOF_neg.tsv` (interactive 3D input).
