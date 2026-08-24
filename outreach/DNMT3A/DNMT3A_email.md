Subject: DNMT3A base-editor scan → a ranked, editable interface map (incl. an R882 proxy patch)

Hi [Name],

I re-analyzed the public Lue/Liau DNMT3A base-editor scanning screen (Nat Chem Biol 2023) through
BE3D, a structure-function tool from the Iqbal Lab (Broad) + Liau Lab (Harvard) that projects screen
signal onto 3D structure and tests it against a randomization null. This is an independent analysis,
not endorsed by the BE3D or Lue authors — but the result is directly relevant to your DNMT3A/AML work
and I wanted to share it.

Why trust the calls: I benchmarked discrimination before anything else. Against a structure/literature
FUNCTIONAL set curated independently of BE3D, its loss-of-function hotspots are enriched ~2.2–2.5×
(odds ratio ~12, p≈1e-11), hit 83% of functional residues vs 24% of tolerant ones, and reach
precision@10 = 0.60 (4.5× over chance). It beats a "just call buried residues" baseline outright (that
baseline has no functional signal here, p=0.69). Honest caveat: BE3D flags ~37% of scored residues at
p<0.05, so I rely on the RANK and per-surface enrichment, not the raw flag. BE3D has recovered known
functional sites zero-shot before (e.g. KBTBD4's medulloblastoma R379–R390 Kelch site), and here it
independently re-found the catalytic/SAM pocket and the R736/R729 DNMT3L (FF) interface — both
documented functional/cancer residues.

Three decision-relevant results:
1. Among protein–protein surfaces, the **DNMT3A homotetramer RD–RD interface is the dominant
   editable vulnerability** (not the FF or the autoinhibitory ADD contact, which is depleted).
2. A clinically **unannotated, base-editable RD patch — T671, M674, V675, R676, H873, D876, V877,
   S878** — sits next to R882 and is rock-solid (p<0.001 in 12/12 parameter runs). Hypothesis: it
   phenocopies the R882H tetramerization defect *without* touching R882 — a clean proxy allele series
   and a candidate tetramerization-disrupting drug epitope.
3. **Honest blind spot: BE3D misses R882 itself** — R882H is a G→A change a cytosine base editor
   cannot make, and no R882 guide exists in the screen. For R882, somatic recurrence tools
   (cBioPortal/3dhotspots) win. BE3D's value is the editable proxy + interface resolution, not R882.

Concrete next step: test S878F / V877 / D876 / H873 / T671 (and M674/V675/R676) in your
methyltransferase + tetramerization assay (SEC / cooperative DNA binding) as R882-proxy alleles.
Attached: the ranked hotspot table (with TP / novel-candidate / false-positive-flag / false-negative
columns) and the G2P table for interactive 3D viewing. Happy to share the Colab/GitHub so you can
re-run it, and if you validate (or refute) any of these we'd genuinely value the feedback — it helps us
tune BE3D's thresholds and null model.

Best,
[Your name]

Attachments: DNMT3A_hotspots.tsv · DNMT3A_G2P_LOF_neg.tsv · DNMT3A_brief.md
