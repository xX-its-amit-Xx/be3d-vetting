Subject: BE3D 3-D map of the public TP53 base-editing screens — recovers 16/17 core p53 residues, plus an honest AlphaMissense comparison

Hi [Name],

I ran the two public TP53 base-editing tiling screens (MaveDB `urn:mavedb:00001245-a-1` ABE8e and
`-a-2` CBE, the etoposide activity-based-selection sets) through **BE3D**, a structure-function tool
from the Iqbal Lab (Broad) + Liau Lab (Harvard) that maps per-guide LFC onto 3-D neighborhoods and
scores spatial hotspots against a randomization null. This is an independent analysis, not endorsed by
the BE3D authors. I thought the p53 result was clean enough to be worth your eye.

**Does it actually discriminate, or just flag everything?** The honest answer first: on TP53 the base
rate is high — BE3D calls ~36% of the scored DNA-binding domain a LOF hotspot, because the DBD really is
mutation-hypersensitive. So I benchmarked it properly. Against an independently curated set of 17
functional residues (DNA-contact, Zn-binding, structural), the meta run (both editors) recovers **16/17**
(enrichment R/E = 2.6, Fisher p = 4e-7), flags **94% of functional vs 14% of ClinVar-/AlphaMissense-benign
residues** (discrimination gap +81%), and recovers **84% of cancerhotspots.org recurrent p53 hotspots
zero-shot** (R175, R248, R273, R282, G245, R249, plus R196, H193, P278, D281, H214…). BE3D has done this
kind of zero-shot recovery before — e.g. re-finding KBTBD4's medulloblastoma R379-R390 Kelch site from
screen signal alone.

**Two caveats I want to be straight about.** (1) **AlphaMissense actually beats BE3D** at picking p53's
functional residues (it needs no screen). BE3D's added value here is orthogonal: it reads the actual
etoposide-selection *phenotype* (LOF = depletion, direction documented), and it shows **editor
accessibility matters** — ABE8e alone recovers 78% of recurrent hotspots, CBE only 49%; you need both.
(2) The one canonical residue it misses is **H179** (Zn site) — sub-threshold in both editors, though its
Zn-pocket neighbors C176/C238/C242 are all strongly hit, so the site is recovered as a cluster.

**What's potentially new for you:** beyond the textbook hotspots, BE3D lights up contiguous LOF surfaces —
the **L3 DNA-binding loop (M243/G244/M246/N247/P250)** and the **loop-sheet-helix (G279/R283/T284)** — as
strong as the canonical sites. Several (G244, P250) are real but under-catalogued cancer positions; these
are hypotheses, not claims. A handful of exposed residues (S261, N263, A129) are likely false positives
(ClinVar/AM-benign) — flagged as such in the table.

If useful: the full hotspot table (with TP/novel/FP flags, cBioPortal counts, AlphaMissense scores), the
reformatted screens, and a G2P interactive-3D file are attached; BE3D is public
(github.com/broadinstitute/BE3D, runs in Colab). If you validate or refute any of the novel L3/LSH
candidates, that feedback would directly help us tune BE3D's thresholds and null model.

Best,
[Your name]

Attachments: TP53_hotspots.tsv · TP53_G2P_LFC3D_neg.tsv · TP53_ABE8e/CBE_screen_BE3D_input.tsv · TP53_brief.md
