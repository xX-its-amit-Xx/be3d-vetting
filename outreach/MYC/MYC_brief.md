# MYC — BE3D structure-function brief (1 page)

**What was run.** The public Coelho/Dincer 2024 (*Nat Genet*) base-editing screen of **MYC** (2,498 guides,
ABE+CBE) was run through **BE3D** (LFC3D + randomization null + 3D clustering; Iqbal Lab/Broad + Liau
Lab/Harvard). Coelho has no MYC-targeted drug, so the functional readout is **essentiality** from the no-drug
(Control) arm — MYC knockouts deplete strongly (nonsense mean z = −2.67; QA D=0.354, p=7e-7). Direction
interpreted = depletion/essential. Structure: AlphaFold P01106. **Numbering note: the AF model and screen use
the 454-aa MYC1 frame (= classic 439-numbering + 15).**

**The question (a deliberate hard test).** MYC is ~80% intrinsically disordered; only the C-terminal
**bHLH-LZ** (AF 365-453) folds and does the work (DNA binding + MAX dimerization). Does BE3D **concentrate**
signal on that small folded domain, or **over-call** across the disordered bulk (a BRCA1-style stress test)?

**Result — it concentrates, and passes at the domain level.**
- **Enrichment vs an independent (1NKP DNA+MAX contact) ground-truth set:** observed overlap ≈ 2× expected;
  **odds ratio 5.8–9.9, p from 8e-6 down to 3e-7** across the parameter sweep.
- **Discrimination gap ≈ 4×:** BE3D flags the folded bHLH-LZ at **65–75%** vs the disordered IDR at **16–20%**.
  It does **not** diffusely light up the IDR.
- **The strongest essentiality signals land exactly on function:** DNA-binding basic arginines **R379/R381/R382**,
  the MAX leucine-zipper **Y417/I418/L419/V421**, and **MYC Box II** core **C148/M149**. precision@10 = 0.6–0.7.
- **Beyond a trivial "is it folded"/AlphaMissense predictor:** on the folded contacts, "is it folded" (R/E 3.4)
  and AlphaMissense (R/E 2.1) match or beat BE3D — but BE3D additionally recovers the **low-pLDDT functional
  boxes MBI (60-65) and all of MBII (145-158)**, real cofactor/degradation sites those baselines miss because
  they are disordered. That MBI/MBII recovery is BE3D's distinctive, screen-driven contribution.

**Honest caveats (please read).**
- Resolution is **domain-level, not residue-level**: within the compact bHLH-LZ, BE3D does **not** rank the
  actual DNA/MAX contacts above their folded neighbours (within-domain enrichment ≈ chance). Use it to say
  "the functional domain is engaged," not to pick single contact residues.
- ~26 scattered **IDR false positives** remain (res 37-47, 92-95, 120-127, 317-319; base rate ~32%). Treat
  p<0.05 as a screen; prioritize by magnitude/rank + the enrichment test, not by the flag alone.
- Direction is essentiality (dropout), not a MYC-drug resistance map (no MYC drug exists in this screen).

**Files.** `MYC_hotspots.tsv` (74 robust residues, flagged TP/functional/FP), `MYC_G2P.tsv` (interactive
3D view input), this brief. Full methods + the A0 benchmark: `research/P2_analysis_MYC.md`.

**Track record.** BE3D has zero-shot recovered known functional hotspots before (e.g. KBTBD4's
medulloblastoma R379–R390 Kelch site from screen signal alone); here it independently re-found MYC's
DNA-binding basic region, MAX leucine-zipper, and MYC Box II. This is an independent analysis, not a claim
endorsed by the BE3D authors. If you validate/refute any of these residues, that feedback helps tune the
null model and thresholds.
