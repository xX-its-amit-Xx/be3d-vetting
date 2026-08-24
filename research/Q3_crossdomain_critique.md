# BE3D — Queen 3 Vetting: Cross-Domain Analogies, Critique & Workflow Adoption

**Scope:** BE3D (`beclust3d` v1.0.0, Broad Institute; MIT) interprets base-editor (BE) tiling-mutagenesis
screens by mapping per-residue log-fold-change (LFC) onto 3D protein structure, spatially smoothing it
(LFC3D), building a randomization null, clustering significant residues into hotspots, and
meta-aggregating across screens. Ground truth: `BE3D_BRIEF.md`, `README.md`, `README_FUNCTIONS.md`.

**Abstracted algorithmic pattern (the key to analogies):**
> Noisy point-wise measurements indexed along a 1D sequence → project onto a 3D spatial embedding →
> kernel/neighborhood smoothing to denoise → permutation/randomization null for significance →
> spatial clustering to call hotspots → meta-aggregate across repeated noisy experiments.

This is, almost exactly, **"spatial kernel smoothing + hotspot/cluster detection + permutation
significance testing on a 3D point cloud, with multi-experiment consensus."** That pattern is mature and
well-studied in at least six other fields. Below, each analogy names the mature technique, gives the
mapping, and states a **concrete improvement BE3D could borrow.**

---

## 1. CROSS-DOMAIN ANALOGIES

### 1.1 Neuroscience / genomics — the *closest* analogues

**(a) Cluster-based permutation tests (Maris & Oostenveld 2007, EEG/MEG/fMRI).**
This is BE3D's method almost line-for-line. The neuro pipeline: compute a per-sensor/voxel statistic
(analogous to per-residue LFC), threshold it, group supra-threshold neighbors into spatial clusters
(analogous to BE3D's agglomerative clustering of significant residues), take a cluster-level summary
statistic, and calibrate it against a **permutation null** built by shuffling condition labels
(analogous to `randomize_data`/`randomize_sequence`). The neuro field solved *exactly* BE3D's
multiple-testing-over-space problem and controls family-wise error at the cluster level.
- **What BE3D can learn:**
  1. **Cluster-mass / cluster-level significance, not just per-residue p<0.05 then cluster.** BE3D
     currently thresholds residues (`znorm_score`) and *then* clusters geometrically. It never assigns a
     p-value to a *cluster*. Borrow the Maris-Oostenveld idea: compute a cluster statistic (sum of LFC3D
     over cluster members) and compare it to the distribution of the *largest* cluster statistic under
     the null. This gives each hotspot a calibrated significance and controls FWER across space.
  2. **Threshold-Free Cluster Enhancement (TFCE; Smith & Nichols 2009)** eliminates the arbitrary
     cluster-forming threshold entirely by integrating cluster support over all thresholds. BE3D's twin
     hard 6 Å parameters (smoothing radius *and* clustering radius) are exactly the kind of arbitrary
     knob TFCE was invented to remove.

**(b) ChIP-seq peak calling — MACS (`macs2/3 callpeak`).**
MACS finds enriched genomic regions by comparing local read pileup to a background. Its signature trick:
a **dynamic local λ** — instead of a single genome-wide Poisson background rate, MACS takes the *max* of
background estimates over several nested window sizes (peak, ~1 kb `slocal`, ~10 kb `llocal`) to absorb
local biases (chromatin, copy-number). Peaks get q-values (Benjamini-Hochberg FDR).
- **What BE3D can learn:** BE3D's null is a *global* label shuffle. Borrow MACS's **local/adaptive
  background** — a residue in a mutation-tolerant, flexible loop should be judged against a *local*
  neighborhood null, not the whole-protein null. A residue-local expected LFC3D (from its structural
  neighborhood under randomization) would reduce false hotspots in intrinsically noisy regions and
  false negatives in quiet, well-ordered cores. And **report q-values**, as MACS does.

**(c) Hi-C domain / TAD calling.** Genomic contact maps are smoothed and segmented into
topologically-associated domains — another "smooth a noisy signal over a spatial neighborhood, then
segment into blocks" problem. Reinforces that block/segment detection with an explicit contact/proximity
kernel is the norm; BE3D's single-linkage merge is a comparatively blunt segmentation.

### 1.2 Geospatial / GIS — the *canonical vocabulary* for this exact task

**(a) Getis-Ord Gi\* hot-spot analysis (ArcGIS, PySAL, GeoDa).** Gi\* is THE standard local hotspot
statistic: for each location it computes a z-score comparing the **local weighted sum** of a variable
(the feature + its neighbors) to the global sum, flagging where high values cluster more than chance.
BE3D's LFC3D — "sum/mean of LFC over residues within radius" — is essentially an *un-normalized,
un-standardized Gi\* local sum*.
- **What BE3D can learn:** Adopt the **Gi\* formulation directly**: it gives a principled z-score and a
  distance-decayed spatial-weights matrix. Critically, Gi\* uses **inverse-distance or Gaussian weights**
  rather than a hard cutoff — see 1.2(c).

**(b) Moran's I / Local Moran's I (LISA) — spatial autocorrelation.** Before you smooth, you should test
*whether* the signal is spatially autocorrelated at all, and at what length scale. Moran's I quantifies
this globally; Local Moran's I flags residues in high-high vs high-low (outlier) neighborhoods.
- **What BE3D can learn:** A **spatial-autocorrelation pre-check** (global Moran's I on LFC over the
  structure) would validate BE3D's core premise per target and *empirically choose the smoothing radius*
  (via a correlogram / semivariogram range) instead of a hard-coded 6 Å for every protein.

**(c) Kriging / semivariogram (geostatistics).** Kriging models the spatial covariance structure
(semivariogram) and produces a best-linear-unbiased *interpolation with uncertainty* at every point.
- **What BE3D can learn:** Replace the **hard 6 Å neighbor set with a Gaussian/exponential kernel** whose
  length scale is fit from the data's own semivariogram range. This removes the brittle binary
  "in/out at 6 Å" decision (a residue at 6.01 Å contributes nothing; at 5.99 Å contributes fully), gives
  smoothly-weighted neighbors, and yields **per-residue uncertainty** — something BE3D currently lacks.

### 1.3 Healthcare / epidemiology — variable-window scanning & study meta-analysis

**(a) Kulldorff spatial scan statistic / SaTScan.** The gold standard for disease-cluster detection.
It scans a **circular window of *varying* radius and center** over the map, computes a likelihood-ratio
for "inside vs outside" at every window, takes the maximum, and calibrates it with Monte-Carlo
permutations — **explicitly correcting for the multiple testing over all window sizes/locations.**
- **What BE3D can learn:** BE3D fixes the neighborhood at 6 Å. Borrow the **variable-radius scan**: scan
  spatial windows of multiple radii over the structure, score each, and take the most-likely cluster with
  a properly multiplicity-corrected p-value. This directly fixes BE3D's single-scale blind spot (it can
  miss a diffuse 12 Å functional patch or over-merge two adjacent 4 Å ones). Caveat from the literature
  ("Spatial scan statistics can be dangerous", 2021): circular windows bias cluster shape — so pair it
  with the flexibly-shaped scan (Tango & Takahashi) for irregular hotspots.

**(b) GWAS multiple-testing & meta-analysis.** GWAS tests millions of SNPs and lives or dies by
multiplicity control (genome-wide significance, Benjamini-Hochberg FDR) and by **inverse-variance /
random-effects meta-analysis** (METAL, GWAMA) across cohorts with heterogeneity statistics (Cochran's Q,
I²).
- **What BE3D can learn:** (i) BE3D tests every residue at p<0.05/0.01/0.001 with **no FDR/q-value across
  residues** — a clear gap. Add Benjamini-Hochberg / Storey q-values. (ii) BE3D's meta step just **SUMs**
  LFC3D across screens. Replace with **inverse-variance weighting** (down-weight noisy screens, informed by
  BE-QA effect sizes) and report **between-screen heterogeneity (I²)** so users see when screens disagree.

### 1.4 Finance — signal denoising, alpha combination, FDR over many bets

**(a) Multi-signal alpha combination.** Combining many noisy predictive signals into one score is the
core of quant equity. Naive summing (BE3D's meta-SUM) is dominated by **information-coefficient / inverse-
variance weighting** and shrinkage combiners.
- **What BE3D can learn:** Weight each screen's contribution by its reliability (BE-QA gives a natural
  quality weight — a screen whose knockout-vs-neutral separation is weak should count less).

**(b) Empirical-Bayes / James-Stein shrinkage.** Estimating many means (per-residue effects) from noisy
data is the textbook shrinkage problem; empirical-Bayes (Efron; Stephens' `ashr` adaptive shrinkage)
"borrows strength" across residues, pulling unreliable estimates toward the global/local mean.
- **What BE3D can learn:** Apply **empirical-Bayes shrinkage to per-residue LFC before 3D smoothing**, and
  to the multi-screen meta-estimate. This is a more principled denoiser than uniform neighborhood
  averaging and naturally handles residues with few guides.

**(c) False-discovery control across many bets & anomaly/hotspot detection.** Quant research formalized
"testing many strategies" (Harvey-Liu-Zhu; deflated Sharpe; Benjamini-Hochberg on backtests) and regime/
anomaly detection on time series.
- **What BE3D can learn:** Same message as GWAS — **q-values and explicit FDR reporting** for the residue-
  and cluster-level calls, plus a "how many hotspots would I expect by chance?" summary (BE3D has the null
  machinery to produce this cheaply).

### 1.5 Astronomy / image processing — source detection & matched filtering

**(a) Source detection (SExtractor, DAOPHOT, `photutils`).** Astronomers detect stars/galaxies by
smoothing an image, estimating a **local background + noise (RMS)**, thresholding at *n·σ above local
background*, deblending merged sources, and cataloging — a direct parallel to hotspot calling.
- **What BE3D can learn:** The **local-background-and-noise + n·σ threshold** discipline (again pointing at
  a local null), and **deblending**: BE3D's single-linkage clustering has no way to split one large merged
  blob into two genuine adjacent hotspots. Astronomy's multi-threshold deblending is a ready recipe.

**(b) Matched filtering.** Optimal detection of a known-shape signal in noise = convolve with a template
matched to the expected signal profile and the noise covariance.
- **What BE3D can learn:** If functional hotspots have a characteristic spatial *size*, a **matched
  (e.g. Gaussian) spatial filter** is provably the optimal linear detector — more sensitive than a flat
  6 Å boxcar sum.

**(c) Sigma-clipping.** Iteratively reject outliers beyond n·σ before estimating background statistics.
- **What BE3D can learn:** Robustify the null and the neighborhood aggregation against a single extreme-LFC
  guide artifact contaminating a whole neighborhood (BE3D's `mean` aggregation is not outlier-robust;
  offering a sigma-clipped or median/Huber option would help).

### 1.6 Art / design / visualization — communicating spatial intensity

Heatmaps, **choropleths**, and kernel-density (hotspot) maps are the mature grammar for showing "where is
the signal hot." Hard-won lessons: (i) use **perceptually-uniform, colorblind-safe diverging colormaps**
(ColorBrewer / viridis / cividis) — critical here because BE3D's signal is *bidirectional* (LFC pos vs
neg), so a diverging map centered at zero is the correct choice, and `matplotlib` defaults are not always
safe; (ii) the **Modifiable Areal Unit Problem** — a choropleth's message changes with binning, exactly
like BE3D's percentile `bin_score`; document sensitivity; (iii) **bivariate maps** to co-encode effect and
confidence (e.g. LFC3D magnitude × pLDDT) rather than hiding low-confidence residues.
- **What BE3D can learn:** Ship a **diverging, colorblind-safe colormap** for on-structure coloring; encode
  **confidence (pLDDT/uncertainty) as saturation or transparency** so viewers can't over-read low-pLDDT
  hotspots; and provide the export as a standard PyMOL/ChimeraX B-factor recolor script (see §3).

---

## 2. CRITIQUE (balanced)

### 2.1 What BE3D does well (genuine strengths)
- **Right idea, real niche.** The central insight — that single-residue BE-screen LFC is noisy and that
  *3D spatial aggregation* surfaces structurally coherent functional hotspots a 1D view misses — is
  sound and matches how every adjacent field (neuro, GIS, ChIP-seq) attacks the same signal shape. It
  occupies a genuine, underserved gap: a **downstream interpretation layer** between functional-genomics
  screens and structural biology, not another BE efficiency/guide/pathogenicity predictor.
- **End-to-end and honest about QC.** BE-QA (Mann-Whitney + KS on knockout-vs-neutral) is a sensible,
  non-optional **quality gate** before interpretation — many tools skip this. Two hypotheses (within-screen
  and vs pooled controls) is thoughtful.
- **Has a null.** It does build a randomization baseline and z-normalizes against it — better than tools
  that eyeball thresholds. It splits +/- effects (bidirectional signal) rather than taking |LFC|, which
  preserves biological direction (LOF vs GOF/resistance).
- **Meta-aggregation is the right lever.** Explicitly combining multiple screens to rescue weak signal,
  with an optional cross-species/isoform alignment step, is exactly where the power is.
- **Interoperable and low-friction to try.** Auto-fetches UniProt + AlphaFold (isoform-aware), runs DSSP,
  supports PPI chains, and exports to the **G2P portal** for interactive 3D viewing. Zero-install Colab +
  YAML-driven local run lowers the barrier. MIT-licensed, modular functions.
- **Structural-feature characterization.** Fisher's-exact enrichment of hits in structural features, and
  the pLDDT/RSA and LFC-vs-LFC3D scatter diagnostics, show awareness that hotspots must be interpreted
  against structure confidence and burial.

### 2.2 Weaknesses / risks (specific, method-grounded)
1. **Twin hard-radius cutoffs (6 Å smoothing + 6 Å clustering) are the biggest methodological weakness.**
   A binary in/out at a fixed radius is scale-blind and discontinuous, and there's no data-driven radius
   selection. Fields that solved this use **distance-decay kernels (Gi\*, kriging), variable-radius scans
   (SaTScan), or threshold-free enhancement (TFCE).** Sensitivity to this knob is likely large and appears
   undocumented.
2. **No FDR / q-values.** Per-residue significance is called at raw p<0.05/0.01/0.001 across *every*
   residue with no multiplicity correction — a real statistical gap for genome/proteome-scale claims.
   Benjamini-Hochberg or Storey q-values are the standard fix (GWAS, MACS, finance all do this).
3. **Clusters get no significance.** Hotspots are found by *geometric* single-linkage clustering of
   already-significant residues; the cluster itself has no calibrated p-value. Cluster-mass permutation
   (Maris-Oostenveld) would fix this.
4. **Single-linkage agglomerative clustering is fragile.** Single-linkage suffers "chaining" — one bridging
   residue can merge two distinct hotspots into one, and results are unstable to the distance threshold.
   No deblending. Consider average/Ward linkage, DBSCAN with a density criterion, or the astronomy
   deblending approach; report cluster stability.
5. **Global null may be mis-specified.** Label/score randomization assumes exchangeability and a
   homogeneous background; it ignores **local** tolerance to mutation (loops vs core) and guide-count
   heterogeneity. A MACS-style **local λ** or covariate-adjusted null would be more defensible.
6. **Naive SUM meta-aggregation.** Summing LFC3D across screens ignores per-screen reliability and
   between-screen heterogeneity. Inverse-variance / empirical-Bayes weighting + an I² heterogeneity report
   are the mature alternatives. Sum also scales with the number of screens (biasing thresholds).
7. **Heavy dependence on AlphaFold confidence.** Hotspots in low-pLDDT / disordered regions or across
   inter-domain interfaces (where AF is least reliable) can be artifacts of geometry, not biology. pLDDT is
   *plotted* but not, apparently, *propagated* into the significance model (e.g., down-weighting neighbors
   by pLDDT). Predicted structures also lack conformational ensembles/allostery.
8. **Fragile external toolchain.** By the authors' own README: **DSSP often fails on large structures**;
   **MUSCLE/CLUSTAL don't build on ARM Macs**, and the MUSCLE API fallback "may also fail." These are
   real, documented adoption blockers (see §3).
9. **Engineering maturity.** No evidence of a test suite or CI; correctness rests on manual example
   reproduction (KBTBD4). Config is a Python/YAML script rather than a documented CLI with schema
   validation. Reproducibility hinges on the `seed` flag (default off) and `nRandom` (500–1000, modest for
   tail p-values like 0.001). Outputs are a sprawl of TSVs across many subdirectories — no single
   machine-readable results manifest.
10. **Usability sharp edges.** Many functions with long, order-dependent, must-match parameters
    (`function_name` must match across `prioritize_by_sequence`/`randomize_sequence`/`calculate_lfc3d`);
    easy to mis-wire silently. This is where LLM/MCP assistance (§3) pays off.

### 2.3 Where to double down + concrete features to add
BE3D should **own the "structural interpretation layer for functional screens"** niche and go deep on
statistical rigor + structural-biology integration rather than broaden. Concrete additions, ranked:
- **Distance-decay kernel (Gaussian) with data-driven length scale** (from a semivariogram/correlogram),
  replacing the hard 6 Å — highest-impact, well-precedented change.
- **Cluster-level permutation p-values + q-values (FDR)** at both residue and cluster level.
- **Local/covariate-adjusted null** (MACS-style local λ; adjust for pLDDT, RSA, guide count).
- **Uncertainty propagation:** per-residue confidence intervals (kriging-style), pLDDT-weighted neighbors.
- **Inverse-variance / empirical-Bayes meta-aggregation + I² heterogeneity** instead of SUM.
- **Robust aggregation option** (median / sigma-clipped) and **deblending / non-single-linkage** clustering.
- **A spatial-autocorrelation pre-flight (global Moran's I)** that tells the user whether 3D aggregation is
  even justified for this target, and suggests the radius.
- **Engineering:** a real CLI (`be3d run config.yaml`), JSON results manifest, unit tests + CI, pinned
  environments, seeded-by-default runs, and a bundled AlphaFold-derived DSSP fallback (mkdssp/`biotite`)
  to kill the DSSP-fails-on-large-structures problem.
- **Visualization:** ship diverging colorblind-safe on-structure coloring + PyMOL/ChimeraX recolor scripts.

---

## 3. ADOPTION INTO SCIENTIST WORKFLOWS (library / Colab / CLI / LLM / MCP)

### 3.1 How a bench/comp scientist folds it in *today*
Typical path: run a BE tiling screen → count guides & compute LFC **upstream in MAGeCK/other** (BE3D does
*not* do read counting; it starts from an LFC TSV with `Mutation_type, Mutation_list, Gene, sgRNA_score`)
→ QC in BE-QA → LFC3D + hotspots in BE-Clust3D → (optional) meta across screens → export to G2P for
interactive 3D viewing / figures.
- **Colab (zero-install):** best for a first look on one gene (KBTBD4 template). Friction: session
  timeouts, re-fetching AlphaFold each run, and slow randomization for large proteins/high `nRandom`.
- **Local + YAML:** recommended for real datasets, but friction is real — conda env with DSSP/MUSCLE/
  CLUSTAL (broken on ARM Macs), DSSP failing on big PDBs, and long chains of functions whose parameters
  must be kept consistent by hand. There is **no single CLI command**; `be3d_local.py <yaml>` is the
  closest thing.
- **Batch/many genes:** possible (functions accept lists, `gene_list=True`) but users must script the
  fan-out and output bookkeeping themselves.

### 3.2 LLM-assisted use (Claude/ChatGPT) — high leverage *now*, no code changes
- **Config authoring:** an LLM can generate/validate the YAML and the parameter-matching (`function_name`,
  `screen_names`, `psig_columns`) that is currently error-prone — resolving weakness 2.2(10).
- **Column mapping:** map a user's arbitrary MAGeCK/screen TSV headers onto BE3D's expected
  `mut_col/val_col/gene_col/edits_col`.
- **Interpretation:** summarize hotspot TSVs into plain-language ("cluster 3 spans residues 210–224, a
  buried β-strand in the Kelch domain, LFC3D-negative → likely loss-of-function core"), and cross-reference
  known domains/mutations. Caution: an LLM will happily over-interpret low-pLDDT hotspots — the tool must
  surface confidence for the LLM to respect.
- **Follow-up generation:** draft PyMOL/ChimeraX coloring scripts, mutagenesis validation lists, and
  methods text.

### 3.3 A "BE3D MCP" server — let an agent *drive* the pipeline
MCP would expose BE3D's modular functions as agent-callable tools with typed schemas + a results manifest,
turning the brittle multi-step script into an orchestrated flow. Sketch of the tool surface:

- `be3d.fetch_structure(uniprot_id, isoform?, chain?)` → paths to AlphaFold/PDB + DSSP (with automatic
  DSSP fallback); returns pLDDT summary.
- `be3d.run_screen_qc(screen_tsv, cases, controls, col_map)` → BE-QA stats (MW-U/KS p-values, effect
  sizes) + a **go/no-go quality verdict** the agent can branch on.
- `be3d.build_features(uniprot_id|fasta+pdb, radius|kernel)` → sequence-structure feature table
  (neighbors, RSA, pLDDT, secondary structure).
- `be3d.compute_lfc3d(features, screen(s), params)` → per-residue LFC/LFC3D + null summary.
- `be3d.call_hotspots(lfc3d, pthr|qthr, clustering_params)` → significant residues + clusters **with
  cluster-level q-values** (the feature to add in §2.3).
- `be3d.meta_aggregate(screens[], method='invvar'|'sum')` → consensus hotspots + heterogeneity (I²).
- `be3d.characterize(hits, features)` → enrichment (Fisher), pLDDT/RSA diagnostics.
- `be3d.export_g2p(hits)` / `be3d.export_pymol(hits)` → G2P TSV / ChimeraX-PyMOL recolor script.
- `be3d.get_hotspots(uniprot_id)` → cached/prior results for a target (agent memory).
Design notes: every tool returns a small **structured JSON summary + artifact paths** (not giant TSVs) so
the agent's context stays clean; long randomization runs should be **async/job-handle** style; include a
`dry_run`/`validate_config` tool so the agent can catch the parameter-mismatch class of errors before a
30-minute run. This maps cleanly onto the compute options in the user's environment (Modal/Explorer for the
randomization-heavy `calculate_lfc3d`/`znorm` steps).

### 3.4 Integrations — current vs next
- **Serves / integrates today:** UniProt (sequence, isoforms), **AlphaFold DB** (structure auto-fetch),
  **DSSP** (secondary structure/RSA), **MUSCLE / CLUSTAL** (conservation alignment), and outputs to the
  **G2P portal** (interactive 3D). Standard PDB input supported.
- **Should integrate next:** (i) **PyMOL / ChimeraX** direct session/recolor export (offline figure
  workflow, complements G2P). (ii) **ESMFold / Boltz** for structures absent from AlphaFold DB and for
  complexes/conformers (the user has monthly Boltz credits — a natural fit). (iii) **MAGeCK / MAGeCK-VISPR
  upstream** so BE3D can ingest counts → LFC directly, closing the gap that today forces users to compute
  LFC elsewhere. (iv) **cBioPortal / ClinVar / COSMIC** to overlay clinical mutation frequency on hotspots
  (does a hotspot recur in patients?). (v) **AlphaMissense / ESM-1v** to compare functional hotspots against
  predicted pathogenicity (orthogonal validation). (vi) A lightweight **dashboard** (Streamlit/Dash or a
  marimo notebook) over the results manifest for non-programmers.

---

## Sources
- [ArcGIS Pro — How Hot Spot Analysis (Getis-Ord Gi*) works](https://pro.arcgis.com/en/pro-app/3.4/tool-reference/spatial-statistics/h-how-hot-spot-analysis-getis-ord-gi-spatial-stati.htm)
- [Wikipedia — Getis–Ord statistics](https://en.wikipedia.org/wiki/Getis%E2%80%93Ord_statistics)
- [SaTScan / Kulldorff spatial scan statistic (NIH grant record)](https://grantome.com/grant/NIH/R01-CA165057-08)
- [A flexibly shaped spatial scan statistic (Int. J. Health Geographics)](https://ij-healthgeographics.biomedcentral.com/articles/10.1186/1476-072X-4-11)
- ["Spatial scan statistics can be dangerous" (PubMed)](https://pubmed.ncbi.nlm.nih.gov/33595399/)
- [MACS3 — Model-based Analysis for ChIP-Seq (docs)](https://macs3-project.github.io/MACS/)
- [MACS2 callpeak manual (dynamic local lambda)](https://manpages.ubuntu.com/manpages/jammy/man1/macs2_callpeak.1.html)
- [Maris & Oostenveld — Nonparametric statistical testing of EEG/MEG data (cluster-based permutation)](https://www.semanticscholar.org/paper/Nonparametric-statistical-testing-of-EEG-and-Maris-Oostenveld/d99ef919dcfa0f089484e505844bb996992459e0)
