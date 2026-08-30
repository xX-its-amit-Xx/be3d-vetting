# BE3D — Cross-Domain Technique Inspiration (Phase 5)

A fan-out across ~25+ fields for mature, named techniques that map onto BE3D's pipeline
(*noisy 1D-sequence measurements → 3D point cloud → kernel smoothing → permutation null → hotspot
clustering → multi-screen consensus*). The striking result is **convergence**: unrelated fields keep
prescribing the *same* fixes for BE3D's known gaps (over-calling / no-FDR, hard 6Å cutoff, naive SUM
meta, no uncertainty, no cluster-level significance, reachability). Citations verified via
Crossref/DOI/Semantic Scholar/ACL.

## Priority stack (by cross-domain convergence)

### 1. Empirical-Bayes / shrinkage — the strongest signal (≥6 independent fields)
Fixes the ~7–50% over-calling (regression to the mean), replaces naive SUM meta, and supplies
per-residue uncertainty — one family, three gaps. Ready-made recipes:
- Highway safety — **Hauer** EB crash-hotspot method (TRR 1784, 2002) — site-level.
- Disease mapping — **Marshall** spatial/local EB rate smoothing (JRSS-C, 1991) — spatial.
- Sports/statistics — **James–Stein / Efron–Morris** shrinkage.
- Actuarial — **Bühlmann credibility** Z = n/(n+k) (closed-form, tuning-free) — cleanest drop-in for
  per-screen weighting.
- Phylogenetics — **FUBAR / Rate4Site** empirical-Bayes posterior per site.
- Deep mutational scanning — **Enrich2** inverse-variance random-effects meta (replaces SUM).
→ Shrink each residue's smoothed score toward a local/global expected background; weight by
measurement variance and pLDDT.

### 2. Spatial-FDR on a local cluster statistic — fixes multiplicity + cluster-level significance
- **Getis–Ord Gi\*** per-neighborhood z-score (1992) + **Caldas de Castro & Singer** spatial-dependence-
  aware BH-FDR (2006) — flagged independently by ecology, cartography, linguistics.
- Single-cell genomics converges: **Milo SpatialFDR** over overlapping kNN neighborhoods — the same idea
  on a 3D point cloud.
→ Replace raw p<0.05 with a Gi\*/Milo-style per-neighborhood statistic under spatially-aware FDR.
(Directly extends PR #20's residue-level BH q-values to the cluster level.)

### 3. Change-point detection for bystander runs
- **CUSUM** (Page, 1954) over residue-indexed scores distinguishes a sustained shift from one smoothing
  neighborhood (bystander run) from genuinely independent 3D hotspots; **Foote** novelty / local
  median+MAD (Bello 2005) for adaptive peak-picking.

### 4. GP / kernel smoothing with calibrated variance — replaces the hard 6Å cutoff + adds uncertainty
- **GP occupancy maps** (O'Callaghan & Ramos 2012) — squared-exponential length-scale by marginal
  likelihood; pLDDT enters as observation noise; posterior variance everywhere.
- **RTS/Kalman smoother** (1960/1965) — along-sequence denoising with per-residue posterior SD.
- MAGIC diffusion kernel; KDE + Silverman bandwidth; **dasymetric weighting** (Mennis 2003) to
  down-weight neighbors by pLDDT; Ripley's K for multi-scale significance.
→ Swap the 6Å indicator for a Gaussian/GP kernel with data-driven bandwidth + pLDDT-scaled noise,
emitting mean ± variance per residue.

### 5. Occupancy models with imperfect detection — a near-exact model for the reachability gap
- **MacKenzie occupancy models** (2002): a residue can be functional (ψ) yet undetected because no
  CBE/ABE installs the needed substitution (detection p→0). Model p per residue from editable-window
  presence + achievable substitutions, so unreachable residues score "missing/uninformative," not
  "no effect." Nothing else surveyed models base-editing coverage this precisely — it's the principled
  successor to PR #19's reachability report.

## Secondary: robust aggregation
Weighted head-banging median smoothing (Mungiole 1999) and local median+MAD (Bello 2005) — outlier-
resistant replacements for mean/SUM, complementary to the EB/credibility methods.

## Roadmap → next PRs (after the 4 already opened)
| Rank | Idea | Source domains | Gap addressed | Effort |
|---|---|---|---|---|
| 1 | Empirical-Bayes shrinkage (Bühlmann/Hauer/Marshall) | safety, actuarial, disease-mapping, DMS, phylo | over-calling + SUM meta + uncertainty | M |
| 2 | Spatial-FDR on Gi*/Milo cluster statistic | GIS, single-cell, ecology | multiplicity + cluster significance | M |
| 3 | GP/RTS kernel + pLDDT-as-noise | robotics, geostats, controls | hard 6Å cutoff + uncertainty | L |
| 4 | Occupancy-model detectability layer | ecology | reachability (principled) | M |
| 5 | CUSUM bystander-run separation | quality control, MIR | bystander vs independent hotspots | S |

```json
[
 {"rank":1,"idea":"Empirical-Bayes / credibility shrinkage of per-residue scores","domains":["highway safety (Hauer)","actuarial (Buhlmann)","disease mapping (Marshall)","DMS (Enrich2)","phylogenetics (FUBAR)"],"gap_addressed":"over-calling + naive SUM meta + no uncertainty","effort":"M","pr_candidate":true},
 {"rank":2,"idea":"Spatial-FDR on a Getis-Ord Gi*/Milo per-neighborhood cluster statistic","domains":["GIS","single-cell (Milo)","ecology"],"gap_addressed":"no multiplicity control + no cluster-level significance","effort":"M","pr_candidate":true},
 {"rank":3,"idea":"GP/RTS kernel smoothing with data-driven bandwidth and pLDDT-as-observation-noise","domains":["robotics (GP occupancy)","geostatistics","control theory (Kalman/RTS)"],"gap_addressed":"hard 6A cutoff + no uncertainty","effort":"L","pr_candidate":true},
 {"rank":4,"idea":"Occupancy-model detectability layer for base-editing reachability","domains":["ecology (MacKenzie)"],"gap_addressed":"reachability / imperfect detection (principled successor to PR#19)","effort":"M","pr_candidate":true},
 {"rank":5,"idea":"CUSUM change-point to separate bystander runs from independent 3D hotspots","domains":["statistical process control","music information retrieval"],"gap_addressed":"bystander-edit contiguous runs","effort":"S","pr_candidate":true}
]
```

## Bottom line
The priority stack — (1) empirical-Bayes shrinkage, (2) spatial-FDR on a Gi*/Milo statistic, (3) GP/RTS
kernel with pLDDT-as-noise, (4) occupancy-model detectability, (5) CUSUM — is a coherent, well-precedented
research roadmap. Notably, the two most novel BE3D-specific angles we shipped (PR #19 reachability, PR #20
FDR) each have a principled cross-domain successor (occupancy models; spatial-FDR at the cluster level),
so the roadmap is continuous with the code already contributed.
