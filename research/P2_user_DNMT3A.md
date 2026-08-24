# BE3D Per-User Package — DNMT3A (Lue et al. 2023 base-editor scan)

Synthesizer: "DNMT3A-outreach". Date 2026-08-24. Target user: **epigenetics / AML lab**.
Source run: `real_output/P2_DNMT3A/` (COMPLETE; see `research/P2_run_DNMT3A.md` for the decision log,
sign convention and sweep). BE3D was NOT re-run. All numbers below are from the **primary run r01**
(4U7T chain A, r=6, mean, nRandom=1000, residue-level), LOF = **negated score / `LFC3D_neg` channel**
(activity reporter: DNMT3A methylation silences citrine, so LOF de-represses citrine → we negate so LOF
is captured by BE3D's negative-direction machinery). Analysis script + tables:
`scratchpad/bench.py`, `scratchpad/master_residues.tsv`.

---

## 1. DISCRIMINATION BENCHMARK (A0) — the primary test

**Denominator.** N = **314 scored residues** (chain-A residues 474–912 with a real `LFC3D_neg` value;
the other 598 rows in the padded 912-length table are unscored `-`). This is the population for every
enrichment/rate statistic below. (Matches the run's "residues scored = 314".)

### 1a. BASE RATE — BE3D flags a LOT (honest caveat up front)
| p-threshold | # LOF-significant (K) | base rate K/N |
|---|---|---|
| p<0.05 | 117 | **37.3%** |
| p<0.01 | 107 | 34.1% |
| p<0.001 | 100 | 31.8% |
| clustered (6 Å, p<0.05) | 108 | 34.4% |

**~1 in 3 scored residues is called LOF-significant.** By the OUTREACH_GUIDE rule (base rate >15–20% →
raw overlap is weak evidence), we CANNOT rest on "BE3D recovered site X". The high base rate is partly
real biology — the only crystallized region (474–912) is the **MTase + ADD catalytic domain**, which is
functionally dense — but it means the binary p<.05 flag alone is loose. **Everything credible below comes
from ENRICHMENT, RANK, and the functional-vs-tolerant GAP, not from the flag count.**

### 1b. ENRICHMENT vs chance (Fisher / hypergeometric, one-sided) — BE3D beats chance decisively
FUNCTIONAL set (M = **42** scored residues), curated INDEPENDENTLY of BE3D from 4U7T contacts +
literature (see §5 for the set & citations): catalytic/SAH-cofactor pocket ∪ RD homotetramer interface
∪ FF (DNMT3A–DNMT3L) interface ∪ clinical R882/R736/R729.

| threshold | obs R | expected E=K·M/N | **R/E** | **odds ratio** | Fisher p |
|---|---|---|---|---|---|
| p<0.05 | 35 | 15.7 | **2.24** | **11.6** | 6.0e-11 |
| p<0.01 | 34 | 14.3 | 2.38 | 11.6 | 2.5e-11 |
| p<0.001 | 34 | 13.4 | **2.54** | **13.3** | 2.0e-12 |

Per functional sub-surface (p<.05 hotspots):
| subset | R / size | R/E | OR | p |
|---|---|---|---|---|
| catalytic/SAH pocket | 19/20 | 2.55 | **38.0** | 3.6e-08 |
| FF interface | 6/7 | 2.30 | 10.6 | 1.2e-02 |
| RD interface | 10/15 | 1.79 | 3.6 | 1.8e-02 |

Verdict: enrichment is **strong and highly significant** (OR ≈ 11–13 overall, p ≈ 1e-11), and it **survives
the high base rate**. Catalytic pocket is nearly saturated (OR 38); the **RD interface is the weakest of the
three but still significantly enriched** (OR 3.6, p=0.018) — consistent with the run's finding that RD is the
dominant *quaternary* vulnerability but a softer signal than the active site.

### 1c. DISCRIMINATION GAP — functional ≫ tolerant (the "not everything is a hotspot" test)
BE3D p<.05 hit-rate:
| set | n | hit-rate |
|---|---|---|
| **FUNCTIONAL (all)** | 42 | **83%** |
| — catalytic pocket | 20 | 95% |
| — FF interface | 7 | 86% |
| — RD interface | 15 | 67% |
| **TOLERANT proxy (ADD linker + H3 reader, non-functional)** | 41 | **24%** |
| — ADD_MTase linker | 26 | 27% |
| — H3-reader surface | 15 | 20% |
| — (curated) DNA face | 9 | 11% |
| all other scored | 231 | 31% |

**Discrimination gap ≈ 83% vs 24% (~3.4×).** Functional residues are hit far more often than tolerant
ones — BE3D is discriminating, not blanket-calling. Note the structural-tolerant proxy uses the
run's *depleted* surfaces (ADD linker, H3 reader); a purpose-built ClinVar-benign set could not be
retrieved residue-by-residue (worker flagged this gap — pull the DNMT3A ClinVar Benign/LB missense
set to harden this).

### 1d. PRECISION@K (rank by strongest |LFC3D_neg|)
Base expectation M/N = 0.134.
| K | precision@K | lift vs base |
|---|---|---|
| 10 | **0.60** | **4.5×** |
| 20 | 0.55 | 4.1× |
| 30 | 0.53 | 4.0× |

The strongest signals concentrate ~4–5× more than chance on functional residues — the ranking is
trustworthy even though the flag is loose.

### 1e. BEAT A TRIVIAL BASELINE — burial has ZERO functional signal; BE3D wins outright
Baseline = "call the most-buried residues" (contact number within 10 Å from 4U7T chain A, top third,
n=95). Burial → FUNCTIONAL: R=14, E=14.9, **R/E=0.94, OR=0.89, p=0.69 (NOT significant)**;
burial precision@10 = 0.10, @20 = 0.25. **Burial does not predict functional residues at all** (the
functional set here is dominated by *surface* interface residues + a partly-exposed cofactor pocket).
Crucially, **BE3D hits 79% of NON-buried functional residues** (p<.05) — exactly the surface interface
residues a burial/conservation predictor misses. AlphaMissense-high baseline: could not be computed —
per-residue AlphaMissense scores for DNMT3A were not retrievable (worker gap; download the
AlphaMissense catalogue filtered to Q9Y6K1 to add this). Qualitatively AlphaMissense is expected to
call the whole 634–912 catalytic domain intolerant, so it would recover the catalytic signal but, being
screen-free, cannot weight the *editable, activity-relevant* interface residues the way BE3D does.

### 1f. HONEST VERDICT
**BE3D discriminates — it does not merely over-call.** Enrichment OR ≈ 11–13 (p ≈ 1e-11), functional-vs-
tolerant gap ≈ 83% vs 24%, precision@10 = 0.60 (4.5× lift), and it **beats the burial baseline outright**
(burial p=0.69). **BUT** the absolute base rate is high (37% at p<.05), so the *binary* significance flag is
loose and should never be quoted alone — use the RANK and the per-surface enrichment. The signal is
**concentrated, not diffuse**: the catalytic/SAH pocket dominates (OR 38, 95% hit-rate — the expected,
paper-level result for an activity screen), and among protein–protein surfaces the **RD homotetramer
interface is the dominant editable vulnerability** (robust but softer, OR 3.6). This directly answers the
run's "SAH pocket dominates overall" — quantified: catalytic ≫ FF ≈ RD ≫ (ADD, H3, DNA depleted).

---

## 2. FALSE POSITIVES / FALSE NEGATIVES — the honest story

### R882 — the headline FALSE NEGATIVE (state prominently)
**R882 is the #1 somatic DNMT3A hotspot in AML (~50–60% of DNMT3A-mutant AML; 53% in AML vs 10.6%
in CHIP — Venugopal 2022, PMID 35296003) and BE3D does NOT flag it** (LFC3D_neg z = +2.9, NS;
significant in only **1/12** sweep runs, and only at r=10 via distant neighbors). This is **not a BE3D failure
of intelligence — it is a hard limit of the assay**: **R882H is a G→A codon change that a cytosine base
editor (CBE) cannot install**, and no R882 guide exists in the Lue screen. BE3D can only see base-editable
positions. So at the residue level R882 is **base-editing-invisible** → a true false negative. What BE3D
*does* recover is the **same RD tetramer interface via editable neighbors** (see below): interface-level
co-localization YES, residue-level NO.

### The NOVEL-CANDIDATE proxy — clinically-unannotated editable RD patch
BE3D flags a compact, robust, **clinically UNANNOTATED** RD-interface patch spatially adjacent to R882 but
distinct from it: **T671, M674, V675, R676, H873, D876, V877, S878** (T671/D876/V877/S878 are p<0.001 in
**12/12** sweep runs; M674/V675/H873 in ≥10/12). None carry documented patient mutations (worker check:
R676K/M674T/N879A appear only as engineered probes). Mechanistically these line the RD–RD
oligomerization surface (the R882–D876 salt-bridge network; R676–E820′; M674 — Emperle 2018 NAR
PMID 29506156). **Hypothesis (not a claim):** editing this patch phenocopies the R882H tetramerization
defect *without* touching R882 — a clean, base-editable **proxy allele series** for the dominant-negative
codon that CBE cannot reach, and a candidate **tetramerization-disrupting drug epitope**.

### Possible false positives / lower-confidence
Genuine FP candidates are **few** — the crystallized region is functionally dense. Of the 117 p<.05
residues, 43 are catalytic pocket + catalytic-core motif residues (630–800, strong z ≤ −9), 8 are FF, 12 are
RD, 8 are real out-of-scope PWWP/N-term LOF hits (G293/E294/S337 — PWWP reads H3K36me3, mutated
in AML/TBRS; real hits, just outside the 4U7T 474–912 window), leaving ~62 residues flagged
**LOWER_CONFIDENCE** (p<.05, weaker z, no annotated element; some are catalytic-core neighbors, others
may be neighbor/bystander-driven). We flag these explicitly rather than assert function. No single hotspot
sits in an obviously benign/disordered context — pLDDT is high across the crystallized catalytic domain.

### Precision/recall FEEL (honest, not a fake number)
Of the robust hotspots: the catalytic pocket + FF + RD interface calls are **TP against structure/literature**;
the T671/674–678 + 873–878 RD patch is the **NOVEL-CANDIDATE** (8 residues); ~62 weaker p<.05 residues
are **lower-confidence / possible FP**; and **R882 is the notable FALSE NEGATIVE** (assay-limited). Recall of
*base-editable* functional residues is high (83% of the functional set); recall of the clinically dominant
residue (R882) is zero for the mechanistic reason above.

---

## 3. COMPETITOR COMPARISON (B)

| Tool | Input signal | What it finds on DNMT3A | Advantage | Disadvantage |
|---|---|---|---|---|
| **cBioPortal / 3dhotspots.org / COSMIC** (somatic recurrence) | Patient tumor mutation frequency (HotMAPS-style 3D clustering) | **R882 (dominant), R736; RD/DNA-binding surface cluster** | Decade of clinical benchmarking; **WINS at R882** — flags it trivially as #1 AML residue; disease-anchored | Recurrence ≠ function; blind to base-editable-but-not-recurrent residues; can't attribute *which interface* drives LOF; needs large cohorts |
| **AlphaMissense** (Cheng 2023, Science) | Sequence + AF structure, unsupervised pathogenicity | Whole 634–912 catalytic domain ~uniformly "likely pathogenic"; interface residues plausibly high | No screen needed; genome-wide; orthogonal predictor | Per-residue DNMT3A scores not retrieved here (gap); **can't rank by editability or activity**; low-resolution across a uniformly-intolerant domain |
| **ProTiler-Mut** (Cell Systems 2026; 3D-RRA) | Any CRISPR *tiling* screen (BE/HDR/PE/SGE) → 3D rank-aggregation + PPI-disruption inference | Would cluster the same catalytic+interface signal; could infer unscreened-mutation effects | Screen-agnostic; if fed HDR/PE covers substitutions CBE can't (potentially **including R882**); explicit PPI-disruption calls | On a CBE screen it inherits the SAME R882 blind spot; method, not a precomputed DNMT3A result |
| **BE3D (this run)** | Base-editor *activity* screen LFC → 3D neighbor null + randomization meta | Catalytic pocket (OR 38) ≫ **RD editable interface** (the dominant PPI vulnerability) ≫ FF; the T671/873–878 R882-proxy patch | Functional screen signal + controls; **resolves the interface** recurrence tools lump together; names *editable* handles | **Misses R882** (CBE can't make R882H); high base rate; interface signal softer than catalytic; needs a screen |

**Fair bottom line:** for R882 itself, **somatic tools win, unambiguously** — it is their textbook strength and a
CBE screen structurally cannot reach it. BE3D's complementary value is **interface resolution + editability**:
it converts the diffuse "interdomain / RD-DNA cluster" that recurrence and AlphaMissense report into a
ranked, per-surface map (RD ≫ FF ≫ depleted ADD/H3/DNA) and names a **base-editable R882 proxy** that
recurrence-based tools cannot attribute and that provides an experimentally cleaner allele than the
confounded R882 codon.

---

## 4. OUTREACH DELIVERABLES (C, D)
Files assembled in `outreach/DNMT3A/`:
- `DNMT3A_hotspots.tsv` — 117 p<.05 LOF hotspots + R882 (false-negative row); columns: residue, aa,
  LFC3D_neg, |LFC3D_neg|, z, p_best, surface, **flag** (TP_catalytic / TP_catalytic_core / TP_FF /
  TP_FF_clinical / TP_RD / **NOVEL_CANDIDATE** / TP_outofscope_PWWP / LOWER_CONFIDENCE /
  **FALSE_NEGATIVE**), n_runs_p001_of12 (robustness), direction, rationale.
- `DNMT3A_G2P_LOF_neg.tsv` — the run's NEG-channel G2P table for the interactive 3D viewer.
- `DNMT3A_email.md`, `DNMT3A_brief.md` — copyable outreach.

---

## 5. FUNCTIONAL / TOLERANT SETS + citations (curated independently of BE3D)
**Numbering:** UniProt Q9Y6K1 canonical (912 aa). Domains (UniProt): PWWP 292–350; ADD 482–614;
SAM-dependent MTase 634–912. **Active-site catalytic cysteine = C710** (UniProt ACT_SITE; note the run
log's "C666/C706" is a numbering slip — C710 is the nucleophile). SAM/SAH cofactor-pocket binding
residues: **641–645, 664, 686–688, 891–893** (UniProt). CpG target-recognition contact R836 (Zhang 2018).

- **FUNCTIONAL (M=42 scored):** catalytic/SAH pocket (4U7T <6 Å SAH; UniProt cofactor residues) ∪ RD
  homotetramer interface (4U7T contacts; R882–D876 salt bridge, R676–E820′, M674 — Emperle 2018 NAR
  PMID 29506156) ∪ FF DNMT3A–DNMT3L interface (F732/R736/E733 anchors — Emperle & Jeltsch 2023
  PMID 36528185; Jia 2007 PMID 17713477) ∪ clinical R882/R736/R729.
- **TOLERANT (structural proxy):** ADD–MTase linker + H3-reader surface residues (the run's depleted
  surfaces), non-functional. NOTE: a curated ClinVar-Benign set was not retrievable per-residue (gap).
- **Clinical anchors:** R882H dominant somatic AML hotspot, dominant-negative by blocking active-tetramer
  formation (**Russler-Germain 2014 Cancer Cell PMID 24656771**; ~80% activity loss); in-vitro nuance —
  Nguyen 2018 Sci Rep PMID 30202094 disputes pure DN; Sandoval 2019 JBC PMID 30705090 (LOF+GOF).
  **R736H = 2nd-most-frequent DNMT3A cancer residue** (FF interface; exact % unconfirmed — do not
  quote "~2%"). R729W abolishes processive catalysis. **PDB 4U7T = Guo et al. 2015 Nature** (autoinhibition;
  DNMT3A–DNMT3L–H3) — NOT Zhang 2018 (that is 5YX2, the CpG-DNA complex; correction vs run log).
- **Germline:** Tatton-Brown–Rahman syndrome = germline LOF DNMT3A across all domains (PMID 24614070).

### Caveats / open items
(i) 4U7T covers 474–912 only — PWWP/N-term hits (G293/E294) are real but out of interface scope.
(ii) DNA face is literature-curated (no DNA in 4U7T) → its low enrichment is soft; a 5YX2/6W89 run would
harden it. (iii) FF-vs-RD ranking is radius-sensitive (FF abuts the catalytic core); the RD call is the more
conformation-independent one (holds in the AF baseline). (iv) Gaps to close before publication: exact
3dhotspots.org DNMT3A cluster list; per-residue AlphaMissense scores; curated ClinVar-benign tolerant set.
