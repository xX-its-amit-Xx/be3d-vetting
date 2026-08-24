# BE3D Outreach Package — KRAS(G12C) resistance (sotorasib vs adagrasib)

Persona: **KRAS(G12C) drug-discovery / med-chem lab.** Source run: completed BE3D grid,
`real_output/P2_KRAS/` (baseline AF_r6_mean_n1000_res, canonical KRAS4B numbering 1–189,
N=189 residues scored). Direction convention: **positive LFC3D = enrichment = drug RESISTANCE**.
Analysis scripts: `scratchpad/discrim.py`, `scratchpad/make_tsv.py` (venv gone; used system Python
3.14 + pandas/scipy). This doc is the benchmark; deliverables in `outreach/KRAS/`.

---

## 1. DISCRIMINATION BENCHMARK (A0 — the primary deliverable)

### 1a. Base rate (is BE3D just flagging everything?)
Positive/resistance direction, AlphaFold baseline run. Denominator = all 189 structured residues
(each residue is assigned to exactly one direction; 67 residues are positive-leaning for sotorasib,
93 for adagrasib).

| Drug | p<0.05 | p<0.01 | p<0.001 |
|---|---|---|---|
| **Sotorasib** #sig / base rate (all 189) | 24 / **12.7%** | 20 / 10.6% | 17 / 9.0% |
| Sotorasib base rate (of 67 pos-direction) | 35.8% | 29.9% | 25.4% |
| **Adagrasib** #sig / base rate (all 189) | 30 / **15.9%** | 27 / 14.3% | 25 / 13.2% |
| Adagrasib base rate (of 93 pos-direction) | 32.3% | 29.0% | 26.9% |

**Read:** BE3D calls ~13–16% of the whole protein a positive resistance hotspot at p<0.05. That is a
**moderate** base rate — not the "everything is a hotspot" pathology (would be >30–50% of the whole
protein), but high enough that raw overlap with known sites is weak evidence on its own. The
enrichment test below is what matters. (Counts reproduce the decision log exactly: Sotor 24, Adag 30.)

### 1b. Enrichment vs chance (the real test) — hypergeometric / Fisher, N=189
Two INDEPENDENT functional sets, curated without reference to BE3D:
- **FUNCTIONAL-LIT** (M=9): validated G12Ci resistance residues **T58, S65, R68, D69, M72, D92, H95,
  Y96, Q99** (Feng 2022 *PNAS* PMID 35471904; Awad 2021 *NEJM* PMID 34161704; Tanaka 2021 *Cancer
  Discov* PMID 33846219).
- **FUNCTIONAL-POCKET/SWII** (M=60): switch-II 59–76 ∪ all residues <8 Å from either bound drug
  (structural SII-P lining, from `analysis/pocket_dist.json`, 6OIM/6UT0 contacts). SII-P ≈ switch-II
  (59–76) + α3 groove; drugs bind covalently at C12 (standard KRAS structural convention).

E = K·M/N (expected overlap). **BE3D discriminates only if R/E ≫ 1 with significant p.**

| Drug @ p<0.05 | FUNCTIONAL set | K | M | obs R | exp E | **R/E** | odds ratio | **p (Fisher=hyperg.)** |
|---|---|---|---|---|---|---|---|---|
| **Adagrasib** | LIT (9 resistance res) | 30 | 9 | **8** | 1.43 | **5.6×** | 57.5 | **1.4×10⁻⁶** |
| **Adagrasib** | POCKET/SWII (60) | 30 | 60 | **27** | 9.52 | **2.83×** | 34.4 | **4.9×10⁻¹³** |
| **Sotorasib** | LIT (9 resistance res) | 24 | 9 | 3 | 1.14 | 2.62× | 3.8 | 0.091 (NS) |
| **Sotorasib** | POCKET/SWII (60) | 24 | 60 | 8 | 7.62 | **1.05×** | 1.09 | **0.51 (NS)** |

Holds across thresholds (adagrasib p<0.01: LIT R/E 6.2× OR 68 p=5×10⁻⁷; POCKET R/E 3.0× OR 98
p=2×10⁻¹⁴). **Sotorasib never reaches enrichment significance for the canonical pocket** (R/E≈1.1,
p≈0.46–0.51 at every threshold; LIT stays borderline, p≈0.06–0.19).

### 1c. Discrimination gap (functional hit-rate vs tolerant hit-rate)
TOLERANT set (M=53): coordination-number surface residues (not functional, not <8 Å to pocket, P-loop
excluded) + HVR 170–189. Hit-rate = fraction of the set BE3D calls at p<0.05.

| Drug | hit-rate FUNCTIONAL(union,60) | hit-rate TOLERANT(53) | **gap** |
|---|---|---|---|
| **Adagrasib** | **45.0%** | **1.9%** | **~24×** (p<0.01: 43% vs 0% → ∞) |
| **Sotorasib** | 13.3% | 9.4% | **1.4×** |

Adagrasib: clean, huge gap → real signal. Sotorasib: functional ≈ tolerant → **weak/no
discrimination on the canonical pocket axis** (its "tolerant" hits are 47/48 loop, 154 = α5, 178/179 =
HVR — i.e. partly the novel candidate, partly artifact; see §2).

### 1d. Precision@K (top residues by |LFC3D_pos|)
| Drug | prec@10 (FUNC union) | prec@10 (LIT) | prec@20 (FUNC) | prec@20 (LIT) | chance (M/N) |
|---|---|---|---|---|---|
| **Adagrasib** | **1.00** (10/10) | 0.30 | **0.95** | 0.30 | 0.32 / 0.05 |
| Sotorasib | 0.40 | 0.10 | 0.35 | 0.15 | 0.32 / 0.05 |

Adagrasib's strongest signals concentrate perfectly on the pocket (top-10 all pocket-lining,
|LFC3D| up to 2.7 at A66). Sotorasib's top signals are split between pocket (58,60,66,67) and the
distal α5/HVR (148,149,152,178,179) → modest precision.

### 1e. Beat a trivial baseline
Naive "call all buried residues" (CA coordination top tertile, 70 residues):
- vs FUNCTIONAL-union: R/E **1.3×**, OR 2.0, p=0.02 (weakly enriched).
- vs FUNCTIONAL-LIT (9 resistance res): R/E **0.6×**, OR 0.47, p=0.91 — **buried MISSES the resistance
  residues entirely** (they are surface pocket residues, not buried core).
- **AlphaMissense** (Cheng 2023 *Science* PMID 37733863): KRAS is highly constrained → uniformly high
  pathogenicity across the core and the activation codons; drug-agnostic, cannot separate resistance
  from activation, and cannot produce the sotorasib≠adagrasib asymmetry.

**BE3D's adagrasib map (LIT R/E 5.6×, OR 57) beats the buried baseline (R/E 0.6×) by ~2 orders of
magnitude** on the residues that actually matter, and adds the drug-specific axis neither baseline has.

### 1f. HONEST VERDICT
**BE3D discriminates — but drug-dependently, and this asymmetry is itself the result.**
- **Adagrasib: unambiguous, strong discrimination.** R/E 3–6×, p down to 2×10⁻¹⁴, 24× functional/tolerant
  gap, precision@10 = 1.0, recovers 8/9 validated resistance residues zero-shot, and crushes both
  trivial baselines. This is a genuine, well-controlled recovery — not over-calling.
- **Sotorasib: BE3D does NOT significantly enrich the canonical switch-II pocket** (R/E≈1.05, p=0.51;
  gap 1.4×). Its positive signal is real but sits **off** the annotated pocket — dominated by the
  distal α5 147–156 cluster (novel candidate) + an HVR artifact + a β2–β3 loop (46–48). Two honest
  readings, both stated to the user: (i) sotorasib's true in-pocket resistance footprint is genuinely
  sparse and its escape geometry is more distal (the scientifically interesting hypothesis), or
  (ii) the sotorasib arm over-calls outside the pocket. The moderate base rate (13%) + null pocket
  enrichment means **for sotorasib we cannot claim discrimination on known sites** — we report the
  distal cluster as a hypothesis, not a validated recovery.

Bottom line for the med-chem lab: trust the **adagrasib orthosteric-shell map** (benchmarked, strong);
treat the **sotorasib distal-α5 finding** as a testable hypothesis, not a benchmarked recovery.

---

## 2. FALSE POSITIVES / NEGATIVES (per-residue, `outreach/KRAS/KRAS_hotspots.tsv`)

**TP (validated), recovered zero-shot:** adagrasib hits **S65, R68, D69, M72, D92, H95, Y96, Q99**
(8/9 LIT resistance residues; H95/Y96/Q99 in the adagrasib-specific His95/Q99 groove, mechanistically
correct — they are the adagrasib-only escapes). Sotorasib hits **T58, S65, R68** (+ pocket 60,66,67).

**NOVEL-CANDIDATE (hypothesis, not TP/FP):** sotorasib **α5/SAK cluster 147–156** (K147,T148,R149,
Q150,R151,V152,E153,D154,A155,F156). High pLDDT (**97.2–98.8**, well-structured), ~13–27 Å from the
drug (distal/allosteric), 148/154 persist 7/7 AF runs. One prior datapoint: **F156L** (Feng 2022, no
assigned mechanism) — our cluster gives it spatial context. Candidate mechanism (unproven): SAK/α5
perturbs nucleotide loading / switch-II equilibrium that sotorasib (inactive-state trap) depends on →
would explain sotorasib-specificity. Caveat: Whaby 2022 *JBC* PMID 36334633 shows most α4–α5 point
mutants don't impair signaling, so a simple dimerization mechanism is disfavored. **Validate, don't claim.**

**LIKELY FALSE POSITIVE / artifact:** **HVR 178–180** (P178,G179 in both drugs) — pLDDT **45–49**
(disordered HVR), no crystal coverage, geometry unreliable → flag as passenger/artifact, not a hotspot.

**UNCERTAIN (possible FP / bystander):** sotorasib **46–48** (β2–β3 loop, 13–18 Å from drug, no known
resistance function) and **142**; adagrasib **90** (F90, pocket-adjacent bystander) and **127**
(α4, single p<0.05). These inflate the sotorasib base rate and depress its enrichment.

**FALSE NEGATIVES:**
- Adagrasib **missed T58** (in LIT set; likely base-editing accessibility / guide coverage — the codon
  change may be ABE/CBE-inaccessible or under-covered).
- Sotorasib "missed" **D92, H95, Y96, Q99** — but these are **adagrasib-specific** (His95/Q99 groove);
  their absence from the sotorasib map is **mechanistically correct (true negatives)**, not misses.
- Sotorasib **Y96**: reported pan-inhibitor (Awad 2021) but sotorasib contacts Y96 only via water
  bridges; its absence is a genuine partial false-negative / weak-signal case worth noting.
- Not recovered (correctly): activating mutations **G12/G13/A59/Q61** fall in BE3D's essential/depleted
  (negative) direction, as expected for in-cis activators, not orthosteric escapes.

Precision/recall FEEL: of the ~30 robust adagrasib hotspots, ~27 are on the validated pocket, ~2–3
uncertain, ~1 artifact → high precision. Of the ~24 sotorasib hotspots, ~6 pocket (validated context),
~10 novel-candidate (α5), ~5 uncertain/loop, ~2 HVR artifact → precision on *known* sites is low but
the α5 cluster is a coherent, high-confidence-geometry hypothesis.

---

## 3. COMPETITOR COMPARISON (B)

| Tool | Input signal | What it finds on KRAS | Advantage | Disadvantage vs BE3D |
|---|---|---|---|---|
| **3dhotspots.org / cBioPortal** (Gao 2017 *Genome Med* PMID 28115009, HotMAPS-style) | Patient tumor mutation **recurrence** | **G12, G13, Q61, A146** — oncogenic **ACTIVATION** codons | Decade of benchmarking; no screen needed; robust for drivers | **Orthogonal axis** — finds activation, **not** inhibitor resistance; blind to sotorasib vs adagrasib escape geometry; no drug dimension |
| **AlphaMissense** (Cheng 2023 *Science* PMID 37733863) | Sequence/structure **conservation** | Uniformly high pathogenicity across KRAS core + activation codons; pocket residues heterogeneous | No screen needed; proteome-wide; strong on constrained core | Drug-agnostic; **cannot separate resistance from activation from LoF**; no drug-specific map; can't nominate α5 |
| **ProTiler-Mut** (Cell Systems 2026; bioRxiv 10.1101/2025.04.17.649336; orig ProTiler Nat Commun 2019 PMID 31586052) | Tiling mutagenesis → **3D-RRA** (robust rank aggregation) | Would find functional 3D hotspots + separation-of-function class from tiling screens | Direct 3D tiling rival; SoF class ClinVar-enriched; rank-based null | Not framed as drug-resistance-vs-sensitivity or inhibitor-specific escape; conceptual only (no KRAS G12Ci output pulled) |
| **BE3D** (this run) | Base-editing **resistance screen** LFC3D + randomization null + 3D clustering | **Drug-specific resistance atlas**: adagrasib orthosteric shell (60–72, 89–103, R/E 5.6×, p 1e-6); sotorasib sparse pocket + **distal α5 147–156 candidate** | **Only tool that reads inhibitor-ESCAPE geometry**; separates sotorasib vs adagrasib (Jaccard 0.12); recovers 8/9 validated adagrasib resistance res zero-shot | Needs a resistance screen; sotorasib arm noisy/under-enriched; per-guide causality unresolved; HVR artifact |

**The clear BE3D advantage (for this med-chem lab):** recurrence tools (3dhotspots/cBioPortal) and
AlphaMissense find *where KRAS is mutated to become oncogenic* (G12/G13/Q61) — an axis **orthogonal**
to what a G12Ci program needs. BE3D is the only one that maps *where the drug's grip fails*, and it
resolves **inhibitor-specific** escape geometry (adagrasib's dense His95/Q99 groove vs sotorasib's
shallow footprint) that no recurrence- or conservation-based tool can see. Be fair: for oncogenic
driver calling and for zero-screen prediction, the competitors win; BE3D needs a functional screen.

---

## 4. OUTREACH DELIVERABLES
- `outreach/KRAS/KRAS_hotspots.tsv` — 54 rows, per-residue: residue, AA, drug, LFC3D_pos, p, sig level,
  min-dist-to-drug, pLDDT, **flag (TP / TP-pocket / NOVEL-CANDIDATE / LIKELY-FP / UNCERTAIN)**, note.
- `outreach/KRAS/KRAS_sotorasib_g2p.tsv`, `KRAS_adagrasib_g2p.tsv` — G2P interactive-3D input (baseline run).
- `outreach/KRAS/KRAS_email.md`, `KRAS_brief.md` — copyable outreach.
