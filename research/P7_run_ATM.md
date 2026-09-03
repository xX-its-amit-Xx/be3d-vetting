# BE3D Run — ATM (huge PIKK; the "big / low-confidence protein" stress-test)

> Runner: Queen "ATM/CHEK2". Real public DDR base-editing screen (Cuella-Martin 2021 *Cell*) →
> BE3D LFC3D on a **3056-aa** PIKK kinase. Purpose: probe **where BE3D struggles on very large,
> partly-flexible multi-domain proteins** — structure availability, runtime/RAM of the neighbor
> search, IDR/over-call false positives, and whether signal concentrates on the folded kinase/FAT
> core. Central deliverable = the **A0 DISCRIMINATION BENCHMARK** + an explicit struggle audit.

## 0. Target, data, direction — and the first struggle
- **Gene/protein:** ATM, UniProt **Q13315**, **3056 aa**, PIKK. Domains (UniProt + cryo-EM): N-terminal
  **HEAT/α-solenoid ≈ 1–1898** (spiral+pincer), **FAT 1940–2566**, **PI3K/PI4K kinase 2686–2998**
  (catalytic core ~2712–2962), **PRD 2999–3023** (R3008), **FATC 3024–3056**.
- **STRUGGLE #1 — no AlphaFold model exists.** `AF-Q13315-F1-model_v6/v4` → **HTTP 404**; the AF-DB API
  returns **`{}`** (proteins >2700 aa are excluded from the AF-DB file store). **BE3D's UniProt→AlphaFold
  auto-fetch fails outright for ATM.** We supplied an experimental structure instead: **cryo-EM 7NI5 chain A**
  (human ATM; RCSB), extracted to a single chain (`ATM_7ni5_chainA.pdb`) — **2791/3056 residues modeled
  (91%, span 3–3056)**; the **265 unmodeled residues (9%) are the flexible/disordered ones** and are simply
  absent from scoring. NOTE: with an experimental PDB, BE3D reads the **B-factor column as "pLDDT"** — this
  is **wrong** (it's a real cryo-EM B-factor). We therefore use B-factor **correctly** ourselves (high B =
  poorly-ordered/flexible) and flag the mislabel as a bug (§4).
- **Premise correction (honest):** ATM is **NOT "mostly low-pLDDT/disordered"** — it is a **mostly-FOLDED
  α-solenoid**: 91% is resolved ordered helical repeats; disorder is limited to ~9% unresolved + flexible
  loops. The relevant challenge is **conformational flexibility + tolerant surface**, not a BRCA1-style IDR.
- **Screen:** Cuella-Martin 2021 *Cell* (PMID 33592168), CBE tiling, Table **S5** `cm_mmc5.xlsx`, sheet
  **`MCF10A treatments Subl1`**. **684 ATM guides**, `AAChg` per-token parsed (recipe = P2_run_BRCA1 §1).
- **Assay direction (documented):** like CHEK2 and **opposite to BRCA1** — **nonsense guides are ENRICHED**
  (CISP nonsense mean **+0.66** vs synonymous **+0.23**, KS D=0.41 p=7e-5). ATM-LOF ablates the DNA-damage
  checkpoint → survival under a genotoxin → **LOF = POSITIVE LFC**. **NEGATED** so LOF → BE3D `neg` channel.
  **Primary arm = cisplatin (CISP)** — the strongest control separation for ATM (KS D=0.41; **OLAP is WEAK
  for ATM**, nonsense vs syn KS p=0.09 — a documented finding: ATM-LOF barely modulates olaparib in this
  non-transformed line, unlike cisplatin/doxorubicin). DOX also separates (p=2e-3); CPT does not.
- **BE-QA (in-run, negated):** Nonsense vs No Mutation **KS D=0.43, p=1.1e-6; MWU p=1.8e-6** → **PASSED.**
  Input `ATM_CISP.tsv` (351 Missense / 91 Silent / 45 Nonsense / 4 Splice-acc / 193 No-Mutation).

## 1. Sweep + the runtime/RAM audit (STRUGGLE #2, measured)
Monomer, 7NI5-A, r4 / **r6 (PRIMARY)** / r8, nRandom **500** (chosen to bound runtime on 2791 residues),
mean/mean. All EXITCODE 0. `C:/Temp/p7_out/ATM_r{4,6,8}/`.

| run | wall time | **peak RAM (sampled)** | modeled | scored |
|---|---|---|---|---|
| ATM r4 | 101 s | **~0.71 GB** | 2791 | 575 |
| ATM r6 | 114 s | ~0.7 GB (≈r8) | 2791 | 1018 |
| ATM r8 | 147 s | **~0.72 GB** | 2791 | 1350 |
(vs CHEK2 543-aa r6 = 37 s.) **The feared O(N²) blowup did NOT materialize.** Reason (code-grounded,
`beclust3d/lfc3d/structure_helpers.py`): the residue neighbor search's **primary path uses Biopython
`NeighborSearch` (a KD-tree)** — O(N log N) — so 2791 residues / 22.5k atoms finish in seconds; the
**pure-Python brute-force O(N_atoms²) double loop** is only a **fallback used if Biopython import fails**
(it *would* be catastrophic on a 22.5k-atom structure, but did not run). So the real ATM bottlenecks are
**upstream (no AF model)** and **statistical (weak signal)**, not compute.

## 2. A0 DISCRIMINATION BENCHMARK — neg (LOF) channel, r6
Universe = 2791 modeled residues; "scored" = 1018. Gold standard = **AlphaMissense** (ATM full-19 table is
unavailable — Q13315 absent from AF-DB — so the **SNV-accessible subset** via dbNSFP/MyVariant, ~5.85
subs/residue, 3055/3056 covered, flagged as a limitation; frac mean-AM>0.564 = **17.6%**) + a **cited
catalytic set** (K2717 ATP-Lys, D2870/N2875 catalytic DxxxxN, D2889 DFG, R3008 PRD; Warren&Pavletich 2018
*Genes Dev*; Baretić 2017; Canman 1998 *Science*).

| Test | Number (p<.05 / p<.001) | Honest reading |
|---|---|---|
| **Base rate** (sig / scored) | **24.6% / 16.4%** | high (of modeled: 8.8% / 5.9%) |
| base rate vs radius (of scored) | r4 27% / r6 24.6% / r8 24.1% | radius-robust |
| **Enrichment** for AM-pathogenic (M=535) | R=66, E=47, **R/E=1.39, OR=1.61, p=1.5e-3** (p<.001: OR 1.98, p 1.3e-4) | **weak** — significant only because N is huge |
| **Discrimination gap** (AM-high vs AM-low) | **12.3% vs 8.4% (+0.04)** | **near-flat** — poor residue-level separation |
| **Precision@10 / @20** (frac AM-pathogenic) | **0.20 (lift 1.04) / 0.40 (lift 2.09)** | **top-10 no better than chance**; @20 recovers as it hits the ATP cleft |
| **domain hit-rate** HEAT / FAT / kinase / FATC | **8.7% / 5.0% / 19.2% / 24.2%** | **real DOMAIN-level concentration on the kinase (+FATC)** — 2–4× the solenoid |
| confidence: hit-rate high-B (flexible) vs ordered | **9.0% vs 8.8%** | flat — no B-factor discrimination |
| **Beat burial baseline** (top-K Naa_count) | burial enr **1.73 > BE3D 1.39** | naive burial edges BE3D |
| **Beat AlphaMissense** (curated catalytic, n=8) | **BE3D 1/8 vs AM 8/8** | **AM dominates** on the catalytic core |

**Honest verdict.** ATM is the case BE3D handles **poorly at residue resolution**. It shows **real
domain-level intelligence** — signal concentrates on the folded **kinase domain (19%) and FATC (24%)** vs
the HEAT solenoid (8.7%) / FAT (5%) — and it does **recover the ATP-cleft neighbourhood** (top cluster
2707/2708/2716–2720/2763/2764, around **K2717**, the one directly-hit catalytic residue). **But** residue-level
discrimination is **near chance**: AM-enrichment 1.4× (gap +0.04, precision@10 no better than chance), it
**loses to both a burial baseline (1.73) and AlphaMissense (8/8 vs 1/8 on catalytic)**, and **57% of its
robust hotspots are likely false positives**. → On a huge, tolerant-surface protein BE3D is at best a
**coarse domain-level flag**, not a residue-level VUS tool.

## 3. FP / FN — including the IDR/over-call audit
- **STRUGGLE #3 — the false-positive pattern (audited).** 230 robust hotspots (sig in ≥2/3 radii):
  **1 TP-functional (K2717), 44 likely-functional (AM-high, kinase), 53 novel, 132 LIKELY-FP** = **57% FP**.
  The FPs are **NOT an IDR explosion** (BE3D does not preferentially call the flexible regions — hit-rate is
  flat vs B-factor). They are **AM-benign, poorly-packed HEAT-solenoid *surface* residues.** Exemplar: the
  **#1–3 robust cluster H600/S601/N602** — all 3 radii, yet **mean-AM 0.09–0.16 (benign)**, **B-factor 86–97
  (flexible)**, **Naa 2–4 (surface, poorly packed)**. A **single depleted guide bleeding over a sparse
  surface neighbourhood** manufactures a "cluster." 109/132 FPs are AM-benign; 23 are high-B/flexible.
- **STRUGGLE #4 — reachability wipes out the catalytic core.** **7 of 8** curated catalytic residues
  (D2870, N2875, D2889, G2891, R2871, H2872, R3008) are **base-editing-UNREACHABLE** (a CBE cannot install
  the catalytic Asp/Asn/Gly codon changes) → **unscored, invisible**. Only **K2717** (ATP-Lys) is both
  reachable and flagged (TP). AlphaMissense scores all 8 at ~1.0 with no experiment. This is the assay's
  ceiling, but it means **BE3D structurally cannot see ATM's catalytic determinants**.
- Files: `outreach/ATM/ATM_hotspots.tsv` (230 rows + flags/mean-AM/ClinVar/B-factor), `ATM_G2P.tsv`.

## 4. WHERE BE3D STRUGGLED (ATM) — summary
1. **No AF model at all** (>2700 aa excluded) → auto-fetch 404/`{}`; must hand-supply a cryo-EM PDB, which
   is a **dimer with 9% unresolved** and whose **B-factor is mislabeled "pLDDT"** by BE3D (inverts the
   confidence semantics of the characterization step).
2. **Weak residue-level discrimination** (OR 1.6, gap +0.04, precision@10 ≈ chance); **loses to burial and
   AlphaMissense**; only a **coarse domain-level** flag survives.
3. **57% of robust hotspots are likely false positives** — AM-benign, flexible, poorly-packed solenoid
   surface (e.g. 600–602). Spatial aggregation over a sparse surface fabricates clusters.
4. **Reachability** removes 7/8 catalytic residues → the functional core is largely invisible.
5. (**Non-issue, tested:** runtime/RAM — 100–147 s, ~0.7 GB; the KD-tree path scales fine.)

## 5. Cited change proposals (patient-impact-ordered)
1. **[highest value] pLDDT/AM confidence gate + residue-level FDR + a packing/occupancy check on hotspot
   calls.** 57% of ATM hotspots are AM-benign flexible surface; the 600–602 FP is high-B + Naa≤4. Gating on
   confidence + local packing (Naa_count) + BH-q would remove most. Cited by §2 base-rate + §3 FP audit.
   *Partly covered by open PRs #20 (FDR) and #23 (shortlist w/ pLDDT); a **local-packing/occupancy** term is
   the new ask (roadmap issue #24 "occupancy-model").*
2. **Handle proteins with no AF model / >2700 aa** — detect the 404/`{}`, and when a `user_pdb` is
   experimental, **stop labeling B-factor as pLDDT** (read `_processed` header or a `structure_source` flag)
   so confidence-aware steps aren't inverted. Cited by §0/§4.1.
3. **Base-editing reachability report** (7/8 ATM catalytic residues unreachable). *Covered by PR #19.*
4. **Domain-restricted analysis mode** — on a 3056-aa protein the honest signal is domain-level; a mode that
   scopes randomization/clustering to a chosen folded domain (kinase 2686–2998) would sharpen the null and
   the FDR vs diluting them across 2791 residues. Cited by §2 domain-hit-rate contrast.
