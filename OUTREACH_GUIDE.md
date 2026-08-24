# Outreach & Benchmarking Guide (Phase 3) — shared spec for all per-user packages

Goal: for each (user persona × real protein), turn a real BE3D run into a shareable package that a
working scientist would actually value: benchmarked against ground truth, compared to competitors,
with a copyable email + downloadable files. Be HONEST and specific; never overclaim.

## A0. DISCRIMINATION TEST — the primary, non-negotiable benchmark (do this FIRST)
"BE3D recovered known hotspot X" is MEANINGLESS if BE3D flags a large fraction of the protein — it could
just be calling everything a hotspot. You MUST test discrimination, not recall:
1. **Base rate.** Report the fraction of residues BE3D calls significant at each p-threshold
   (#sig / #scored) AND the fraction that end up in a cluster. If the base rate is high (say >15–20%),
   any overlap with known sites is weak evidence — SAY SO plainly.
2. **Enrichment vs chance (the real test).** Curate, INDEPENDENTLY of BE3D, two residue sets:
   - `FUNCTIONAL` = experimentally-validated functional residues (active/catalytic site, binding/PPI
     interface, ClinVar-pathogenic clusters, validated resistance/hotspot residues) — with citations.
   - `TOLERANT` = residues known to tolerate mutation (ClinVar-benign, silent/neutral positions, exposed
     non-functional loops).
   Build a 2×2 (BE3D-hotspot vs not) × (FUNCTIONAL vs not) and run **Fisher's exact / hypergeometric**:
   under the null of drawing K hotspots at random from N residues, expected overlap E = K·M/N (M=#FUNCTIONAL);
   report observed overlap R, **enrichment R/E, odds ratio, and p-value**. BE3D is only "intelligent" here
   if R/E is clearly >1 with a significant p — quantify it, don't assert it.
3. **Discrimination gap.** Report BE3D's hit-rate on FUNCTIONAL residues vs on TOLERANT residues. A real
   signal means functional ≫ tolerant. If BE3D flags TOLERANT residues at a similar rate → it is NOT
   discriminating (the "everything is a hotspot" failure) — report this honestly.
4. **Precision@K (rank).** Take the top-K residues by |LFC3D| (K = 10, 20). What fraction sit on/near
   FUNCTIONAL sites vs the base expectation M/N? The strongest signals should concentrate more than chance.
5. **Beat a trivial baseline.** Compare BE3D's functional-site enrichment to a naive baseline on the SAME
   protein: e.g. "call all buried (low-RSA) residues" or "call all AlphaMissense-high residues". Does BE3D
   add discrimination BEYOND these trivial predictors, or merely recapitulate burial/conservation?
Write these numbers up front. The credibility of everything below depends on BE3D beating chance here.

## A. Ground-truth benchmark (false-positive / false-negative flagging) — REQUIRED
- Compile the EXPERIMENTALLY VALIDATED functional hotspots for this protein from literature + databases
  (cite: PubMed/PMID, UniProt features, ClinVar, cBioPortal/3dhotspots.org, COSMIC, specific papers).
- Classify each ROBUST BE3D hotspot residue (persistent across the sweep) as one of:
  - **TP (validated)** — matches a known experimentally-validated functional/interface/pocket residue.
  - **NOVEL-CANDIDATE** — no prior validation, but structurally plausible (state why); a hypothesis, not a claim.
  - **LIKELY FALSE POSITIVE** — flag explicitly: strong BE3D signal with no known function AND a benign/tolerant
    context (e.g. surface loop, low pLDDT, disordered, ClinVar-benign region). Say WHY it may be spurious
    (e.g. single high-LFC guide dominating a neighborhood; low-pLDDT geometry; bystander-edit artifact).
- Also list **FALSE NEGATIVES**: validated hotspots BE3D MISSED, and WHY (e.g. base-editing-inaccessible —
  the required codon change isn't reachable by CBE/ABE, as with DNMT3A R882H; or no guide coverage).
- Give an honest precision/recall FEEL (not a fake exact number): "of N robust hotspots, X match known
  functional sites, Y are novel candidates, Z look like false positives."

## B. Competitor comparison — REQUIRED
Compare BE3D's hotspots on THIS protein against, where feasible, actually pulling competitor output:
- **Somatic 3D-hotspot tools** (3dhotspots.org / cBioPortal 3D hotspots / HotMAPS-style) — these use
  patient mutation RECURRENCE, not screen LFC. Fetch this protein's precomputed 3D hotspots and compare.
- **AlphaMissense** per-residue pathogenicity (via Ensembl/AlphaMissense lookup for key residues) — a
  sequence/structure PREDICTOR, orthogonal. Do BE3D hotspots have high AlphaMissense? Where do they differ?
- **ProTiler-Mut** (Cell Systems 2026) — the one direct tiling→3D rival; compare conceptually (3D-RRA vs
  BE3D's randomization null + meta), note what each would/wouldn't catch.
- Produce a small table: tool | input signal | what it finds on THIS protein | advantage | disadvantage.
- Be fair: state clearly where a competitor is BETTER (e.g. somatic tools have decade-long benchmarking;
  AlphaMissense needs no screen) and where BE3D wins (functional screen signal, controls, allosteric/
  resistance sites recurrence-based tools can't see, editor-specific escape geometry).

## C. Copyable outreach email — REQUIRED (outreach/<GENE>/<GENE>_email.md)
Tone: a knowledgeable peer sharing a genuinely useful, hypothesis-generating result. NOT a sales pitch.
Structure:
1. One-line who/why: "I ran your (or the public <paper>) base-editing screen of <GENE> through BE3D, a
   structure-function tool from the Iqbal Lab (Broad) + Liau Lab (Harvard)."
2. Credibility with HONEST track record: BE3D has recovered known functional hotspots zero-shot before —
   e.g. it re-found KBTBD4's medulloblastoma R379–R390 Kelch site from screen signal alone, and on YOUR
   protein it independently recovered [the validated positive controls found in this run]. Cite them.
3. The result relevant to them: the 1–3 most decision-relevant hotspots (validated + the novel candidate),
   in plain language tied to their biology. Include the explicit caveats (direction, pLDDT, false-positive flags).
4. How to apply it: the concrete next step (residues to mutate/validate; how to run it themselves — the
   Colab/GitHub link; the G2P interactive view; attached files).
5. Optional collaboration note: "If useful, we'd value your validation results back — experimental
   confirmation (or refutation) of these hotspots helps us tune BE3D's thresholds and null model."
6. Attachments list (the downloadable files).
Keep it ~200–300 words, skimmable, honest about uncertainty. Also write a 1-page brief
(outreach/<GENE>/<GENE>_brief.md) they can forward.

## D. Downloadable files — REQUIRED (assemble into outreach/<GENE>/)
- `<GENE>_hotspots.tsv` — the robust hotspot table (residue, AA, LFC3D, direction, significance, TP/novel/FP flag).
- the G2P-formatted TSV from the run (interactive 3D view input).
- `<GENE>_brief.md` (1-page) + `<GENE>_email.md`.
- Keep it lean; no giant randomization dumps.

## Honesty rules (non-negotiable)
- Never claim a novel hotspot is validated. Label hypotheses as hypotheses.
- Always show the false-positive flags and false negatives — that candor is the credibility.
- The "zero-shot track record" claim must reference REAL recoveries (KBTBD4 MB site; this run's positive controls).
- Don't imply endorsement by the BE3D authors; this is an independent analysis.
