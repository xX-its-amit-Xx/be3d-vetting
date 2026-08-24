# BE3D Vetter — Master Spec (re-anchor here on every wake-up)

## Meta-goal
Comprehensively "vet" **BE3D** (broadinstitute/BE3D, pip pkg `beclust3d` v1.0.0) — a Python package
for structure-function analysis of base-editor (BE) tiling-mutagenesis screens — for **every possible
use case and user**, and produce a high-value comparison/evaluation dashboard + deep-dive artifacts.

Authors: Calvin XiaoYang Hu (Harvard, xiaohu@g.harvard.edu), Yoochan Myung (Broad). Ships with the
Genomics 2 Proteins Portal (G2P, g2p.broadinstitute.org). ("APAR lab" per user — confirm exact lab.)

## Orchestration model (user-requested)
- I am the **meta-king** orchestrator. I delegate to **queens** (subagents), who may delegate to **workers**.
- I MUST be resilient to usage limits / context resets: on every wake-up, re-read this SPEC.md + STATE.md
  and continue the next unfinished item. Never lose the thread.
- Keep scheduling wake-ups to check on delegated work until all deliverables are DONE.

## Deliverables (the 11 asks, consolidated)
- D1. Comparison/evaluation **dashboard (single HTML)** ranking/comparing BE3D vs equivalent tools across
      exhaustive use cases → **deploy to Vercel**.
- D2. **Use-case deep dives** showing BE3D's *actual output*, key scientist decisions, deep analysis.
- D3. **Design deep-dive** of the BE3D repo (architecture, modules, data model, choices, smells).
- D4. **Exhaustive personas & use cases** — academic groups + industry groups + everyone who'd use BE3D;
      decision criteria (BE3D vs alternatives); **quantify impact / market**.
- D5. **Cross-domain analogies** (finance, healthcare, art, etc.).
- D6. **Critique + compliment + double-down-on-niche + improvements**; adoption into scientist workflow
      (LLM, MCP, what tools BE3D serves).

## What BE3D is (ground truth — from repo)
3 modules: **BE-QA** (Mann-Whitney U + KS tests: knockout vs neutral LFC distributions, 2 hypotheses),
**BE-Clust3D** (map LFC by residue onto AlphaFold/PDB structure → LFC3D spatial score, default 6Å radius →
randomization baseline → z-norm/p-values → agglomerative clustering of hotspots → characterization:
enrichment tests, pLDDT/RSA scatter, feature barplots), **BE-MetaClust3D** (aggregate multiple screens /
cross-species / cross-isoform via MUSCLE/CLUSTAL alignment). Output → G2P-compatible TSV.
Inputs: BE screen scores TSV (Mutation_type, Mutation_list, Gene, sgRNA_score) + UniProt ID (auto AlphaFold)
or custom FASTA/PDB. Run via Colab or local (YAML config). ~6100 LOC, MIT license.
Full function list: see BE3D_BRIEF.md and the cloned repo README_FUNCTIONS.md.

## PHASE 2 MANDATE (user, 2026-08-23) — USE BE3D on REAL, NEW problems
The Phase-1 KBTBD4 run was real but only reproduced the shipped example. Phase 2:
- **Queens orchestrate workers** to ACTUALLY USE BE3D to solve real problems / real gaps —
  NOT reproducing the source papers' findings. Propose genuinely NEW, actionable next steps for other groups.
- Workers must: (1) RESEARCH real public base-editor tiling-mutagenesis datasets suitable for BE3D
  (NOT KBTBD4/HDAC1/MORC2), reformat to BE3D input; (2) run BE3D EXHAUSTIVELY (parameter sweeps);
  (3) keep a CLEAR DECISION LOG (every parameter: radius, aggregation, p-thr, nRandom, structure choice,
  QA, meta strategy, residue vs atom) and data provenance; (4) extract NOVEL hotspots/proposals.
- Throughout, LOG improvements to BE3D in BE3D_IMPROVEMENTS.md as real usage surfaces gaps.
- Deliverable: new-findings report(s) + decision logs + dashboard update + redeploy.

## Reusable harness (Phase-1, keep using)
- venv: scratchpad/be3dvenv (py3.14; beclust3d editable-installed from scratchpad/BE3D).
- Patches applied to CLONE source (inspectable): DSSP-optional, pandas sep=, filename sanitize, find_union.
- scratchpad/gen_dssp.py: generate a placeholder classic-DSSP from a PDB (bypasses missing mkdssp).
- Run: cd scratchpad/BE3D/examples && ../../be3dvenv/Scripts/python.exe be3d_local.py <yaml>
- Deploy: copy site/* → scratchpad/be3d-vetting/ ; `npx vercel deploy --prod --yes` (authed).

## Status
See STATE.md for the live checklist.
