# KRAS(G12C) 3D resistance atlas — sotorasib vs adagrasib (BE3D, 1-page brief)

**What / source.** Independent BE3D analysis of the public KRAS(G12C) base-editing drug-resistance
screen (Coelho, Dincer et al. 2024, *Nat Genet*). BE3D (Iqbal Lab, Broad; Liau Lab, Harvard) projects
per-guide resistance log-fold-changes onto AlphaFold KRAS (P01116, res 1–189), builds a per-residue
randomization null, and 3D-clusters the significant residues. Positive direction = enrichment = drug
RESISTANCE. Not endorsed by the BE3D authors.

**Headline.** The two approved G12C covalent inhibitors have **quantitatively different 3D resistance
signatures**, and BE3D reads inhibitor-*escape* geometry that recurrence/conservation tools cannot.

## Benchmark (honest, done before any claims)
Curated the validated G12Ci resistance residues and the switch-II pocket *independently* of BE3D, then
tested enrichment (hypergeometric, N=189 residues):

| | Adagrasib | Sotorasib |
|---|---|---|
| Base rate (p<0.05, % of protein) | 15.9% | 12.7% |
| Enrichment vs 9 validated resistance residues | **R/E 5.6×, OR 57, p=1.4×10⁻⁶** | 2.6×, p=0.09 (NS) |
| Enrichment vs switch-II pocket (60 res) | **2.8×, OR 34, p=5×10⁻¹³** | 1.05×, p=0.51 (NS) |
| Functional vs tolerant hit-rate (discrimination gap) | **45% vs 1.9% (~24×)** | 13% vs 9% (1.4×) |
| Precision@10 (top |LFC3D|) | **1.00** | 0.40 |
| Validated resistance residues recovered zero-shot | **8 / 9** | 3 / 9 (pocket) |

Beats trivial baselines: "call all buried residues" fails on the resistance residues (R/E 0.6×);
AlphaMissense is drug-agnostic and can't separate resistance from oncogenic activation.
**Verdict:** the **adagrasib map is benchmarked and strong**; the **sotorasib arm does not enrich the
known pocket** — its signal is distal (see α5 candidate), so we present that as a hypothesis, not a
recovery.

## The two findings
1. **Different escape maps.** Adagrasib = broad orthosteric shell (switch-II 60–72 + α3 groove 89–103,
   incl. His95/Tyr96/Gln99); sotorasib = sparse pocket footprint. Shared hits only {60,65,66,67,68};
   Jaccard ≈ 0.12. Consistent with pharmacology: adagrasib H-bonds Tyr96 and engages the His95/Q99
   groove; sotorasib contacts Tyr96 only via water and is His95-independent.
   → **Design lever:** a sotorasib-class binder extended toward the His95/Q99 exit should have a higher
   genetic barrier to resistance; adagrasib's dense pocket dependence makes single-residue escapes
   (Y96D, R68, H95, Q99) efficient → pair with SOS1/SHP2 or a distinct-pocket agent.
2. **Sotorasib-specific distal candidate, α5/SAK 147–156** (K147–F156; ~13–27 Å from drug; well-
   structured, pLDDT 97–99; 148/154 robust across all runs). No prior resistance annotation except
   F156L (Feng 2022, mechanism unassigned). **Hypothesis:** SAK/α5 perturbs nucleotide loading /
   switch-II equilibrium that sotorasib (an inactive-state trap) depends on. **To validate**, not a claim.

## Flags (in KRAS_hotspots.tsv)
- **TP:** T58, S65, R68, D69, M72, D92, H95, Y96, Q99 (validated; Feng 2022 PMID 35471904; Awad 2021
  PMID 34161704; Tanaka 2021 PMID 33846219).
- **NOVEL-CANDIDATE:** α5/SAK 147–156 (sotorasib).
- **LIKELY FALSE POSITIVE:** HVR 178–180 (pLDDT 45–49, disordered — artifact).
- **UNCERTAIN:** sotorasib 46–48 loop, 142; adagrasib 90, 127.
- **False negatives:** adagrasib missed T58 (edit accessibility); sotorasib's absence of H95/Y96/Q99 is
  mechanistically correct (adagrasib-specific), not a miss.

## How to use
- `KRAS_hotspots.tsv` — per-residue table with flags. `KRAS_{sotorasib,adagrasib}_g2p.tsv` — interactive
  3D (G2P) view input.
- Test α5: mutate T148/R149/V152/D154 in KRAS(G12C), measure sotorasib vs adagrasib IC50 + GTP-loading;
  Boltz co-fold as in-silico pre-screen. Run BE3D on your own screens via the Colab/GitHub release.

Competitors: 3dhotspots.org / cBioPortal (Gao 2017 PMID 28115009) report KRAS G12/G13/Q61/A146 — those
are oncogenic **activation** codons from tumor recurrence, an axis orthogonal to inhibitor resistance.
BE3D is the only one of these that maps where the drug's grip fails, per drug.
