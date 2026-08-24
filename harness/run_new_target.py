#!/usr/bin/env python
"""
run_new_target.py -- thin wrapper that runs BE3D (monomer mode) on ANY target.

Given a reformatted screen TSV (or several, comma-separated, for meta-aggregation),
the four screen column names, a gene/UniProt, and the AlphaFold PDB + DSSP produced
by prep_target.py, this:

  1. WRITES a self-contained BE3D monomer YAML (conservation OFF, user_pdb/user_dssp
     wired in, nRandom / radii / QA cases+controls configurable), and
  2. INVOKES examples/be3d_local.py on that YAML with this interpreter (the BE3D venv).

The pipeline runs LFC -> LFC3D -> significance -> clustering -> characterization and
writes LFC/, LFC3D/, cluster_LFC/, cluster_LFC3D/, etc. under --outdir.

Screen TSV format expected (tab-separated), with column names supplied on the CLI:
  * gene_col      : gene symbol column; rows are filtered to --gene
  * mut_col       : per-guide mutation category (Missense / Silent / Nonsense / No Mutation ...)
  * mut_list_col  : edit list column, e.g. "M1V;K2G;" (refAA + pos + altAA, ';'-joined)
  * val_col       : numeric guide score (log fold-change / sgRNA score)

Example (single screen):
  python run_new_target.py \
      --screen ./plumbing_test_LYZ.tsv \
      --gene LYZ --uniprot P61626 \
      --pdb ./targets/LYZ/AF-P61626-F1-model_v4.pdb \
      --dssp ./targets/LYZ/AF-P61626-F1-model_v4.dssp \
      --outdir ./output/LYZ_smoke \
      --mut-list-col Mutation_list --mut-col Mutation_type \
      --val-col sgRNA_score --gene-col Gene \
      --nrandom 100

Multiple screens (meta-aggregation): --screen a.tsv,b.tsv
"""
import argparse
import os
import shutil
import subprocess
import sys

import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
BE3D_ROOT = os.path.join(HERE, "BE3D")
BE3D_LOCAL = os.path.join(BE3D_ROOT, "examples", "be3d_local.py")


def build_config(args, screen_dir, screen_basenames):
    """Assemble the monomer YAML dict. Kept explicit so it is easy to audit/tweak."""
    cfg = {
        "mode": "monomer",
        # beclust3d is editable-installed in the venv; this path only feeds a
        # harmless sys.path.append in be3d_local.py. Point it at the clone root.
        "beclust3d_path": os.path.abspath(BE3D_ROOT),
        "screen_dir": screen_dir,
        "output_dir": os.path.abspath(args.outdir),
        "muscle_path": "muscle",  # unused: conservation is OFF
        "input_gene": args.gene,
        "input_uniprot": args.uniprot,
        "input_chain": args.chain,
        "screens": ",".join(screen_basenames),
        "user_pdb": os.path.abspath(args.pdb),
        "user_dssp": os.path.abspath(args.dssp),
        # monomer mode ignores PPI config, but keep the keys present/None-safe
        "ppi_chain_gene_dict": None,
        "ppi_gene_edits_dict": {},
        "mutation_category": {
            "missense": ["Missense"],
            "silent": ["Silent"],
            "nonsense": ["Nonsense"],
            "no_mutation": ["No Mutation"],
            "splice": ["Splice-donor", "Splice-acceptor"],
            "intron": ["Intron"],
        },
        # collapse a guide's ';'-joined multi-category edits to one category
        "mutation_priority": [
            "Nonsense", "Splice-donor", "Splice-acceptor",
            "Missense", "Silent", "Intron", "No Mutation",
        ],
        "qa": {
            "qa_passed_only": False,
            "qa_only": False,
            "cases": args.qa_cases,
            "controls": args.qa_controls,
        },
        "conservation": {
            "run": False,
            "v_score_threshold": 3,
            "alt_gene_name": None,
            "alt_uniprot_id": None,
            "alt_screen_start": None,
        },
        "priority_on_alternative": False,
        "database": {
            "mut_list_col": None,
            "mut_col": args.mut_col,
            "val_col": args.val_col,
            "gene_col": args.gene_col,
            "edits_col": args.mut_list_col,
            "mut_delimiter": args.mut_delimiter,
            "gRNA_col": None,
        },
        "user_fasta": os.path.abspath(args.fasta) if args.fasta else None,
        "function_for_lfc": "mean",
        "function_for_lfc3d": args.function_for_lfc3d,
        "function_for_meta": "SUM",
        "nRandom": args.nrandom,
        "pthr": {
            "single_screen": args.single_screen_pthr,
            "multi_screen": args.multi_screen_pthr,
        },
        "structure_radius": args.structure_radius,
        "clustering_radius": args.clustering_radius,
        "atom_level_naa": args.atom_level,
    }
    return cfg


def main():
    ap = argparse.ArgumentParser(
        description="Write a BE3D monomer YAML and run it on an arbitrary target.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    # --- required inputs ---
    ap.add_argument("--screen", required=True,
                    help="Screen TSV path, or comma-separated list for meta-aggregation.")
    ap.add_argument("--gene", required=True, help="Gene symbol (must match values in gene_col).")
    ap.add_argument("--uniprot", required=True, help="UniProt accession (used for structureid).")
    ap.add_argument("--pdb", required=True, help="AlphaFold/user PDB (from prep_target.py).")
    ap.add_argument("--dssp", required=True, help="Classic-DSSP file (from prep_target.py).")
    ap.add_argument("--outdir", required=True, help="Output directory for the BE3D run.")
    # --- column names ---
    ap.add_argument("--mut-list-col", dest="mut_list_col", required=True,
                    help="Edit-list column (BE3D edits_col), e.g. Mutation_list.")
    ap.add_argument("--mut-col", dest="mut_col", required=True,
                    help="Mutation-category column, e.g. Mutation_type.")
    ap.add_argument("--val-col", dest="val_col", required=True,
                    help="Numeric guide-score column, e.g. sgRNA_score.")
    ap.add_argument("--gene-col", dest="gene_col", required=True,
                    help="Gene-symbol column, e.g. Gene.")
    # --- tunables ---
    ap.add_argument("--chain", default="A", help="Target chain id (AF monomer = A).")
    ap.add_argument("--nrandom", type=int, default=500, help="Randomization iterations.")
    ap.add_argument("--structure-radius", dest="structure_radius", type=float, default=6.0)
    ap.add_argument("--clustering-radius", dest="clustering_radius", type=float, default=6.0)
    ap.add_argument("--single-screen-pthr", dest="single_screen_pthr", type=float, default=0.05)
    ap.add_argument("--multi-screen-pthr", dest="multi_screen_pthr", type=float, default=0.05)
    ap.add_argument("--mut-delimiter", dest="mut_delimiter", default=";")
    ap.add_argument("--function-for-lfc3d", dest="function_for_lfc3d", default="mean",
                    choices=["mean", "sum", "median", "min", "max"],
                    help="Neighborhood aggregation for LFC3D.")
    ap.add_argument("--atom-level", dest="atom_level", action="store_true",
                    help="Use atom-level (not residue-level) neighbor search.")
    ap.add_argument("--qa-cases", nargs="+", default=["Nonsense"],
                    help="Mutation categories treated as QA 'cases'.")
    ap.add_argument("--qa-controls", nargs="+", default=["No Mutation"],
                    help="Mutation categories treated as QA 'controls'.")
    ap.add_argument("--fasta", default=None,
                    help="Optional UniProt FASTA; omit to let BE3D fetch it online.")
    ap.add_argument("--dry-run", action="store_true", help="Write the YAML but do not run BE3D.")
    args = ap.parse_args()

    for p in (args.pdb, args.dssp):
        if not os.path.isfile(p):
            raise SystemExit(f"[run_new_target] ERROR: file not found: {p}")

    os.makedirs(args.outdir, exist_ok=True)

    # Stage every screen into <outdir>/screens/ so a single screen_dir works even
    # when the input screens live in different directories.
    staged_dir = os.path.join(os.path.abspath(args.outdir), "screens")
    os.makedirs(staged_dir, exist_ok=True)
    screen_basenames = []
    for s in [x.strip() for x in args.screen.split(",") if x.strip()]:
        if not os.path.isfile(s):
            raise SystemExit(f"[run_new_target] ERROR: screen TSV not found: {s}")
        base = os.path.basename(s)
        shutil.copy2(s, os.path.join(staged_dir, base))
        screen_basenames.append(base)
    print(f"[run_new_target] staged {len(screen_basenames)} screen(s) into {staged_dir}")

    cfg = build_config(args, staged_dir, screen_basenames)
    # IMPORTANT: be3d_local.py does shutil.copy2(config_yaml, output_dir/basename),
    # so the source YAML must NOT live inside output_dir or copy2 hits src==dst
    # (WinError 32 on Windows / SameFileError elsewhere). Keep it in a sibling dir.
    cfg_dir = os.path.join(os.path.dirname(os.path.abspath(args.outdir)), "be3d_configs")
    os.makedirs(cfg_dir, exist_ok=True)
    yaml_path = os.path.join(cfg_dir, f"{args.gene}_run.yaml")
    with open(yaml_path, "w") as fh:
        yaml.safe_dump(cfg, fh, sort_keys=False, default_flow_style=False)
    print(f"[run_new_target] wrote config: {yaml_path}")

    if args.dry_run:
        print("[run_new_target] --dry-run set; not invoking BE3D.")
        return

    cmd = [sys.executable, BE3D_LOCAL, yaml_path]
    print(f"[run_new_target] invoking: {' '.join(cmd)}")
    # be3d_local.py resolves screen_dir/output_dir relative to CWD when not absolute;
    # everything here is absolute, but run from examples/ for parity with the repo.
    res = subprocess.run(cmd, cwd=os.path.join(BE3D_ROOT, "examples"))
    if res.returncode != 0:
        raise SystemExit(f"[run_new_target] BE3D exited with code {res.returncode}")
    print(f"[run_new_target] BE3D completed. Outputs under: {os.path.abspath(args.outdir)}")


if __name__ == "__main__":
    main()
