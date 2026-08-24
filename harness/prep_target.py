#!/usr/bin/env python
"""
prep_target.py -- plug-and-play AlphaFold + DSSP fetcher for BE3D.

Given a UniProt accession, downloads the AlphaFold predicted structure and
generates a matching classic-DSSP file (so the BE3D pipeline never has to shell
out to the missing `mkdssp` binary). Both artifacts land in --outdir and their
paths are printed (machine-parseable KEY=VALUE lines + human summary).

The generated DSSP is a placeholder (coil SS, filler ASA/torsions) that is only
consumed by BE3D's final *characterization* step -- it does NOT affect
LFC3D / significance / clustering. See gen_dssp.py for the column layout.

Example:
    python prep_target.py --uniprot P61626 --gene LYZ --outdir ./targets/LYZ

Notes on AlphaFold accessions:
    AlphaFold DB only hosts the canonical sequence (fragment F1). Pass the BASE
    accession (e.g. P61626, or Q9NVX7 -- NOT the isoform Q9NVX7-2). If you need a
    specific isoform's numbering you must supply your own PDB to run_new_target.py.
"""
import argparse
import os
import subprocess
import sys
import urllib.request
import urllib.error

# DSSP generator lives next to this script; reused verbatim via subprocess so
# the two stay in lockstep (single source of truth for the DSSP column layout).
GEN_DSSP = os.path.join(os.path.dirname(os.path.abspath(__file__)), "gen_dssp.py")

# AlphaFold model versions to try, newest first. v6 is current (2024+), v4 is the
# long-standing default, v2 is a fallback for older/rarely-updated entries.
AF_VERSIONS = ("v6", "v4", "v2")
AF_URL = "https://alphafold.ebi.ac.uk/files/AF-{uniprot}-F1-model_{ver}.pdb"


def download_alphafold(uniprot, outdir):
    """Try v6 -> v4 -> v2; return the local path of the first that downloads."""
    last_err = None
    for ver in AF_VERSIONS:
        url = AF_URL.format(uniprot=uniprot, ver=ver)
        dest = os.path.join(outdir, f"AF-{uniprot}-F1-model_{ver}.pdb")
        try:
            print(f"[prep_target] trying {url}")
            req = urllib.request.Request(url, headers={"User-Agent": "BE3D-prep/1.0"})
            with urllib.request.urlopen(req, timeout=60) as resp:
                data = resp.read()
            if not data or not data.lstrip().startswith((b"HEADER", b"ATOM", b"CRYST", b"TITLE", b"REMARK")):
                raise ValueError("downloaded content does not look like a PDB file")
            with open(dest, "wb") as fh:
                fh.write(data)
            print(f"[prep_target] downloaded {ver}: {dest} ({len(data)} bytes)")
            return dest, ver
        except (urllib.error.HTTPError, urllib.error.URLError, ValueError, TimeoutError) as e:
            print(f"[prep_target]   {ver} failed: {e}")
            last_err = e
    raise SystemExit(
        f"[prep_target] ERROR: could not download AlphaFold model for '{uniprot}' "
        f"(tried {', '.join(AF_VERSIONS)}). Last error: {last_err}. "
        f"Confirm the BASE accession resolves at https://alphafold.ebi.ac.uk/entry/{uniprot}"
    )


def generate_dssp(pdb_path, dssp_path):
    """Reuse gen_dssp.py (same interpreter) to build a parseable classic-DSSP."""
    cmd = [sys.executable, GEN_DSSP, pdb_path, dssp_path]
    print(f"[prep_target] generating DSSP: {' '.join(cmd)}")
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.stdout:
        print(res.stdout.rstrip())
    if res.returncode != 0:
        print(res.stderr, file=sys.stderr)
        raise SystemExit(f"[prep_target] ERROR: gen_dssp.py failed (exit {res.returncode})")
    return dssp_path


def main():
    ap = argparse.ArgumentParser(description="Download AlphaFold model + generate DSSP for a BE3D target.")
    ap.add_argument("--uniprot", required=True, help="BASE UniProt accession (e.g. P61626). NOT an isoform id.")
    ap.add_argument("--outdir", required=True, help="Directory to write the PDB + DSSP into.")
    ap.add_argument("--gene", default=None, help="Optional gene symbol, used only in the printed summary.")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    pdb_path, ver = download_alphafold(args.uniprot, args.outdir)
    dssp_path = os.path.join(args.outdir, f"AF-{args.uniprot}-F1-model_{ver}.dssp")
    generate_dssp(pdb_path, dssp_path)

    pdb_path = os.path.abspath(pdb_path)
    dssp_path = os.path.abspath(dssp_path)

    print("\n[prep_target] DONE" + (f" for {args.gene}" if args.gene else ""))
    # Machine-parseable lines (run_new_target.py / scripts can grep these):
    print(f"PDB_PATH={pdb_path}")
    print(f"DSSP_PATH={dssp_path}")
    print(f"AF_VERSION={ver}")


if __name__ == "__main__":
    main()
