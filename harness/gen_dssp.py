"""
Generate a minimal but VALID classic-DSSP-format file for a processed PDB so that
DSSPparser.parseDSSP (used by beclust3d structure_helpers.parse_dssp) can parse it
WITHOUT the mkdssp binary.

The DSSPparser reads fixed-width columns (0-indexed slices):
    [0:5]   resnum (sequential)
    [5:10]  inscode  <- must equal the uniprot/PDB residue number (matched in parse_dssp)
    [10:12] chain
    [12:14] aa (1-letter)
    [14:17] struct (SS9)   -> blank/coil is fine (only used in final characterization)
    [34:38] acc (numeric)  -> must be float-parseable
    [103:109] phi (numeric)
    [109:115] psi (numeric)
Only chain == target_chainid rows are used downstream. SS codes are placeholder coil
because DSSP secondary-structure only affects the final characterization step.
"""
import sys
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa
from Bio.SeqUtils import seq1

pdb_path = sys.argv[1]
out_path = sys.argv[2]

parser = PDBParser(QUIET=True)
structure = parser.get_structure("s", pdb_path)
model = structure[0]

def build_line(seq, resnum, chain, aa1, acc, phi, psi):
    line = [" "] * 140
    def put(s, start, end, right=True):
        s = str(s)
        w = end - start
        s = s[:w]
        s = s.rjust(w) if right else s.ljust(w)
        for i, ch in enumerate(s):
            line[start + i] = ch
    put(seq, 0, 5)
    put(resnum, 5, 10)
    put(chain, 10, 12)          # e.g. " B"
    put(aa1, 12, 14)            # e.g. " M"
    # struct [14:17] left blank -> coil placeholder
    put(acc, 34, 38)           # numeric ASA placeholder
    put(f"{phi:.1f}", 103, 109)
    put(f"{psi:.1f}", 109, 115)
    return "".join(line).rstrip() + "\n"

lines = []
# classic-DSSP header line containing '#': parser skips until it sees '#'
lines.append("  #  RESIDUE AA STRUCTURE BP1 BP2  ACC     N-H-->O    O-->H-N    N-H-->O    O-->H-N    TCO  KAPPA ALPHA  PHI   PSI    X-CA   Y-CA   Z-CA\n")

seq = 0
nres = 0
for chain in model:
    cid = chain.id
    for res in chain:
        if not is_aa(res, standard=True):
            continue
        resnum = res.id[1]
        try:
            aa1 = seq1(res.get_resname())
        except Exception:
            aa1 = "X"
        if aa1 == "" or aa1 == "X":
            continue
        seq += 1
        nres += 1
        # placeholder solvent accessibility + torsions (coil-ish)
        lines.append(build_line(seq, resnum, cid, aa1, 100, -60.0, -45.0))

with open(out_path, "w", newline="\n") as fh:
    fh.writelines(lines)

print(f"Wrote {nres} residues to {out_path}")
