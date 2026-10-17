# How much of each key residue is exposed on the molecule's surface.
from Bio.PDB import MMCIFParser
from Bio.PDB.DSSP import residue_max_acc
from Bio.PDB.SASA import ShrakeRupley

LARGEST = residue_max_acc[
    "Sander"
]  # each amino acid's area fully exposed
RESIDUES = (6, 42, 63, 85, 88, 92)

for code in ("2HHB", "2HBS"):
    model = MMCIFParser(QUIET=True).get_structure(
        code, f"../data/{code}.cif"
    )[0]
    for chain in list(model):  # keep one hemoglobin molecule, A to D
        if chain.id not in "ABCD":
            model.detach_child(chain.id)
    ShrakeRupley().compute(model, level="R")
    print(code)
    for number in RESIDUES:
        residue = model["B"][number]
        share = residue.sasa / LARGEST[residue.get_resname()]
        print(
            f"  {residue.get_resname()} {number:3}  {residue.sasa:5.0f}"
            f" square angstroms, {share:4.0%} exposed"
        )
