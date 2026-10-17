# How close is AlphaFold's prediction of beta globin to the crystal?
from Bio.PDB import MMCIFParser, PDBParser, Superimposer

crystal = MMCIFParser(QUIET=True).get_structure(
    "2HHB", "../data/2HHB.cif"
)[0]["B"]
predicted = PDBParser(QUIET=True).get_structure(
    "AF", "../data/AF-P68871-F1.pdb"
)[0]["A"]

# AlphaFold stores its confidence, pLDDT, in the B-factor column
confidence = [residue["CA"].bfactor for residue in predicted]
low = [n for n, c in enumerate(confidence, start=1) if c < 90]
print(
    f"mean pLDDT {sum(confidence) / len(confidence):.1f}, "
    f"below 90 at residues {low}"
)

# The prediction keeps the first methionine, so its numbers are one higher
fixed = [crystal[n]["CA"] for n in range(1, 147)]
moving = [predicted[n + 1]["CA"] for n in range(1, 147)]
fit = Superimposer()
fit.set_atoms(fixed, moving)
print(f"after superposition, {len(fixed)} alpha carbons differ")
print(f"  by {fit.rms:.2f} A on average (root mean square)")
