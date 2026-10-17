# What is inside a structure file: normal and sickle deoxyhemoglobin.
from Bio.PDB import MMCIFParser

parser = MMCIFParser(QUIET=True)
for code in ("2HHB", "2HBS"):
    model = parser.get_structure(code, f"../data/{code}.cif")[0]
    print(f"{code}: {len(model)} chains")
    for chain in model:
        residues = [r for r in chain if r.id[0] == " "]  # amino acids
        others = sorted(
            {r.get_resname() for r in chain if r.id[0] != " "}
        )
        kind = "beta" if len(residues) == 146 else "alpha"
        line = f"  {chain.id}  {kind:5} {len(residues)} residues"
        if kind == "beta":
            line += f", residue 6 is {chain[6].get_resname()}"
        print(f"{line}; also {', '.join(others)}")
