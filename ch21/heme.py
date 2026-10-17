# How the heme's iron sits between the two histidines in each beta chain
# of normal deoxyhemoglobin.
from Bio.PDB import MMCIFParser

model = MMCIFParser(QUIET=True).get_structure("2HHB", "../data/2HHB.cif")[
    0
]
for chain_id in ("B", "D"):
    chain = model[chain_id]
    heme = next(r for r in chain if r.get_resname() == "HEM")
    iron = heme["FE"]
    proximal = (
        iron - chain[92]["NE2"]
    )  # subtracting atoms gives angstroms
    distal = iron - chain[63]["NE2"]
    print(
        f"chain {chain_id}: iron to His92 {proximal:.2f} A, "
        f"to His63 {distal:.2f} A"
    )
