# In the sickle hemoglobin crystal, which residues of another hemoglobin
# molecule touch valine 6?
from Bio.PDB import MMCIFParser, NeighborSearch

model = MMCIFParser(QUIET=True).get_structure("2HBS", "../data/2HBS.cif")[
    0
]
molecule = {chain: 1 if chain in "ABCD" else 2 for chain in "ABCDEFGH"}
search = NeighborSearch(list(model.get_atoms()))

for chain_id in ("B", "D", "F", "H"):
    valine = model[chain_id][6]
    closest = {}
    for atom in valine:
        if atom.get_id() in ("N", "CA", "C", "O"):  # side chain only
            continue
        for near in search.search(atom.coord, 4.5):
            residue = near.get_parent()
            chain = residue.get_parent().id
            if (
                molecule[chain] == molecule[chain_id]
                or residue.id[0] != " "
            ):
                continue
            key = (chain, residue.get_resname(), residue.id[1])
            closest[key] = min(closest.get(key, 99), atom - near)
    print(f"valine 6 of chain {chain_id}: {len(closest)} contacts")
    for (chain, name, number), distance in sorted(
        closest.items(), key=lambda item: item[0][2]
    ):
        print(f"  chain {chain} {name} {number:3}  {distance:.1f} A")
