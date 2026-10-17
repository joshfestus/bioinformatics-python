# A neighbor joining tree of sixteen globins from their alignment.
from Bio import AlignIO, Phylo
from Bio.Phylo.TreeConstruction import (
    DistanceCalculator,
    DistanceTreeConstructor,
)

alignment = AlignIO.read("globins_aligned.fasta", "fasta")
for record in alignment:
    record.id = record.id.split("|")[
        2
    ]  # sp|P68871|HBB_HUMAN -> HBB_HUMAN

calculator = DistanceCalculator("blosum62")
distances = calculator.get_distance(alignment)
print("distances from human beta globin:")
for name in (
    "HBB_PANTR",
    "HBB1_MOUSE",
    "HBB_CHICK",
    "HBA_HUMAN",
    "MYG_HUMAN",
):
    print(f"  {name:11} {distances['HBB_HUMAN', name]:.3f}")

constructor = DistanceTreeConstructor(calculator, "nj")
tree = constructor.build_tree(alignment)
tree.root_with_outgroup("MYG_HUMAN")
tree.ladderize()
Phylo.draw_ascii(tree, column_width=70)
Phylo.write(tree, "globins.nwk", "newick")
