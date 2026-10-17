# Bootstrap support: how often each group comes out together when the
# alignment's columns are resampled 100 times.
import random

from Bio import AlignIO
from Bio.Phylo.Consensus import bootstrap_trees
from Bio.Phylo.TreeConstruction import (
    DistanceCalculator,
    DistanceTreeConstructor,
)

GROUPS = {
    "human and chimpanzee": ["HBB_HUMAN", "HBB_PANTR"],
    "mouse and rat": ["HBB1_MOUSE", "HBB1_RAT"],
    "the alpha globins": ["HBA_HUMAN", "HBA_MOUSE", "HBA_CHICK"],
    "dog with pig, cow, horse": [
        "HBB_CANLF",
        "HBB_PIG",
        "HBB_BOVIN",
        "HBB_HORSE",
    ],
    "dog with the primates": ["HBB_CANLF", "HBB_HUMAN", "HBB_PANTR"],
}

random.seed(1)  # the resampling is random; this fixes it for the book
alignment = AlignIO.read("globins_aligned.fasta", "fasta")
for record in alignment:
    record.id = record.id.split("|")[2]

constructor = DistanceTreeConstructor(
    DistanceCalculator("blosum62"), "nj"
)
trees = list(bootstrap_trees(alignment, 100, constructor))
for tree in trees:
    tree.root_with_outgroup("MYG_HUMAN")


def together(tree, names: list[str]) -> bool:
    """True if these leaves, and no others, share a branch of the tree."""
    leaves = [tree.find_any(name) for name in names]
    return tree.is_monophyletic(leaves) is not False


for label, names in GROUPS.items():
    count = sum(together(tree, names) for tree in trees)
    print(f"{label:26} together in {count:3} of {len(trees)} trees")
