"""A tree of the sampled genomes, and a check that the lineages hold up.

The distances are computed with NumPy, because comparing 115 genomes to
one another one letter at a time in Python would take minutes: every pair
of genomes is 29,903 comparisons, and there are 6,555 pairs.

uv run tree.py             reads aligned.fasta, writes genomes.nwk
"""

import csv

import numpy as np
from Bio import Phylo, SeqIO
from Bio.Phylo.TreeConstruction import (
    DistanceMatrix,
    DistanceTreeConstructor,
)

REFERENCE = "NC_045512.2"


def load(path):
    names, rows = [], []
    for record in SeqIO.parse(path, "fasta"):
        names.append(record.id)
        rows.append(np.frombuffer(str(record.seq).upper().encode(),
                                  dtype="S1"))
    return names, np.array(rows)


def distances(letters):
    """The share of compared positions at which two genomes differ."""
    known = np.isin(letters, [b"A", b"C", b"G", b"T"])
    count = len(letters)
    matrix = np.zeros((count, count))
    for i in range(count):
        both = known[i] & known
        differs = (letters[i] != letters) & both
        shared = both.sum(axis=1)
        matrix[i] = np.where(shared > 0, differs.sum(axis=1) / shared, 0)
    return (matrix + matrix.T) / 2  # make it exactly symmetric


def main():
    names, letters = load("aligned.fasta")
    matrix = distances(letters)
    print(f"{len(names)} genomes, distances from the reference:")
    reference = names.index(REFERENCE)
    order = np.argsort(matrix[reference])
    for index in list(order[:2]) + list(order[-2:]):
        print(f"  {names[index]:12} {matrix[reference, index]:.5f}")

    lower = [[matrix[i, j] for j in range(i + 1)]
             for i in range(len(names))]
    tree = DistanceTreeConstructor().nj(DistanceMatrix(names, lower))
    tree.root_with_outgroup(REFERENCE)
    tree.ladderize()
    Phylo.write(tree, "genomes.nwk", "newick")

    with open("samples.tsv") as f:
        lineage = {row["accession"]: row["lineage"]
                   for row in csv.DictReader(f, delimiter="\t")}
    print("\ndoes each lineage form one branch of the tree?")
    for name in dict.fromkeys(lineage.values()):
        members = [tip for tip in tree.get_terminals()
                   if lineage.get(tip.name) == name]
        if len(members) < 2:
            continue
        ancestor = tree.common_ancestor(members)
        under = ancestor.get_terminals()
        extra = [tip.name for tip in under
                 if lineage.get(tip.name) != name]
        verdict = ("one branch" if not extra
                   else f"{len(extra)} others among it")
        print(f"  {name:10} {len(members):3} genomes  {verdict}")


if __name__ == "__main__":
    main()
