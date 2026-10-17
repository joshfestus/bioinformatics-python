# How strongly each column of the globin alignment is conserved, and which
# columns are identical in all sixteen proteins.
import math
from collections import Counter

from Bio import AlignIO

alignment = AlignIO.read("globins_aligned.fasta", "fasta")
human = alignment[0].seq  # human beta globin, the reference for numbering


def entropy(column: str) -> float:
    """Shannon entropy in bits: 0 if all letters agree, higher if not."""
    counts = Counter(column)
    total = len(column)
    return -sum(n / total * math.log2(n / total) for n in counts.values())


position = 0  # traditional count: valine after the first methionine is 1
rows = []
for index in range(alignment.get_alignment_length()):
    column = alignment[:, index]
    if human[index] != "-":
        position += 1
    rows.append((index + 1, position - 1, column, entropy(column)))

invariant = [r for r in rows if len(set(r[2])) == 1 and r[2][0] != "-"]
mean = sum(r[3] for r in rows) / len(rows)
print(f"{len(rows)} columns, mean entropy {mean:.2f} bits")
print(f"{len(invariant)} columns identical in all sixteen proteins:")
for column, beta, letters, _ in invariant:
    print(f"  column {column:3}  beta {letters[0]}{beta}")
