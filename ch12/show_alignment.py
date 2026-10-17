# Print part of the globin alignment with a conservation line under it.
import sys

from Bio import AlignIO

alignment = AlignIO.read("globins_aligned.fasta", "fasta")
start, end = (int(n) for n in sys.argv[1:3])  # columns, counted from 1

for record in alignment:
    name = record.id.split("|")[2]
    print(f"{name:11} {record.seq[start - 1 : end]}")

marks = []
for column in range(start - 1, end):
    letters = set(alignment[:, column])
    beta_letters = set(alignment[:12, column])
    if len(letters) == 1 and "-" not in letters:
        marks.append("*")  # identical in all sixteen
    elif len(beta_letters) == 1 and "-" not in beta_letters:
        marks.append(".")  # identical in the twelve beta globins
    else:
        marks.append(" ")
print(" " * 12 + "".join(marks))
