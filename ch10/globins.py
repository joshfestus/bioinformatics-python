# Human beta globin against mouse beta globin and human alpha globin.
# Sequences: UniProt P68871, P02088 and P69905 (CC BY 4.0).
from Bio import Align, SeqIO
from Bio.Align import substitution_matrices

BLOSUM62 = substitution_matrices.load("BLOSUM62")


def mark(x: str, y: str) -> str:
    """| identical, : similar (scores above zero), . different."""
    if "-" in (x, y):
        return " "
    if x == y:
        return "|"
    return ":" if BLOSUM62[x][y] > 0 else "."


def show(alignment: Align.Alignment, width: int = 60) -> None:
    top, bottom = alignment[0], alignment[1]
    for start in range(0, len(top), width):
        a, b = top[start : start + width], bottom[start : start + width]
        print(a, "".join(mark(x, y) for x, y in zip(a, b)), b, sep="\n")
        print()


proteins = {}
for record in SeqIO.parse("globins.fasta", "fasta"):
    proteins[record.id.split("|")[2]] = record.seq

aligner = Align.PairwiseAligner(mode="global")
aligner.substitution_matrix = BLOSUM62
aligner.open_gap_score = -10
aligner.extend_gap_score = -0.5

human = proteins["HBB_HUMAN"]
for name in ("HBB1_MOUSE", "HBA_HUMAN"):
    best = aligner.align(human, proteins[name])[0]
    counts = best.counts(BLOSUM62)
    print(name)
    print(f"  score      {best.score:6.1f}")
    print(f"  identical  {counts.identities:4} of {best.length}")
    print(f"  similar    {counts.positives:4} of {best.length}")
    print(f"  gaps       {counts.gaps:4}")
    print()

show(best)
