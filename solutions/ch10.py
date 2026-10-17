"""Solutions to the exercises of Chapter 10.

Run from this folder: uv run ch10.py
Exercise 10.5 downloads one protein from UniProt; the rest are local.
"""

import sys
import urllib.error
import urllib.request
from pathlib import Path

from Bio import Align, SeqIO
from Bio.Align import substitution_matrices

sys.path.insert(0, str(Path("../ch10").resolve()))

import nw  # noqa: E402

CH10 = Path("../ch10")
BLOSUM62 = substitution_matrices.load("BLOSUM62")
UNIPROT = "https://rest.uniprot.org/uniprotkb"


def globins():
    return {record.id.split("|")[2]: str(record.seq)
            for record in SeqIO.parse(CH10 / "globins.fasta", "fasta")}


def align_pair(a, b, gap=nw.GAP):
    """Chapter 10's own aligner, with a gap penalty you can change."""
    kept = nw.GAP
    nw.GAP = gap
    try:
        grid = nw.fill(a, b)
        top, bottom = nw.traceback(a, b, grid)
        return grid[len(a)][len(b)], top, bottom
    finally:
        nw.GAP = kept


def exercise_10_1():
    print("10.1  ATGAGCA against AGACCA, filled by nw.py:")
    for a, b in (("ATGAGCA", "AGACCA"), ("GATTACA", "GCATGCA")):
        score, top, bottom = align_pair(a, b)
        matches = sum(x == y for x, y in zip(top, bottom))
        print(f"       {a} against {b}: score {score}")
        print(f"         {top}")
        marks = "".join("|" if x == y else " "
                        for x, y in zip(top, bottom))
        print(f"         {marks}")
        print(f"         {bottom}")
        print(f"         {matches} identical columns of {len(top)}")
    print("      Filling the grid by hand is worth doing once. Every")
    print("      cell is the best of three: the diagonal plus a match")
    print("      or mismatch, or either neighbor plus a gap. The last")
    print("      cell holds the score of the best alignment, and the")
    print("      path back through the choices is the alignment.")


def exercise_10_2():
    a, b = "GCATGCT", "GATTACA"
    print("\n10.2  The same pair under two gap penalties:")
    for gap in (-2, -1):
        score, top, bottom = align_pair(a, b, gap)
        gaps = top.count("-") + bottom.count("-")
        matches = sum(x == y for x, y in zip(top, bottom))
        print(f"       gap {gap}: score {score}, {gaps} gap columns, "
              f"{matches} identical")
        print(f"         {top}")
        print(f"         {bottom}")
    print("      Cheap gaps buy matches with insertions and deletions,")
    print("      so the alignment grows more gapped and scores higher.")
    print("      For two sequences this short and this different, the")
    print("      gapped version is not more believable: it is the")
    print("      aligner being allowed to invent history. Real gap")
    print("      penalties come from how often insertions actually")
    print("      happen, which is rarely, and that is why Chapter 10's")
    print("      default makes a gap cost twice a mismatch.")


def smith_waterman(a, b, match=1, mismatch=-1, gap=-2):
    """Needleman-Wunsch with a floor of zero: local alignment."""
    grid = [[0] * (len(b) + 1) for _ in range(len(a) + 1)]
    best, best_at = 0, (0, 0)
    for i in range(1, len(a) + 1):
        for j in range(1, len(b) + 1):
            step = match if a[i - 1] == b[j - 1] else mismatch
            grid[i][j] = max(0,
                             grid[i - 1][j - 1] + step,
                             grid[i - 1][j] + gap,
                             grid[i][j - 1] + gap)
            if grid[i][j] > best:
                best, best_at = grid[i][j], (i, j)
    top, bottom = [], []
    i, j = best_at
    while i > 0 and j > 0 and grid[i][j] > 0:
        step = match if a[i - 1] == b[j - 1] else mismatch
        if grid[i][j] == grid[i - 1][j - 1] + step:
            top.append(a[i - 1])
            bottom.append(b[j - 1])
            i, j = i - 1, j - 1
        elif grid[i][j] == grid[i - 1][j] + gap:
            top.append(a[i - 1])
            bottom.append("-")
            i -= 1
        else:
            top.append("-")
            bottom.append(b[j - 1])
            j -= 1
    return best, "".join(reversed(top)), "".join(reversed(bottom)), (i, j)


def exercise_10_3():
    peptide = "VHLTPVEK"
    beta = globins()["HBB_HUMAN"][:30]
    score, top, bottom, start = smith_waterman(peptide, beta)
    print(f"\n10.3  Smith-Waterman, {peptide} against the first 30")
    print(f"      residues of human beta globin: score {score}")
    print(f"        {top}")
    print(f"        {bottom}")
    print(f"      It begins at residue {start[1] + 1} of the globin.")
    print("      The change from global to local is three small ones:")
    print("      a floor of zero in every cell, a traceback that starts")
    print("      at the best cell anywhere rather than the corner, and")
    print("      a stop at the first zero. That is the whole difference")
    print("      between Needleman-Wunsch and Smith-Waterman, and it is")
    print("      why a local aligner finds a short peptide inside a")
    print("      long protein while a global one would not.")


def exercise_10_4():
    beta = globins()["HBB_HUMAN"]
    sickle = beta[:6] + "V" + beta[7:]
    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = BLOSUM62
    aligner.open_gap_score, aligner.extend_gap_score = -10, -0.5
    same = aligner.score(beta, beta)
    changed = aligner.score(beta, sickle)
    print(f"\n10.4  Beta globin against itself scores {same:.0f}.")
    print(f"      Against the sickle protein it scores {changed:.0f}, "
          f"{same - changed:.0f} less.")
    same_score = BLOSUM62["E"]["E"]
    swap_score = BLOSUM62["E"]["V"]
    print(f"      BLOSUM62 scores E against E as {same_score:.0f} and E")
    print(f"      against V as {swap_score:.0f}, a difference of "
          f"{same_score - swap_score:.0f}, which is the whole gap.")
    print("      One residue in 147, and the alignment barely notices.")
    print("      Alignment scores measure similarity, not consequence:")
    print("      this difference of five points is a disease, and")
    print("      Chapter 21 shows why by looking at where the residue")
    print("      sits rather than how much the score moved.")


def exercise_10_5():
    print("\n10.5  Chicken beta globin, P02112, against the others:")
    path = Path("P02112.fasta")
    sequences = globins()
    if not path.exists():
        try:
            request = urllib.request.Request(
                f"{UNIPROT}/P02112.fasta",
                headers={"User-Agent": "python"})
            with urllib.request.urlopen(request, timeout=120) as reply:
                path.write_text(reply.read().decode())
        except (urllib.error.URLError, TimeoutError) as problem:
            print(f"      Could not download it ({problem}).")
            return
    chicken = str(SeqIO.read(path, "fasta").seq)
    path.unlink()
    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = BLOSUM62
    aligner.open_gap_score, aligner.extend_gap_score = -10, -0.5
    human = sequences["HBB_HUMAN"]
    for name, other in (("mouse", sequences["HBB1_MOUSE"]),
                        ("chicken", chicken)):
        alignment = aligner.align(human, other)[0]
        identical = sum(a == b
                        for a, b in zip(alignment[0], alignment[1]))
        print(f"       human against {name:8}: score "
              f"{alignment.score:6.1f}, {identical} of "
              f"{len(alignment[0])} identical "
              f"({identical / len(alignment[0]):.0%})")
    print("      Mouse is closer. Mammals and birds last shared an")
    print("      ancestor about 310 million years ago, while humans and")
    print("      mice shared one about 90 million years ago, so the")
    print("      sequences have had three times as long to diverge.")
    print("      The number of differences follows the time, which is")
    print("      the assumption Chapter 13 builds trees on.")


def identities(alignment):
    """Percent identity counted three ways."""
    top, bottom = alignment[0], alignment[1]
    identical = sum(a == b for a, b in zip(top, bottom))
    ungapped = sum(1 for a, b in zip(top, bottom) if "-" not in (a, b))
    shorter = min(len(top.replace("-", "")), len(bottom.replace("-", "")))
    return {"over the alignment": identical / len(top),
            "over the shorter sequence": identical / shorter,
            "over ungapped columns": identical / ungapped}


def exercise_10_6():
    sequences = globins()
    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = BLOSUM62
    aligner.open_gap_score, aligner.extend_gap_score = -10, -0.5
    print("\n10.6  Percent identity, counted three ways:")
    pairs = (("HBB_HUMAN", "HBB1_MOUSE"), ("HBB_HUMAN", "HBA_HUMAN"))
    for first, second in pairs:
        alignment = aligner.align(sequences[first], sequences[second])[0]
        print(f"       {first} against {second}:")
        for how, value in identities(alignment).items():
            print(f"         {how:26} {value:.1%}")
    print("      For two proteins of the same length with no gaps the")
    print("      three agree. For beta against alpha, which needs gaps,")
    print("      they differ by several points, and a paper that says")
    print("      only 'percent identity' has not said which it means.")
    print("      State the denominator, as this book's figures do.")


if __name__ == "__main__":
    exercise_10_1()
    exercise_10_2()
    exercise_10_3()
    exercise_10_4()
    exercise_10_5()
    exercise_10_6()
