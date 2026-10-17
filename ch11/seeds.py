# BLAST's first step by hand: seed words shared by human beta and alpha
# globin, scored with BLOSUM62.
from Bio import SeqIO
from Bio.Align import substitution_matrices

WORD = 3  # BLAST's word length for proteins
THRESHOLD = 11  # a pair of words must score at least this to seed

blosum62 = substitution_matrices.load("BLOSUM62")
beta = str(SeqIO.read("hbb.fasta", "fasta").seq)
alpha = str(SeqIO.read("hba_human.fasta", "fasta").seq)


def word_score(a: str, b: str) -> float:
    return sum(blosum62[x][y] for x, y in zip(a, b))


seeds = []
for i in range(len(beta) - WORD + 1):
    for j in range(len(alpha) - WORD + 1):
        score = word_score(beta[i : i + WORD], alpha[j : j + WORD])
        if score >= THRESHOLD:
            seeds.append((i, j, score))

pairs = (len(beta) - WORD + 1) * (len(alpha) - WORD + 1)
print(f"{pairs:,} pairs of words compared, {len(seeds)} seeds")
# The true alignment (Chapter 10) keeps i and j within a few positions
on_diagonal = [s for s in seeds if abs(s[0] - s[1]) <= 8]
print(
    f"{len(on_diagonal)} of them lie within 8 positions of the diagonal"
)
for i, j, score in on_diagonal[:8]:
    print(
        f"  beta {i + 1:3} {beta[i : i + WORD]}   alpha {j + 1:3} "
        f"{alpha[j : j + WORD]}   score {score:4.0f}"
    )
