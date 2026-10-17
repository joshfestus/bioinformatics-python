# Ingram's sickle cell peptide, found in beta globin by local alignment.
from Bio import Align, SeqIO
from Bio.Align import substitution_matrices

beta = next(SeqIO.parse("globins.fasta", "fasta")).seq
sickle_peptide = "VHLTPVEK"

aligner = Align.PairwiseAligner()
aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
aligner.open_gap_score = -10
aligner.extend_gap_score = -0.5

for mode in ("global", "local"):
    aligner.mode = mode
    best = aligner.align(beta, sickle_peptide)[0]
    print(f"{mode:6}  score {best.score:6.1f}")

print(best)
