"""Small tools for DNA sequences, written from scratch."""

PAIRS = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}

# The standard genetic code, as NCBI publishes it: the amino acids for
# all 64 codons, with the letters of each codon taken in the order TCAG.
BASES = "TCAG"
AMINO_ACIDS = (
    "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
)

GENETIC_CODE = {}
for first in BASES:
    for second in BASES:
        for third in BASES:
            GENETIC_CODE[first + second + third] = AMINO_ACIDS[
                len(GENETIC_CODE)
            ]


def reverse_complement(dna: str) -> str:
    """The other strand, read 5' to 3'."""
    return "".join(PAIRS[base] for base in reversed(dna))


def gc_content(dna: str) -> float:
    """The share of the sequence that is G or C, from 0 to 1."""
    return (dna.count("G") + dna.count("C")) / len(dna)


def translate(dna: str) -> str:
    """Translate from the first letter until the first stop codon."""
    protein = []
    for start in range(0, len(dna) - 2, 3):
        amino_acid = GENETIC_CODE[dna[start : start + 3]]
        if amino_acid == "*":
            break
        protein.append(amino_acid)
    return "".join(protein)
