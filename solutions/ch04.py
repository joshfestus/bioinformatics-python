"""Solutions to the exercises of Chapter 4.

Run from this folder: uv run ch04.py
The chapter's own module does the work where it can, because reusing
what you wrote last week is the point of writing it as a module.
"""

import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path("../ch04").resolve()))

from seqtools import GENETIC_CODE, gc_content, translate  # noqa: E402

HBB = Path("../ch01/hbb.txt").read_text().strip()
IUPAC = {"R": "Y", "Y": "R", "S": "S", "W": "W", "K": "M", "M": "K"}


def frames(dna):
    """The translation of each of the three forward reading frames."""
    return [translate(dna[offset:]) for offset in (0, 1, 2)]


def exercise_4_1():
    piece = HBB[:60]
    print("4.1  The first 60 bases of HBB in three frames:")
    for offset, protein in enumerate(frames(piece)):
        bases = len(piece) - offset
        complete = bases // 3
        stopped = len(protein) < complete
        print(f"       frame {offset + 1}: {protein}")
        print(f"                {'stops early' if stopped else 'no stop'}"
              f", {len(protein)} amino acids of {complete} codons")
    print("     Only frame 1 translates all twenty codons without")
    print("     hitting a stop, which is what tells you it is the frame")
    print("     the ribosome uses.")


def exercise_4_2():
    codons = Counter(HBB[i:i + 3] for i in range(0, len(HBB) - 2, 3))
    leucine = {codon: count for codon, count in codons.items()
               if GENETIC_CODE[codon] == "L"}
    favorite = max(leucine, key=leucine.get)
    print(f"\n4.2  HBB uses {len(codons)} different codons.")
    print(f"     Its leucine codons: {leucine}")
    print(f"     The commonest is {favorite}, used {leucine[favorite]} "
          "times.")
    print("     Six codons mean leucine and this gene does not use them")
    print("     evenly, which is codon bias: the same protein can be")
    print("     written many ways, and species and genes have habits.")


def exercise_4_3():
    """The coding sequence inside the whole messenger RNA."""
    mrna = read_record("../ch04/globin_cds.fasta", "NM_000518.5")
    start = mrna.find("ATG")
    protein = translate(mrna[start:])
    print(f"\n4.3  The file here holds the coding sequence itself, "
          f"{len(mrna)} bases,")
    print(f"     so its first ATG is at position {start + 1} and the")
    print(f"     protein is {len(protein)} amino acids beginning "
          f"{protein[:6]}.")
    print("     Downloading the whole mRNA instead, the header says")
    print("     [location=51..494] and the first ATG is at 51, so the")
    print("     two agree. UniProt P68871 is 147 amino acids with the")
    print("     initial methionine, which this matches.")
    print("     Beware: the first ATG is not always the real start. A")
    print("     messenger RNA can carry an upstream ATG in another")
    print("     frame, which is why the record states the location.")


def read_record(path, wanted):
    sequence, keeping = [], False
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            if keeping:
                break
            keeping = wanted in line
        elif keeping:
            sequence.append(line.strip())
    return "".join(sequence)


def reverse_complement(dna):
    """Chapter 4's function with the ambiguity codes added."""
    pairs = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N", **IUPAC}
    return "".join(pairs[base] for base in reversed(dna))


def exercise_4_4():
    print("\n4.4  Each code stands for a set of bases, so its complement")
    print("     is the set of their complements:")
    meaning = {"R": "A or G", "Y": "C or T", "S": "G or C",
               "W": "A or T", "K": "G or T", "M": "A or C"}
    for code, partner in IUPAC.items():
        print(f"       {code} ({meaning[code]:7}) pairs with {partner} "
              f"({meaning[partner]})")
    print("     S and W are their own complements, because G and C swap")
    print("     to each other and so do A and T.")
    print(f"     reverse_complement('ARYGS') is "
          f"{reverse_complement('ARYGS')}")


def strict_read(path):
    """read_fasta, refusing anything that is not a real base."""
    allowed = set("ACGTN")
    records, name = {}, None
    for number, line in enumerate(Path(path).read_text().splitlines(), 1):
        if line.startswith(">"):
            name = line[1:].split()[0]
            records[name] = ""
            continue
        for letter in line.strip().upper():
            if letter not in allowed:
                raise ValueError(
                    f"{path} line {number}, record {name}: "
                    f"{letter!r} is not a base")
        records[name] += line.strip().upper()
    return records


def exercise_4_5():
    print("\n4.5  The reader now names the record and the letter:")
    damaged = Path("damaged.fasta")
    damaged.write_text(">test_record\nACGTACGTXACGT\n")
    try:
        strict_read(damaged)
    except ValueError as problem:
        print(f"       ValueError: {problem}")
    damaged.unlink()
    print("     Saying which record and which letter is the whole point.")
    print("     An error that says only 'invalid sequence' leaves you")
    print("     searching a file of millions of letters by hand.")


def exercise_4_6():
    """Where the mouse coding sequence differs from UniProt P02088."""
    mouse = read_record("../ch04/globin_cds.fasta", "Hbb-bt")
    protein = translate(mouse)
    # the real P02088 record, as Chapter 4's check_uniprot.py reads it,
    # never a sequence typed from memory
    uniprot = read_record("../ch10/globins.fasta", "P02088")
    differences = [(i + 1, a, b) for i, (a, b)
                   in enumerate(zip(protein, uniprot)) if a != b]
    print(f"\n4.6  The translated mouse sequence is {len(protein)} amino")
    print(f"     acids and differs from UniProt P02088 at "
          f"{len(differences)} of them:")
    for position, ours, theirs in differences:
        codon = mouse[(position - 1) * 3:(position - 1) * 3 + 3]
        options = [c for c, a in GENETIC_CODE.items() if a == theirs]
        close = [c for c in options
                 if sum(x != y for x, y in zip(c, codon)) == 1]
        print(f"       {position:4} {ours} instead of {theirs}, codon "
              f"{codon}")
        nearby = ", ".join(close) or "none"
        print(f"            one letter away: {nearby}")
    print("     These are the documented differences between the S and D")
    print("     alleles of mouse beta globin, not errors. The strain the")
    print("     sequence came from carries a different allele from the")
    print("     one UniProt chose, which is a thing to check before")
    print("     calling any difference a mutation.")


if __name__ == "__main__":
    exercise_4_1()
    exercise_4_2()
    exercise_4_3()
    exercise_4_4()
    exercise_4_5()
    exercise_4_6()
