"""Solutions to the exercises of Chapter 1.

Run from this folder: uv run ch01.py
Each answer prints what the exercise asks for, and the comments say why.
"""

from pathlib import Path

HBB = Path("../ch01/hbb.txt").read_text().strip()
CODONS = [HBB[i:i + 3] for i in range(0, len(HBB), 3)]


def exercise_1_1():
    """How many GAG codons, and why that is not gene.count("GAG")."""
    print("1.1  GAG as a codon:", CODONS.count("GAG"))
    print("     GAG anywhere in the sequence:", HBB.count("GAG"))
    print("     The second number is larger because it also counts a GAG")
    print("     that straddles two codons, such as the AG of one codon")
    print("     followed by the G of the next. Only a GAG that starts at")
    print("     a position divisible by three is an actual codon.")


def exercise_1_2():
    """The last codon, and the three stop codons."""
    print(f"\n1.2  Last codon: {CODONS[-1]}")
    print("     The three stop codons are TAA, TAG and TGA. This one")
    print(f"     ends in {CODONS[-1]}, so the ribosome stops there.")


def exercise_1_3():
    """GC content of the third position of every codon."""
    third = HBB[2::3]
    gc = (third.count("G") + third.count("C")) / len(third)
    whole = (HBB.count("G") + HBB.count("C")) / len(HBB)
    print(f"\n1.3  Third positions: {len(third)} bases, GC {gc:.1%}")
    print(f"     The whole gene is {whole:.1%}, so the third position is")
    print("     the GC rich one, as Figure 1.2 shows. A change there is")
    print("     often synonymous, so it is the least constrained.")


def exercise_1_4():
    """The alpha chain, NM_000558.5, beside the beta chain."""
    hba = read_fasta("../ch04/globin_cds.fasta", "NM_000558.5")
    gc = (hba.count("G") + hba.count("C")) / len(hba)
    beta_gc = (HBB.count("G") + HBB.count("C")) / len(HBB)
    print(f"\n1.4  Alpha globin: {len(hba)} bases, {len(hba) // 3 - 1} "
          "amino acids")
    print(f"     Beta globin:  {len(HBB)} bases, {len(HBB) // 3 - 1} "
          "amino acids")
    print(f"     GC: alpha {gc:.1%} against beta {beta_gc:.1%}, so alpha")
    print("     is the higher. Alpha is the shorter chain, and the one")
    print("     hemoglobin carries two of alongside two betas.")


def read_fasta(path, wanted):
    """The sequence of the first record whose header holds `wanted`."""
    sequence, keeping = [], False
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            if keeping:
                break
            keeping = wanted in line
        elif keeping:
            sequence.append(line.strip())
    return "".join(sequence)


def exercise_1_5():
    """Every position of GAG, and which of them begin a codon."""
    places = []
    at = HBB.find("GAG")
    while at != -1:
        places.append(at)
        at = HBB.find("GAG", at + 1)
    starts = [p for p in places if p % 3 == 0]
    print(f"\n1.5  GAG occurs at {len(places)} positions: "
          f"{', '.join(str(p) for p in places[:8])} ...")
    print(f"     First occurrence: {places[0]}, which is "
          f"{'' if places[0] % 3 == 0 else 'not '}the start of a codon.")
    print(f"     Codon 7 begins at position 18, so the first GAG is "
          f"{'' if places[0] == 18 else 'not '}the sickle codon.")
    print(f"     {len(starts)} of the {len(places)} begin a codon, the")
    print("     same count as Exercise 1.1 found.")


if __name__ == "__main__":
    exercise_1_1()
    exercise_1_2()
    exercise_1_3()
    exercise_1_4()
    exercise_1_5()
