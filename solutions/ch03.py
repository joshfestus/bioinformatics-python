"""Solutions to the exercises of Chapter 3.

Run from this folder: uv run ch03.py
These are pencil and paper exercises, so each answer gives what you
should have written and then checks it with Biopython, which is the
habit the rest of the book is built on.
"""

from pathlib import Path

from Bio.Seq import Seq

HBB = Path("../ch01/hbb.txt").read_text().strip()


def exercise_3_1():
    """Translating the first eight codons by hand."""
    codons = "ATG GTG CAT CTG ACT CCT GAG GAG"
    by_hand = "M V H L T P E E"
    checked = Seq(codons.replace(" ", "")).translate()
    print("3.1  By hand, reading T as U:")
    print(f"       {codons}")
    print(f"       {by_hand}")
    print(f"     Biopython agrees: {' '.join(str(checked))}")
    print("     Methionine, valine, histidine, leucine, threonine,")
    print("     proline, then two glutamic acids. The second of those")
    print("     glutamic acids is the one sickle cell disease changes.")


def exercise_3_2():
    """Reverse complements of the normal and sickle sequences."""
    normal, sickle = "CCTGAGGAG", "CCTGTGGAG"
    rc_normal = str(Seq(normal).reverse_complement())
    rc_sickle = str(Seq(sickle).reverse_complement())
    differs = [i for i, (a, b) in enumerate(zip(rc_normal, rc_sickle))
               if a != b]
    print(f"\n3.2  {normal} becomes {rc_normal}")
    print(f"     {sickle} becomes {rc_sickle}")
    print(f"     They differ at position {differs[0] + 1} of nine, where")
    print(f"     the normal strand has {rc_normal[differs[0]]} and the")
    print(f"     sickle has {rc_sickle[differs[0]]}.")
    print("     On the coding strand the mutation reads A to T. On the")
    print("     other strand the same mutation reads T to A. A read from")
    print("     either strand carries the same change, which is why a")
    print("     variant caller reports the position once, against the")
    print("     reference strand, however the read happened to be made.")


def exercise_3_3():
    """Synonymous, missense or nonsense."""
    print("\n3.3  GAG codes for glutamic acid, E.")
    for codon in ("GAA", "GAC", "TAG", "GCG"):
        amino = Seq(codon).translate()
        if str(amino) == "*":
            verdict = "nonsense, it is a stop codon"
        elif str(amino) == "E":
            verdict = "synonymous, still glutamic acid"
        else:
            verdict = f"missense, now {amino}"
        print(f"       GAG to {codon}: {verdict}")
    print("     GAA differs in the third position, which is where the")
    print("     genetic code is most forgiving.")


def exercise_3_4():
    """The same letters read in three frames."""
    piece = "ATGGTGCATCTGACTCCT"
    print(f"\n3.4  {piece}, read three ways:")
    for offset in (0, 1, 2):
        frame = piece[offset:]
        frame = frame[:len(frame) // 3 * 3]
        codons = [frame[i:i + 3] for i in range(0, len(frame), 3)]
        protein = str(Seq(frame).translate())
        stop = protein.find("*")
        ends = f"stops at amino acid {stop + 1}" if stop >= 0 \
            else "runs to the end of the piece with no stop"
        print(f"       frame {offset + 1}: {' '.join(codons)}")
        print(f"                {'   '.join(protein)}   {ends}")
    print("     Only the first frame is the real one. The other two are")
    print("     what the same letters would say if the ribosome began")
    print("     one or two bases late, which is why a deletion that is")
    print("     not a multiple of three is so destructive.")


def exercise_3_5():
    """Deleting CTTT at 126 to 129, the frameshift of codons 41 and 42."""
    removed = HBB[:125] + HBB[129:]
    protein = str(Seq(removed[:len(removed) // 3 * 3]).translate())
    stop = protein.find("*")
    print(f"\n3.5  Removing CTTT leaves {len(removed)} bases, "
          f"{len(HBB) - len(removed)} fewer.")
    print(f"     The first stop codon is now amino acid {stop + 1}, so")
    print(f"     protein would be {stop} amino acids instead of "
          f"{len(HBB) // 3 - 1}.")
    print(f"     It still begins {protein[:7]}, because the deletion is")
    print("     at codon 41 and everything before it is unchanged. After")
    print("     that point every codon is read in the wrong frame.")
    print("     This is a real beta thalassemia mutation, and a protein")
    print("     this short does not work at all.")


if __name__ == "__main__":
    exercise_3_1()
    exercise_3_2()
    exercise_3_3()
    exercise_3_4()
    exercise_3_5()
