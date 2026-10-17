"""Solutions to the exercises of Chapter 12.

Run from this folder: uv run ch12.py
Needs MAFFT on the PATH for the exercises that realign.
"""

import math
import shutil
import subprocess
import tempfile
from collections import Counter
from pathlib import Path

from Bio import AlignIO, SeqIO

CH09, CH12 = Path("../ch09"), Path("../ch12")


def entropy(column):
    """Shannon entropy in bits, as conservation.py defines it."""
    counts = Counter(column)
    total = len(column)
    return -sum(n / total * math.log2(n / total)
                for n in counts.values())


def invariant(alignment):
    """Columns where every sequence has the same residue, no gaps."""
    found = []
    for index in range(alignment.get_alignment_length()):
        column = alignment[:, index]
        if "-" not in column and len(set(column)) == 1:
            found.append((index + 1, column[0]))
    return found


def numbered(alignment, reference=0):
    """{alignment column: position in the reference sequence}."""
    places, position = {}, 0
    for index, letter in enumerate(alignment[reference].seq):
        if letter != "-":
            position += 1
        places[index + 1] = position - 1  # the chapter's count
    return places


def mafft(records, into):
    """Align with MAFFT, as Chapter 12 does."""
    with tempfile.TemporaryDirectory() as work:
        source = Path(work) / "in.fasta"
        SeqIO.write(records, source, "fasta")
        done = subprocess.run(["mafft", "--quiet", "--auto", str(source)],
                              capture_output=True, text=True, check=True)
    Path(into).write_text(done.stdout)
    return AlignIO.read(into, "fasta")


def exercise_12_1():
    column = "LLLLIIVV"
    counts = Counter(column)
    by_hand = -sum((n / 8) * math.log2(n / 8) for n in counts.values())
    print("12.1  A column of four leucines, two isoleucines and two")
    print("      valines, by hand:")
    for letter, n in sorted(counts.items()):
        share = n / 8
        print(f"        {letter}: {n}/8 = {share:.2f}, "
              f"-{share:.2f} x log2({share:.2f}) = "
              f"{-share * math.log2(share):.3f} bits")
    print(f"      Adding them gives {by_hand:.3f} bits.")
    print(f"      The entropy function returns {entropy(column):.3f}.")
    print("      One and a half bits out of a possible 4.32 for twenty")
    print("      equally likely residues. The column is not invariant,")
    print("      but all three residues are hydrophobic and similar in")
    print("      size, which entropy cannot see: it counts letters, not")
    print("      chemistry.")


def exercise_12_2():
    everything = AlignIO.read(CH12 / "globins_aligned.fasta", "fasta")
    all_invariant = invariant(everything)
    betas = [record for record in
             SeqIO.parse(CH12 / "globins.fasta", "fasta")
             if record.id.split("|")[2].startswith("HBB")]
    if not shutil.which("mafft"):
        print("\n12.2  MAFFT is not on the PATH; see Appendix A.")
        return
    beta_alignment = mafft(betas, "betas_aligned.fasta")
    beta_invariant = invariant(beta_alignment)
    kept = {letter for _, letter in all_invariant}
    print(f"\n12.2  All 16 globins: {len(all_invariant)} invariant "
          "columns.")
    print(f"      The {len(betas)} beta globins alone: "
          f"{len(beta_invariant)}.")
    print("      Every column invariant across all sixteen is still")
    print("      invariant among the twelve, because a smaller and")
    print("      more closely related set can only agree more often.")
    print("      The sixteen protein columns are the interesting ones:")
    print("      they are what a globin must have to be a globin, and")
    print("      the chapter named them, the two histidines and the")
    print("      phenylalanine that hold the heme among them.")
    print(f"      The extra columns among the betas, "
          f"{len(beta_invariant)} against")
    print(f"      {len(all_invariant)}, are not functional")
    print("      constraints. They are shared ancestry: these twelve")
    print("      proteins are the same gene in different mammals and")
    print("      have not had time to differ.")
    print(f"      Residues conserved in all sixteen: "
          f"{''.join(sorted(kept))}")
    Path("betas_aligned.fasta").unlink()


def exercise_12_3():
    if not shutil.which("mafft"):
        print("\n12.3  MAFFT is not on the PATH; see Appendix A.")
        return
    globins = list(SeqIO.parse(CH12 / "globins.fasta", "fasta"))
    rhodopsin = SeqIO.read(CH09 / "rhodopsin.fasta", "fasta")
    alignment = mafft(globins + [rhodopsin], "intruder.fasta")
    length = alignment.get_alignment_length()
    before = AlignIO.read(CH12 / "globins_aligned.fasta", "fasta")
    row = next(i for i, record in enumerate(alignment)
               if "OPSD" in record.id or "P08100" in record.id)
    gaps = str(alignment[row].seq).count("-")
    globin_gaps = [str(record.seq).count("-")
                   for i, record in enumerate(alignment) if i != row]
    print(f"\n12.3  Adding rhodopsin stretches the alignment from "
          f"{before.get_alignment_length()} to {length} columns.")
    print(f"      Rhodopsin is {len(rhodopsin.seq)} residues against a")
    print("      globin's 141 to 154, so it is the globins that fill")
    print(f"      with gaps: {min(globin_gaps)} to {max(globin_gaps)} "
          "gap characters each,")
    print(f"      while rhodopsin's own row has {gaps}.")
    print(f"      {len(invariant(alignment))} columns stay invariant, "
          f"against {len(invariant(before))} before.")
    print("      That is how you spot an intruder without being told.")
    print("      One sequence forces the alignment to double in")
    print("      length, every other row becomes mostly gaps, and the")
    print("      conserved columns collapse from fourteen to three.")
    print("      An alignment of unrelated sequences is still an")
    print("      alignment and the program will not complain; the")
    print("      evidence that it means nothing is in its shape.")
    Path("intruder.fasta").unlink()


def exercise_12_4():
    alignment = AlignIO.read(CH12 / "globins_aligned.fasta", "fasta")
    places = numbered(alignment)
    wanted = [column for column, position in places.items()
              if position == 63]
    around = wanted[0] if wanted else 66
    first = last = around
    while first > 1 and "-" not in alignment[:, first - 2]:
        first -= 1
    while (last < alignment.get_alignment_length()
           and "-" not in alignment[:, last]):
        last += 1
    print(f"\n12.4  The distal histidine is at column {around}, which "
          f"is position {places[around]}")
    print("      in human beta globin's traditional count.")
    print(f"      The longest gap free stretch containing it runs from")
    print(f"      column {first} to {last}, {last - first + 1} columns.")
    consensus = "".join(Counter(alignment[:, i - 1]).most_common(1)[0][0]
                        for i in range(first, last + 1))
    print(f"      Its consensus is {consensus}")
    print("      motif.py draws that stretch as a logo, where the")
    print("      height of each letter is how much the column has to")
    print("      say. The histidine stands full height because every")
    print("      globin needs it to hold the heme iron, and its")
    print("      neighbors are shorter because they vary.")


def exercise_12_5():
    alignment = AlignIO.read(CH12 / "globins_aligned.fasta", "fasta")
    AlignIO.write(alignment, "globins.aln", "clustal")
    lines = Path("globins.aln").read_text().splitlines()
    print("\n12.5  Biopython writes the file, and the first surprise is")
    print("      what is missing:")
    print(f"        {lines[0]}")
    print(f"        {len(lines)} lines, "
          f"{sum(1 for line in lines if line.strip())} of them not blank")
    print(f"        {len(alignment)} sequences, so each block is "
          f"{len(alignment)} rows")
    marks = [line for line in lines if line
             and not line[0].isalnum()
             and set(line) <= set(" *:.")]
    print(f"        conservation lines found: {len(marks)}")
    print("      Biopython writes valid Clustal format and leaves the")
    print("      conservation line out, because that line is Clustal's")
    print("      own annotation rather than part of the alignment.")
    print("      Nothing is lost: the alignment is the columns, and")
    print("      the marks are a summary anyone can recompute.")
    stars = "".join("*" if "-" not in alignment[:, i]
                    and len(set(alignment[:, i])) == 1 else " "
                    for i in range(alignment.get_alignment_length()))
    print(f"\n      Computing the identical columns ourselves gives "
          f"{stars.count('*')},")
    print(f"      the same {len(invariant(alignment))} columns "
          "conservation.py found.")
    print("      Clustal's other two marks, ':' for strongly similar")
    print("      and '.' for weakly similar, come from its own tables")
    print("      of residue groups rather than from a score, which is")
    print("      why two programs can disagree about them while")
    print("      agreeing exactly about the stars.")
    print("      show_alignment.py prints a number per column instead,")
    print("      which reads less well and can be sorted, thresholded")
    print("      and plotted. That is usually the better trade.")
    Path("globins.aln").unlink()


if __name__ == "__main__":
    exercise_12_1()
    exercise_12_2()
    exercise_12_3()
    exercise_12_4()
    exercise_12_5()
