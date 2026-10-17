"""Solutions to the exercises of Chapter 6.

Run from this folder: uv run ch06.py
Exercise 6.4 needs the network and is skipped unless you pass --net.
"""

import random
import sys
import time
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq

sys.path.insert(0, str(Path("../ch04").resolve()))
sys.path.insert(0, str(Path("../ch06").resolve()))

from orfs import find_orfs  # noqa: E402
from seqtools import translate as ours  # noqa: E402

REGION = SeqIO.read("../ch06/lac_region.gb", "genbank")
OFFSET = 360_000  # the region begins here in NC_000913.3


def coding_sequence(name):
    """The CDS of one gene, cut out of the region and on its strand."""
    for feature in REGION.features:
        if (feature.type == "CDS"
                and feature.qualifiers.get("gene", [""])[0] == name):
            return str(feature.extract(REGION.seq))
    raise KeyError(name)


def exercise_6_1():
    lacz = coding_sequence("lacZ")  # 3,075 bases, one stop at the end
    rounds = 1000
    started = time.perf_counter()
    for _ in range(rounds):
        ours(lacz)
    mine = time.perf_counter() - started
    sequence = Seq(lacz)
    started = time.perf_counter()
    for _ in range(rounds):
        sequence.translate()
    theirs = time.perf_counter() - started
    print(f"6.1  {rounds} translations of {len(lacz)} bases:")
    print(f"       Chapter 4's translate: {mine:6.2f} seconds")
    print(f"       Biopython's:           {theirs:6.2f} seconds")
    faster = mine / theirs if theirs < mine else theirs / mine
    winner = "Biopython" if theirs < mine else "ours"
    print(f"     {winner} is {faster:.1f} times faster.")
    print("     Biopython translates with a lookup table in C rather")
    print("     than a Python loop. Writing it yourself in Chapter 4 was")
    print("     to understand it; using theirs afterwards is to get work")
    print("     done, and that is the right order.")


def starts_anywhere(dna, min_length=100):
    """find_orfs, but accepting GTG and TTG as starts as well as ATG."""
    found = []
    for strand, sequence in ((+1, dna), (-1, dna.reverse_complement())):
        for frame in range(3):
            usable = (len(sequence) - frame) // 3 * 3
            piece = sequence[frame:frame + usable]
            protein = str(piece.translate(11))
            codons = [str(piece[i * 3:i * 3 + 3])
                      for i in range(len(protein))]
            position = 0
            while position < len(protein):
                if codons[position] not in ("ATG", "GTG", "TTG"):
                    position += 1
                    continue
                stop = protein.find("*", position)
                if stop == -1:
                    break
                if stop - position >= min_length:
                    first = frame + 3 * position
                    last = frame + 3 * stop + 3
                    if strand == -1:
                        first, last = len(dna) - last, len(dna) - first
                    found.append((first + 1, last, strand,
                                  stop - position, codons[position]))
                    position = stop + 1
                else:
                    position += 1
    return sorted(found)


def annotated():
    """{gene name: (start, end, strand)} from the GenBank record."""
    genes = {}
    for feature in REGION.features:
        if feature.type == "CDS" and "gene" in feature.qualifiers:
            genes[feature.qualifiers["gene"][0]] = (
                int(feature.location.start) + 1,
                int(feature.location.end),
                feature.location.strand)
    return genes


def stop_codon_at(start, end, strand):
    """Where the stop codon sits, which is the end an ORF shares.

    A gene on the plus strand ends at its highest coordinate, but one on
    the minus strand is read the other way, so its stop codon is at the
    lowest. Matching an ORF to a gene means matching the stop, because
    that is the end both agree on; the start is what an ORF finder has
    to guess.
    """
    return end if strand == 1 else start


def exercise_6_2():
    genes = annotated()
    with_atg = find_orfs(REGION.seq)
    with_any = starts_anywhere(REGION.seq)
    print("\n6.2  Accepting GTG and TTG as starts as well as ATG:")
    for name in ("lacZ", "lacY", "lacA", "lacI"):
        start, end, strand = genes[name]
        stop = stop_codon_at(start, end, strand)
        real_start = end if strand == -1 else start
        real_codon = coding_sequence(name)[:3]
        before = [o for o in with_atg
                  if stop_codon_at(o[0], o[1], o[2]) == stop]
        after = [o for o in with_any
                 if stop_codon_at(o[0], o[1], o[2]) == stop]
        print(f"       {name}: annotated start {real_codon}, "
              f"{len(before)} ORF with ATG only, {len(after)} with all")
        for label, orfs in (("ATG only", before), ("all three", after)):
            if not orfs:
                print(f"            {label:9}: nothing at this length")
                continue
            orf = orfs[0]
            found = orf[1] if strand == -1 else orf[0]
            codon = orf[4] if len(orf) > 4 else "ATG"
            # positive means later in the reading direction, which for a
            # minus strand gene is a lower coordinate
            along = (real_start - found) if strand == -1 \
                else (found - real_start)
            codons_out = along // 3
            if codons_out == 0:
                where = "exact"
            elif codons_out > 0:
                where = f"{codons_out} codons inside the gene"
            else:
                where = f"{-codons_out} codons upstream"
            print(f"            {label:9}: begins {codon}, {where}")
    print("     The two failures of Chapter 6 are different failures.")
    print("     lacA begins TTG, and the ATG only finder ran back to an")
    print("     earlier ATG in the same frame, so its ORF was too long.")
    print("     Allowing TTG does not fix that, because the earlier ATG")
    print("     is still the first start the finder meets.")
    print("     lacI begins GTG, which the ATG only finder could not")
    print("     see, so it took the first ATG inside the gene and")
    print("     reported an ORF too short. Allowing GTG moves it almost")
    print("     to the real start, missing it by a GTG a little further")
    print("     upstream in the same frame.")
    print("     More start codons means more candidates, not a better")
    print("     choice between them. A real annotation uses the")
    print("     ribosome binding site a few bases upstream and the")
    print("     similarity to known proteins, neither of which an ORF")
    print("     finder looks at.")


def not_annotated(min_length):
    """ORFs at this threshold not sharing a stop with a real gene."""
    stops = {stop_codon_at(*gene) for gene in annotated().values()}
    found = find_orfs(REGION.seq, min_length)
    return [orf for orf in found
            if stop_codon_at(orf[0], orf[1], orf[2]) not in stops]


def exercise_6_3():
    print("\n6.3  ORFs that match no annotated gene in the same frame:")
    for threshold in (50, 100, 200):
        spurious = not_annotated(threshold)
        total = len(find_orfs(REGION.seq, threshold))
        print(f"       {threshold:3} amino acids: {total:3} ORFs, "
              f"{len(spurious):3} of them spurious")
    print("     At 50 the region is full of chance ORFs. At 200 almost")
    print("     none survive, but neither does lacI, a real gene of 360")
    print("     amino acids whose ORF here is shorter than the gene.")
    print("     For this region 100 is the reasonable choice, and the")
    print("     cost is the small real genes: anything under 100 amino")
    print("     acids is invisible, and bacteria have plenty of those.")


def exercise_6_4():
    print("\n6.4  The 13 mitochondrial coding sequences, NC_012920.1.")
    print("     Translate each with table 2, the vertebrate")
    print("     mitochondrial code, as Chapter 6 did for COX1.")
    print("     Several headers carry a note such as")
    print("     [transl_except=(pos:...,aa:TERM)], which says the record")
    print("     ends in T or TA rather than a full stop codon.")
    print("     The cell completes it: the messenger RNA is cut and a")
    print("     tail of A's is added, and the first two of those A's")
    print("     turn T or TA into TAA, a real stop. The stop codon is")
    print("     made after transcription, so it is not in the DNA.")
    if "--net" not in sys.argv:
        print("     Run with --net to download and translate them.")
        return
    import urllib.request
    url = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
           "?db=nuccore&id=NC_012920.1&rettype=fasta_cds_na"
           "&retmode=text")
    with urllib.request.urlopen(url, timeout=300) as reply:
        text = reply.read().decode()
    Path("mito_cds.fasta").write_text(text)
    short = 0
    for record in SeqIO.parse("mito_cds.fasta", "fasta"):
        protein = record.seq.translate(table=2)
        if not str(protein).endswith("*"):
            short += 1
        ending = "with a stop" if str(protein).endswith("*") else "short"
        print(f"       {record.id.split('|')[-1][:20]:22} "
              f"{len(record.seq):5} bases, ends {ending}")
    print(f"     {short} of them end without a full stop codon.")


def exercise_6_5():
    random.seed(1)
    length = 100_000
    made_up = Seq("".join(random.choice("ACGT") for _ in range(length)))
    chance = find_orfs(made_up, 100)
    real_region = find_orfs(REGION.seq, 100)
    per_kb = len(chance) / (length / 1000)
    print(f"\n6.5  {length:,} random bases hold {len(chance)} ORFs of")
    print(f"     100 amino acids or more, {per_kb:.2f} per kb.")
    print(f"     The lac region, {len(REGION.seq):,} real bases, holds "
          f"{len(real_region)},")
    real_per_kb = len(real_region) / (len(REGION.seq) / 1000)
    print(f"     {real_per_kb:.2f} per kb.")
    print("     Chance alone puts long ORFs in DNA that codes for")
    print("     nothing, so finding an ORF is not finding a gene. Real")
    print("     annotation uses the start signals, the codon usage, the")
    print("     similarity to known proteins and the evidence from RNA,")
    print("     which is why Chapter 24 took the annotation from NCBI")
    print("     rather than calling genes itself.")


if __name__ == "__main__":
    exercise_6_1()
    exercise_6_2()
    exercise_6_3()
    exercise_6_4()
    exercise_6_5()
