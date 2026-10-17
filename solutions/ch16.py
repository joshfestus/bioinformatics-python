"""Solutions to the exercises of Chapter 16.

Run from this folder: uv run ch16.py
Needs ../ch15/reads.bam and the VCF files Chapter 16 wrote.
"""

import subprocess
from collections import Counter
from pathlib import Path

import pysam
from Bio import SeqIO
from Bio.Seq import Seq

CH15, CH16, DATA = Path("../ch15"), Path("../ch16"), Path("../data")
BAM, CALLS = CH15 / "reads.bam", CH16 / "haploid.vcf.gz"
PYKF = 1_733_343  # the D127N mutation Chapter 16 followed
IS1 = 263_235


def confident(min_quality=30, min_depth=10):
    """Chapter 16's confident calls: quality 30 or more, depth 10."""
    kept = []
    with pysam.VariantFile(CALLS) as vcf:
        for record in vcf:
            depth = record.info.get("DP", 0)
            if record.qual and record.qual >= min_quality \
                    and depth >= min_depth:
                kept.append(record)
    return kept


def exercise_16_1():
    print("16.1  samtools tview prints the pileup as text:\n")
    print("       samtools tview ../ch15/reads.bam "
          "../data/rel606.fasta \\")
    print(f"           -d T -p NC_012967.1:{PYKF}\n")
    with pysam.AlignmentFile(BAM) as bam:
        counts = Counter()
        for column in bam.pileup("NC_012967.1", PYKF - 1, PYKF,
                                 truncate=True, min_base_quality=0):
            for read in column.pileups:
                if read.query_position is not None:
                    counts[read.alignment.query_sequence[
                        read.query_position]] += 1
    reference = SeqIO.read(DATA / "rel606.fasta", "fasta").seq
    print(f"      At {PYKF:,} the reference has "
          f"{reference[PYKF - 1]}, and the reads carry:")
    for base, count in counts.most_common():
        print(f"        {base}: {count:3} reads "
              f"({count / sum(counts.values()):.0%})")
    print("      Every read agrees on the same change, which is what")
    print("      a real mutation in a clonal population looks like.")
    print("      The picture tview draws is worth the trouble: a")
    print("      column of identical letters is convincing in a way")
    print("      that a quality score is not, and Chapter 16's rule")
    print("      is never to trust a call you have not looked at.")


def pileup_bases(bam, position):
    """The base each read shows at a position, and its mapping quality."""
    bases, qualities = Counter(), Counter()
    for column in bam.pileup("NC_012967.1", position - 1, position,
                             truncate=True, min_base_quality=0):
        for read in column.pileups:
            qualities[read.alignment.mapping_quality] += 1
            if read.query_position is not None:
                bases[read.alignment.query_sequence[
                    read.query_position]] += 1
    return bases, qualities


def exercise_16_2():
    calls = confident()
    print(f"\n16.2  The alternate allele's share in each of the "
          f"{len(calls)} calls:")
    low = []
    with pysam.AlignmentFile(BAM) as bam:
        for record in calls:
            alt = record.alts[0]
            bases, qualities = pileup_bases(bam, record.pos)
            total = sum(bases.values())
            share = bases[alt[0]] / total if total else 0
            ambiguous = (qualities[0] / sum(qualities.values())
                         if qualities else 0)
            if share < 0.9:
                low.append((record.pos, record.ref, alt, share, total,
                            ambiguous))
    print(f"       {len(calls) - len(low)} calls are at 90 percent or "
          f"more. {len(low)} are below,")
    print("       shown with the share of their reads the aligner")
    print("       could not place uniquely:")
    print(f"         {'position':>11} {'change':8} {'alt':>6} "
          f"{'reads':>6} {'quality 0':>10}")
    for position, ref, alt, share, total, ambiguous in low:
        print(f"         {position:>11,} {ref + '>' + alt:8} "
              f"{share:>6.0%} {total:>6} {ambiguous:>10.0%}")
    worst = max(row[5] for row in low)
    least = min(row[5] for row in low)
    with pysam.AlignmentFile(BAM) as bam:
        _, clean = pileup_bases(bam, PYKF)
    clean_share = clean[0] / sum(clean.values())
    print("      Every one of them sits where the aligner is unsure.")
    print(f"      Between {least:.0%} and {worst:.0%} of the reads "
          "under these calls")
    print("      have mapping quality 0, against "
          f"{clean_share:.0%} under the pykF")
    print("      call of Exercise 16.1, where every read carries the")
    print("      change.")
    print("      That explains the share directly. These positions")
    print("      are repeats, so reads from several copies of the")
    print("      same sequence pile onto one place. Some copies")
    print("      carry the difference and some do not, and the")
    print("      mixture comes out as an allele share near a half.")
    print("      In a clone of a haploid organism there is no such")
    print("      thing as half an allele: a share near 50 percent is")
    print("      not a heterozygote, it is two places being read as")
    print("      one. Exercise 16.3 looks at the clearest of them.")


def exercise_16_3():
    print(f"\n16.3  The reads at {IS1:,}, by mapping quality:")
    qualities = Counter()
    with pysam.AlignmentFile(BAM) as bam:
        for read in bam.fetch("NC_012967.1", IS1 - 1, IS1):
            qualities[read.mapping_quality] += 1
    total = sum(qualities.values())
    zero = qualities[0]
    print(f"       {total} reads cover it, {zero} of them at quality 0 "
          f"({zero / total:.0%})")
    print(f"       the rest: "
          f"{dict(sorted(q for q in qualities.items() if q[0]))}")
    print("      Mapping quality 0 means the aligner found more than")
    print("      one equally good place for the read, so it cannot")
    print("      say this read came from here. When most of the reads")
    print("      under a call say that, the call is about a repeat")
    print("      rather than about a position: IS1 has several copies")
    print("      in this genome, and a difference between the copies")
    print("      shows up as a variant at whichever copy the reads")
    print("      were piled onto. Chapter 16 found these calls sit in")
    print("      genes that the annotation says are frameshifted,")
    print("      which is the same story from the other direction.")


def exercise_16_4():
    print("\n16.4  HBB's pathogenic ClinVar variants by consequence:")
    kinds = Counter()
    with open(CH16 / "hbb_clinvar.vcf") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            info = line.split("\t")[7]
            fields = dict(part.split("=", 1) for part in info.split(";")
                          if "=" in part)
            significance = fields.get("CLNSIG", "")
            if "Pathogenic" not in significance:
                continue
            consequence = fields.get("MC", "")
            for entry in consequence.split(","):
                if "|" in entry:
                    kinds[entry.split("|")[1]] += 1
    for kind, count in kinds.most_common():
        print(f"       {kind:28} {count:4}")
    commonest = kinds.most_common(1)[0]
    print(f"      The commonest is {commonest[0]}, with "
          f"{commonest[1]}, ahead of")
    print(f"      missense at {kinds['missense_variant']} and nonsense "
          f"at {kinds['nonsense']}.")
    print("      That ordering is the disease. Most pathogenic HBB")
    print("      variants cause beta thalassemia, where the problem")
    print("      is that too little beta globin is made, and the")
    print("      surest way to make none is to shift the reading")
    print("      frame or stop the protein early. Sickle cell")
    print("      disease is the famous exception: one missense change")
    print("      that makes a protein which folds and carries oxygen")
    print("      and then sticks to its neighbors.")
    print("      Counting by consequence turns a list of 1,784")
    print("      records into a statement about mechanism, which is")
    print("      what the MC field is for.")


def cds_for(gene):
    """(start, end, strand) of a gene's CDS, from the annotation."""
    for line in (CH15 / "rel606.gff3").read_text().splitlines():
        if line.startswith("#") or f"gene={gene};" not in line + ";":
            continue
        parts = line.split("\t")
        if len(parts) > 6 and parts[2] == "CDS":
            return int(parts[3]), int(parts[4]), parts[6]
    raise KeyError(gene)


def codon_change(reference, start, end, strand, position, alt):
    """The codon before and after a single base change."""
    gene = reference[start - 1:end]
    offset = position - start
    if strand == "-":
        gene = gene.reverse_complement()
        offset = end - position
        alt = str(Seq(alt).reverse_complement())
    number = offset // 3
    first = number * 3
    before = str(gene[first:first + 3])
    after = before[:offset % 3] + alt + before[offset % 3 + 1:]
    return number + 1, before, after


def exercise_16_5():
    reference = SeqIO.read(DATA / "rel606.fasta", "fasta").seq
    start, end, strand = cds_for("pykF")
    number, before, after = codon_change(reference, start, end, strand,
                                         PYKF, "A")
    print("\n16.5  The codon before and after, for pykF's change:")
    print(f"       pykF runs {start:,} to {end:,} on the "
          f"{strand} strand")
    print(f"       codon {number}: {before} to {after}")
    print(f"       amino acid: {Seq(before).translate()}{number}"
          f"{Seq(after).translate()}")
    print("      Checked by hand. The reference base at "
          f"{PYKF:,} is {reference[PYKF - 1]},")
    print(f"      the call says this clone has A, and {PYKF:,} minus")
    print(f"      the gene's start {start:,} is {PYKF - start}, which "
          f"is {PYKF - start} / 3")
    print(f"      = codon {number} with the change in its first base.")
    print(f"      GAT to AAT is aspartic acid to asparagine, the")
    print("      D127N that Chapter 16 reports and that Lenski's")
    print("      populations found again and again.")
    print("      Reading the gene's coordinates and strand from the")
    print("      annotation rather than writing them into the program")
    print("      is not tidiness. Typing them by hand put this answer")
    print("      in the wrong gene and the wrong frame the first time")
    print("      it was written, and the codon printed beside the")
    print("      amino acid is what showed it.")


def exercise_16_extra():
    print("\n      For reference, samtools agrees about the depth here:")
    done = subprocess.run(
        ["samtools", "depth", "-r", f"NC_012967.1:{PYKF}-{PYKF}",
         str(BAM)], capture_output=True, text=True)
    print(f"        {done.stdout.strip() or 'samtools not available'}")


if __name__ == "__main__":
    exercise_16_1()
    exercise_16_2()
    exercise_16_3()
    exercise_16_4()
    exercise_16_5()
    exercise_16_extra()
