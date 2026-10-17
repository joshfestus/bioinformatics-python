"""Solutions to the exercises of Chapter 15.

Run from this folder: uv run ch15.py
Needs ../ch15/reads.bam, which Chapter 15's map.sh makes. Exercise 15.3
maps 100,000 pairs again and needs minimap2, samtools and the reads in
../data; it is skipped if any of those is missing.
"""

import gzip
import shutil
import statistics
import subprocess
import tempfile
from collections import Counter
from pathlib import Path

import matplotlib
import numpy as np
import pysam

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

CH15, DATA = Path("../ch15"), Path("../data")
BAM = CH15 / "reads.bam"
DELETION = (3_895_001, 3_901_452)  # the gap Chapter 15 found


def depth_array(path=CH15 / "depth.npy"):
    """Chapter 15's saved depth, or counted again if it is missing."""
    if path.exists():
        return np.load(path)
    with pysam.AlignmentFile(BAM) as bam:
        contig, length = bam.references[0], bam.lengths[0]
        counts = bam.count_coverage(contig, 0, length,
                                    quality_threshold=0)
    return np.array(counts).sum(axis=0)


def genes():
    """{name: (start, end)} from the REL606 annotation, one based."""
    found = {}
    for line in (CH15 / "rel606.gff3").read_text().splitlines():
        if line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 9 or parts[2] != "gene":
            continue
        for field in parts[8].split(";"):
            if field.startswith("Name="):
                found[field[5:]] = (int(parts[3]), int(parts[4]))
    return found


def exercise_15_1():
    qualities = Counter()
    with pysam.AlignmentFile(BAM) as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_secondary or read.is_supplementary:
                continue
            qualities[read.mapping_quality] += 1
    total = sum(qualities.values())
    figure, axes = plt.subplots(figsize=(6, 3.2))
    axes.bar(list(qualities), [qualities[q] for q in qualities],
             color="#555555")
    axes.set_xlabel("Mapping quality")
    axes.set_ylabel("Reads")
    axes.set_yscale("log")
    figure.savefig("ch15-mapq.png", dpi=110, bbox_inches="tight")
    plt.close(figure)
    print(f"15.1  {total:,} primary reads, "
          f"{len(qualities)} different mapping qualities.")
    print(f"       quality 60: {qualities[60]:,} reads "
          f"({qualities[60] / total:.1%})")
    print(f"       quality  0: {qualities[0]:,} reads "
          f"({qualities[0] / total:.1%})")
    middle = total - qualities[60] - qualities[0]
    print(f"       everything else: {middle:,} ({middle / total:.1%})")
    print("      Mapping quality is the aligner's confidence that the")
    print("      read came from where it put it, not that the read is")
    print("      correct. 60 means it found one clear best place. 0")
    print("      means it found two or more equally good ones and is")
    print("      telling you so, which in a bacterial genome usually")
    print("      means a repeat: the insertion sequences of Chapter")
    print("      15, the ribosomal RNA operons, and the transposase")
    print("      copies. Filtering on quality 20 throws those reads")
    print("      away, which is right for variant calling and wrong")
    print("      if the repeats are what you are studying. Saved as")
    print("      ch15-mapq.png.")


def exercise_15_2():
    lengths = []
    with pysam.AlignmentFile(BAM) as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_secondary or read.is_supplementary:
                continue
            if not read.is_read1 or read.is_unmapped:
                continue
            if read.mate_is_unmapped or read.next_reference_id != \
                    read.reference_id:
                continue
            if read.template_length:
                lengths.append(abs(read.template_length))
    proper = []
    with pysam.AlignmentFile(BAM) as bam:
        for read in bam.fetch(until_eof=True):
            if (read.is_read1 and read.is_proper_pair
                    and read.template_length):
                proper.append(abs(read.template_length))
    figure, axes = plt.subplots(figsize=(6, 3.2))
    axes.hist([x for x in lengths if x <= 2000], bins=100,
              color="#777777", edgecolor="#222222", linewidth=0.3)
    cutoff = max(proper) if proper else 0
    axes.axvline(cutoff, color="#222222", lw=1, ls=(0, (4, 3)))
    axes.set_xlabel("Fragment length (TLEN)")
    axes.set_ylabel("Pairs")
    figure.savefig("ch15-fragments.png", dpi=110, bbox_inches="tight")
    plt.close(figure)
    sensible = [x for x in lengths if x <= 2000]
    huge = [x for x in lengths if x > 100_000]
    circular = [x for x in huge if x > 4_000_000]
    print(f"\n15.2  {len(lengths):,} pairs on the same sequence.")
    print(f"       median {statistics.median(lengths):,.0f}, mean of")
    print(f"       the {len(sensible):,} under 2,000 bases "
          f"{statistics.mean(sensible):,.0f}")
    print(f"       the longest minimap2 still called properly paired: "
          f"{cutoff:,}")
    within = sum(1 for x in lengths if x <= cutoff)
    print(f"       {within / len(lengths):.1%} of pairs are at or below "
          "that length.")
    print(f"      The mean over every pair is "
          f"{statistics.mean(lengths):,.0f}, which is nonsense, and")
    print(f"      the reason is the {len(huge):,} pairs whose reads sit "
          "more than")
    print("      100,000 bases apart. They are two different faults.")
    print(f"      {len(circular)} of them stretch past 4 Mb, up to "
          f"{max(lengths):,},")
    print("      which is the chromosome's own circle: a fragment")
    print("      across the point where the circle was cut into a")
    print("      line has its reads at opposite ends of the file.")
    print(f"      The other {len(huge) - len(circular):,} are pairs "
          "where one read landed")
    print("      on the wrong copy of a repeat, which is the same")
    print("      ambiguity Exercise 15.1 saw as mapping quality zero.")
    print("      Either way, report the median for a distribution")
    print("      like this, or say what range you took the mean over.")
    print("      The shape itself is wide with a long tail, which is")
    print("      what a transposase made library looks like: the")
    print("      enzyme cuts where it likes rather than to a size.")
    print("      minimap2 sets its proper pair cutoff from the")
    print("      distribution it sees, so everything past about 820")
    print("      bases is called an unusual pair although most of")
    print("      those fragments are real. That is why Chapter 15's")
    print("      properly paired share is 68.55% while 97.25% of")
    print("      reads mapped. Saved as ch15-fragments.png.")


def exercise_15_3():
    print("\n15.3  Mapping the first 100,000 pairs instead of 600,000:")
    tools = all(shutil.which(tool) for tool in ("minimap2", "samtools"))
    reads = [DATA / "SRR2584863_1.fastq.gz",
             DATA / "SRR2584863_2.fastq.gz"]
    if not tools or not all(path.exists() for path in reads):
        print("      Needs minimap2, samtools and the reads in")
        print("      ../data; run data/fetch.py ch15 first.")
        return
    with tempfile.TemporaryDirectory() as work:
        small = []
        for source in reads:
            target = Path(work) / source.name.replace(".gz", "")
            with gzip.open(source, "rt") as handle, \
                    open(target, "w") as out:
                for number, line in enumerate(handle):
                    if number >= 400_000:  # four lines a read
                        break
                    out.write(line)
            small.append(target)
        bam = Path(work) / "small.bam"
        mapped = subprocess.run(
            ["minimap2", "-ax", "sr", "-t", "2",
             str(DATA / "rel606.fasta"), *map(str, small)],
            capture_output=True, check=True)
        subprocess.run(["samtools", "sort", "-o", str(bam)],
                       input=mapped.stdout, check=True,
                       capture_output=True)
        subprocess.run(["samtools", "index", str(bam)], check=True)
        with pysam.AlignmentFile(bam) as handle:
            contig, length = handle.references[0], handle.lengths[0]
            counts = handle.count_coverage(contig, 0, length,
                                           quality_threshold=0)
        small_depth = np.array(counts).sum(axis=0)
    full = depth_array()
    first, last = DELETION
    print(f"       600,000 pairs: mean depth {full.mean():.1f}")
    print(f"       100,000 pairs: mean depth {small_depth.mean():.1f}")
    print(f"       the ratio is {full.mean() / small_depth.mean():.1f}, "
          "as the read count says it should be.")
    gap = small_depth[first - 1:last]
    print(f"       Across {first:,} to {last:,}, the deleted region,")
    print(f"       {(gap == 0).sum():,} of its {len(gap):,} bases still")
    print(f"       have no read at all, and the deepest point in it is")
    print(f"       {gap.max()}.")
    zeros = int((small_depth == 0).sum())
    print(f"       In the whole genome {zeros:,} bases now have no")
    print(f"       read at all, against {int((full == 0).sum()):,} "
          "before.")
    print("      A deletion is still a deletion at a sixth of the")
    print("      depth, because zero reads is zero reads. What falls")
    print("      apart at low depth is everything that needs counting:")
    print("      a variant seen in three reads of five is noise, and")
    print("      Chapter 16's depth filter would drop most of the")
    print("      genome. Deep sequencing buys confidence in small")
    print("      differences, not the ability to see large ones.")


def exercise_15_4():
    depth = depth_array()
    where = genes().get("lacZ")
    if where is None:
        print("\n15.4  No lacZ in the annotation.")
        return
    first, last = where
    across = depth[first - 1:last]
    figure, axes = plt.subplots(figsize=(6, 3.2))
    axes.plot(range(first, last + 1), across, color="#333333", lw=0.6)
    axes.axhline(depth.mean(), color="#bbbbbb", lw=0.8, ls=(0, (4, 3)))
    axes.set_xlabel("Position in REL606")
    axes.set_ylabel("Depth")
    figure.savefig("ch15-lacz.png", dpi=110, bbox_inches="tight")
    plt.close(figure)
    print(f"\n15.4  lacZ runs {first:,} to {last:,}, "
          f"{last - first + 1:,} bases.")
    print(f"       depth across it: mean {across.mean():.1f}, lowest "
          f"{across.min()}, highest {across.max()}")
    print(f"       the genome's mean is {depth.mean():.1f}, so lacZ "
          f"sits at {across.mean() / depth.mean():.2f} of it.")
    print(f"       bases with no read at all: {(across == 0).sum()}")
    print("      The clone still has lacZ, covered evenly at about the")
    print("      genome's average depth from one end to the other.")
    print("      That is what present looks like, and it is worth")
    print("      seeing beside the deleted region: there the depth")
    print("      goes to nothing for six thousand bases and comes")
    print("      back. Depth is how you tell a gene that is there")
    print("      from a gene that is gone, and neither needs a")
    print("      variant caller. Saved as ch15-lacz.png.")


def is150_copies(gap=200):
    """Merged intervals of the CDS features called IS150 transposase.

    Each IS150 element carries two overlapping reading frames, so the
    annotation has two CDS lines per copy; merging the ones that touch
    turns nine CDS features into the elements they belong to.
    """
    pieces = []
    for line in (CH15 / "rel606.gff3").read_text().splitlines():
        if line.startswith("#") or "IS150" not in line:
            continue
        parts = line.split("\t")
        if len(parts) > 4 and parts[2] == "CDS":
            pieces.append((int(parts[3]), int(parts[4])))
    merged = []
    for start, end in sorted(pieces):
        if merged and start - merged[-1][1] <= gap:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return [tuple(piece) for piece in merged], len(pieces)


def exercise_15_5():
    depth = depth_array()
    average = depth.mean()
    copies, features = is150_copies()
    print(f"\n15.5  {features} CDS features are called IS150 family")
    print(f"      transposase, which merge into {len(copies)} elements:")
    ratios = []
    for start, end in copies:
        across = depth[start - 1:end]
        ratio = across.mean() / average
        ratios.append(ratio)
        print(f"       {start:>9,} to {end:>9,}: mean depth "
              f"{across.mean():6.1f}, {ratio:4.1f} times the genome")
    estimate = statistics.mean(ratios) * len(copies)
    print(f"      The copies average {statistics.mean(ratios):.1f} "
          "times the genome's depth,")
    print(f"      so the clone holds roughly {estimate:.0f} of them")
    print(f"      against the reference's {len(copies)}.")
    print("      Reads from every copy in the clone pile onto the few")
    print("      the reference has, because the copies are identical")
    print("      and the mapper cannot tell them apart. That is why")
    print("      their mapping quality is zero, and it is what makes")
    print("      the depth a copy count at all.")
    print("      The estimate assumes four things, and each can be")
    print("      wrong. That the sequencing covered the genome evenly,")
    print("      which Chapter 15 checked and it does. That every copy")
    print("      in the clone is close enough to a reference copy to")
    print("      map there. That the copies are full length, so depth")
    print("      scales with their number. And that nothing else in")
    print("      the genome resembles IS150 enough to donate reads.")
    print("      A number like this has a direction and an order of")
    print("      magnitude, not a value: to count the copies you would")
    print("      assemble the clone's own genome and look.")


if __name__ == "__main__":
    exercise_15_1()
    exercise_15_2()
    exercise_15_3()
    exercise_15_4()
    exercise_15_5()
