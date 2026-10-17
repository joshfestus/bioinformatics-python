"""Solutions to the exercises of Chapter 14.

Run from this folder: uv run ch14.py
Reads ../ch14/reads.fastq.gz, the 48,000 reads Chapter 14 works on.
Figures are saved here as PNG files.
"""

import gzip
import statistics
from pathlib import Path

import matplotlib
from Bio import SeqIO

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

CH14 = Path("../ch14")
ADAPTER = "CTGTCTCTTATACACATCT"  # Nextera, the kit this run used
WINDOW, MIN_LENGTH = 4, 50


def reads():
    with gzip.open(CH14 / "reads.fastq.gz", "rt") as handle:
        yield from SeqIO.parse(handle, "fastq")


def exercise_14_1():
    columns = []
    for record in reads():
        qualities = record.letter_annotations["phred_quality"]
        while len(columns) < len(qualities):
            columns.append([])
        for position, quality in enumerate(qualities):
            columns[position].append(quality)
    means = [statistics.mean(column) for column in columns]
    medians = [statistics.median(column) for column in columns]
    first_mean = next((i for i, value in enumerate(means) if value < 30),
                      None)
    first_median = next((i for i, value in enumerate(medians)
                         if value < 30), None)
    print(f"14.1  {len(columns)} positions over "
          f"{len(columns[0]):,} reads.")
    print(f"       mean first below 30 at position "
          f"{first_mean + 1 if first_mean is not None else 'never'}")
    print(f"       median first below 30 at position "
          f"{first_median + 1 if first_median is not None else 'never'}")
    print(f"       at the last position: mean {means[-1]:.1f}, "
          f"median {medians[-1]:.0f}")
    figure, axes = plt.subplots(figsize=(6, 3.2))
    axes.plot(range(1, len(means) + 1), means, color="#333333",
              label="mean")
    axes.plot(range(1, len(medians) + 1), medians, color="#888888",
              label="median")
    axes.axhline(30, color="#bbbbbb", lw=0.8, ls=(0, (4, 3)))
    axes.set_xlabel("Position in the read")
    axes.set_ylabel("Quality")
    axes.legend(frameon=False)
    figure.savefig("ch14-median.png", dpi=110, bbox_inches="tight")
    plt.close(figure)
    print("      The median holds up longer than the mean, because the")
    print("      mean is dragged down by a minority of bad reads while")
    print("      the median follows the typical read. Both are worth")
    print("      having: the mean tells you what trimming will cost,")
    print("      the median what a good read looks like. Saved as")
    print("      ch14-median.png.")


def exercise_14_2():
    starts = []
    for record in reads():
        at = str(record.seq).find(ADAPTER)
        if at != -1:
            starts.append(at)
    figure, axes = plt.subplots(figsize=(6, 3.2))
    axes.hist(starts, bins=50, color="#777777", edgecolor="#222222",
              linewidth=0.4)
    axes.set_xlabel("Position where the Nextera adapter begins")
    axes.set_ylabel("Reads")
    figure.savefig("ch14-adapter-start.png", dpi=110,
                   bbox_inches="tight")
    plt.close(figure)
    print(f"\n14.2  {len(starts):,} reads contain the full adapter.")
    print(f"       It starts anywhere from {min(starts)} to "
          f"{max(starts)}, median {statistics.median(starts):.0f}.")
    print("      The adapter appears when the fragment is shorter than")
    print("      the read, so the sequencer runs off the end of the")
    print("      DNA and into the adapter. Where it starts is the")
    print("      length of that fragment. The spread of these starts")
    print("      is the short tail of the library's fragment length")
    print("      distribution, and Chapter 15 measures the rest of it")
    print("      from the mapped pairs. Saved as ch14-adapter-start.png.")


def trimmed(record, min_quality, partial=0):
    """Chapter 14's trim, with the threshold as an argument."""
    sequence = str(record.seq)
    cut = sequence.find(ADAPTER)
    if cut == -1 and partial:
        for length in range(len(ADAPTER) - 1, partial - 1, -1):
            if sequence.endswith(ADAPTER[:length]):
                cut = len(sequence) - length
                break
    if cut != -1:
        record = record[:cut]
    qualities = record.letter_annotations["phred_quality"]
    end = len(qualities)
    for start in range(len(qualities) - WINDOW + 1):
        if sum(qualities[start:start + WINDOW]) / WINDOW < min_quality:
            end = start
            break
    return record[:end]


def exercise_14_3():
    print("\n14.3  The same trimming at three thresholds:")
    print(f"       {'quality':>8} {'reads kept':>11} {'bases kept':>12} "
          f"{'mean length':>12}")
    for threshold in (15, 20, 30):
        kept, bases, total_in = 0, 0, 0
        for record in reads():
            total_in += len(record)
            cut = trimmed(record, threshold)
            if len(cut) >= MIN_LENGTH:
                kept += 1
                bases += len(cut)
        print(f"       {threshold:>8} {kept:>11,} {bases:>12,} "
              f"{bases / kept:>12.1f}")
    print("      A gentler threshold keeps more bases and more of the")
    print("      poor ends with them; a stricter one throws away good")
    print("      sequence to be sure of what remains. At 30 the reads")
    print("      are shorter and cleaner, and a mapper will place")
    print("      fewer of them. Chapter 14 chose 20 because mappers")
    print("      handle a few wrong bases better than they handle")
    print("      missing sequence.")


def exercise_14_4():
    print("\n14.4  Adapters that run off the end of the read:")
    full = partial_only = 0
    for record in reads():
        sequence = str(record.seq)
        if ADAPTER in sequence:
            full += 1
            continue
        for length in range(len(ADAPTER) - 1, 7, -1):
            if sequence.endswith(ADAPTER[:length]):
                partial_only += 1
                break
    print(f"       {full:,} reads hold the whole adapter.")
    print(f"       {partial_only:,} more end in at least eight letters")
    print("       of it, which the chapter's trim leaves in place.")
    print(f"       That is {partial_only / (full + partial_only):.1%} "
          "more reads trimmed.")
    print("      Eight letters is a judgment: short enough to catch")
    print("      the real tails, long enough that it happens by")
    print("      chance about once in 65,000 reads at each length.")
    print("      Go much shorter and you start cutting real sequence")
    print("      off the end of every read, which is worse than")
    print("      leaving a few adapter bases in.")


def exercise_14_5():
    print("\n14.5  The FastQC report on the untrimmed reads fails per")
    print("      base sequence content, and the plot shows why: the")
    print("      four base lines are far apart over the first ten or")
    print("      so positions and level off afterwards.")
    print("      This is not a fault in the run. The library was made")
    print("      with the Nextera kit, where a transposase cuts the")
    print("      DNA and adds the adapters in one step, and that")
    print("      enzyme prefers certain sequences at its cut site. So")
    print("      the first bases of every fragment are biased by the")
    print("      chemistry that made it.")
    print("      FastQC has no way to know that, so it reports a")
    print("      failure, and Chapter 14 keeps the reads anyway. A")
    print("      tool's red cross is evidence, not a verdict: the")
    print("      question is always whether the pattern it found")
    print("      matters for what you are about to do, and a known")
    print("      bias in the first ten bases does not stop a mapper.")


if __name__ == "__main__":
    exercise_14_1()
    exercise_14_2()
    exercise_14_3()
    exercise_14_4()
    exercise_14_5()
