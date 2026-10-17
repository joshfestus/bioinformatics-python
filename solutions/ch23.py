"""Solutions to the exercises of Project One, Chapter 23.

Run from this folder: uv run ch23.py
Reads the project's own files in ../projects/virus. Exercise 23.1 needs
the network; the rest work from what the project already downloaded.
"""

import csv
import statistics
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib
from Bio import SeqIO
from Bio.Seq import Seq

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

VIRUS = Path("../projects/virus").resolve()
BASES = "ACGT"


def read(name):
    with open(VIRUS / name) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def genomes():
    return {record.id: str(record.seq).upper()
            for record in SeqIO.parse(VIRUS / "aligned.fasta", "fasta")}


def exercise_23_1():
    print("23.1  Adding a fourteenth lineage.")
    print("      One line in choose.py's LINEAGES, giving the Pango")
    print("      name and the WHO label if it has one:\n")
    print('        "XEC": "Omicron",\n')
    print("      Then rerun choose.py, download.py and the pipeline.")
    print("      XEC was circulating in late 2024 and early 2025, so")
    print("      it should sit beyond KP.3 on the clock, with more")
    print("      substitutions than anything else in the sample.")
    print("      Two things to watch. NCBI adds genomes every week,")
    print("      so your sample will not match the book's even for")
    print("      the thirteen lineages it already has, which is why")
    print("      samples.tsv is committed. And a lineage with few")
    print("      complete genomes will come back short of the ten")
    print("      the chooser asks for, which is fine and worth")
    print("      reporting rather than hiding.")


def exercise_23_2():
    print("\n23.2  Loosening max_unknown from 0.01 to 0.05.")
    kept = {}
    for name in ("clean.fasta",):
        kept[name] = sum(1 for _ in SeqIO.parse(VIRUS / name, "fasta"))
    samples = {row["accession"]: row for row in read("samples.tsv")}
    counts = []
    for record in SeqIO.parse(VIRUS / "genomes.fasta", "fasta"):
        sequence = str(record.seq).upper()
        called = sum(sequence.count(base) for base in BASES)
        unknown = 1 - called / len(sequence)
        counts.append((record.id, len(sequence), unknown))
    for threshold in (0.01, 0.05):
        survivors = [row for row in counts
                     if row[1] >= 29_000 and row[2] <= threshold]
        print(f"       max_unknown {threshold:.2f}: "
              f"{len(survivors)} of {len(counts)} genomes")
    extra = [row for row in counts
             if row[1] >= 29_000 and 0.01 < row[2] <= 0.05]
    print(f"       the {len(extra)} genomes this lets back in:")
    for accession, length, unknown in sorted(extra,
                                             key=lambda r: -r[2])[:10]:
        lineage = samples.get(accession, {}).get("lineage", "?")
        print(f"         {accession:12} {lineage:10} {unknown:.1%} "
              "uncalled")
    print("      They are the ones the stricter rule was written to")
    print("      exclude: a genome that is five percent N has 1,500")
    print("      positions where nothing can be said. Those positions")
    print("      do not become wrong, they become invisible, so the")
    print("      count of differences per genome falls and the clock")
    print("      flattens slightly. The honest way to use the looser")
    print("      threshold is to report coverage alongside the")
    print("      counts, which is what the project's own quality step")
    print("      does from the other direction.")


def exercise_23_3():
    print("\n23.3  Marking every substitution synonymous or not.")
    annotation = SeqIO.read(VIRUS / "reference.gb", "genbank")
    genes = {}
    for feature in annotation.features:
        if feature.type != "CDS":
            continue
        name = feature.qualifiers["gene"][0]
        product = feature.qualifiers["product"][0]
        if name == "ORF1ab" and "ORF1a " in product:
            continue
        positions = []
        for part in feature.location.parts:
            positions += list(range(int(part.start), int(part.end)))
        genes[name] = positions
    where = {}
    for name, positions in genes.items():
        for index, position in enumerate(positions):
            where.setdefault(position, (name, index))
    sequences = genomes()
    reference = sequences.pop(annotation.id)
    counts = defaultdict(lambda: [0, 0])
    for sequence in sequences.values():
        first = len(sequence) - len(sequence.lstrip("-"))
        last = len(sequence.rstrip("-")) - 1
        for position, (ours, theirs) in enumerate(zip(reference,
                                                      sequence)):
            if not first <= position <= last or theirs not in BASES:
                continue
            if ours == theirs or position not in where:
                continue
            name, index = where[position]
            codon_start = index - index % 3
            codon_positions = genes[name][codon_start:codon_start + 3]
            if len(codon_positions) < 3:
                continue
            old = "".join(reference[p] for p in codon_positions)
            new = list(old)
            new[index % 3] = theirs
            new = "".join(new)
            if any(base not in BASES for base in old + new):
                continue
            same = Seq(old).translate() == Seq(new).translate()
            counts[name][0 if same else 1] += 1
    print(f"       {'gene':10} {'synonymous':>11} {'changes amino':>14} "
          f"{'share changing':>15}")
    def changing(item):
        return -(item[1][1] / max(sum(item[1]), 1))

    for name, (same, changed) in sorted(counts.items(), key=changing):
        total = same + changed
        if total < 20:
            continue
        print(f"       {name:10} {same:>11,} {changed:>14,} "
              f"{changed / total:>14.0%}")
    envelope = counts.get("E", [0, 0])
    spike = counts.get("S", [0, 0])
    orf1ab = counts.get("ORF1ab", [0, 0])
    print(f"      E tops the table at {envelope[1]}/{sum(envelope)}, and")
    print("      that is the small denominator trap again: the")
    print("      envelope protein is 228 bases long and barely")
    print("      changes, so a handful of changes decide its share.")
    print(f"      The spike is the meaningful one, "
          f"{spike[1] / sum(spike):.0%} of "
          f"{sum(spike):,} changes")
    print("      altering the protein. That is what selection on a")
    print("      protein looks like: a change that alters nothing is")
    print("      no use to a virus escaping antibodies.")
    print(f"      ORF1ab is the other end, {orf1ab[0]:,} synonymous")
    print(f"      changes of {sum(orf1ab):,}, the most of any gene.")
    print("      It is two thirds of the genome and encodes the")
    print("      replication machinery, where most of what happens")
    print("      is drift that selection never sees.")


def exercise_23_4():
    print("\n23.4  One mutation through time.")
    wanted = "E484K"
    rows = [row for row in read("mutations.tsv")
            if row["gene"] == "S" and row["protein"] == wanted]
    print(f"       {wanted}: {len(rows)} genomes carry it\n")
    print(f"       {'date':12} {'lineage':10} {'accession':12}")
    for row in sorted(rows, key=lambda r: r["date"])[:12]:
        print(f"       {row['date']:12} {row['lineage']:10} "
              f"{row['accession']:12}")
    if len(rows) > 12:
        print(f"       ... and {len(rows) - 12} more")
    lineages = Counter(row["lineage"] for row in rows)
    print(f"\n       by lineage: "
          + ", ".join(f"{name} {count}" for name, count in
                      lineages.most_common()))
    print("      E484K appears first in Beta, then in Gamma, and")
    print("      again in later Omicron lineages. It is not inherited")
    print("      down one line: the same change arose more than once")
    print("      in lineages that had already parted. That is")
    print("      convergent evolution, and it is the strongest")
    print("      evidence a sequence can give that a change is")
    print("      useful to the virus rather than an accident.")


def exercise_23_5():
    print("\n23.5  The clock fitted inside single lineages.")
    sys.path.insert(0, str(VIRUS))
    from timeline import fit
    rows = [dict(row, days=int(row["days"]),
                 substitutions=int(row["substitutions"]))
            for row in read("clock.tsv")]
    everything = fit([r["days"] for r in rows],
                     [r["substitutions"] for r in rows])
    print(f"       across the whole sample: "
          f"{everything[0] * 365:.1f} substitutions a year")
    print(f"\n       {'lineage':10} {'genomes':>8} {'span, days':>11} "
          f"{'rate a year':>12}")
    by_lineage = defaultdict(list)
    for row in rows:
        by_lineage[row["lineage"]].append(row)
    for lineage, found in by_lineage.items():
        days = [r["days"] for r in found]
        if len(found) < 4 or max(days) - min(days) < 60:
            continue
        rate, _ = fit(days, [r["substitutions"] for r in found])
        print(f"       {lineage:10} {len(found):>8} "
              f"{max(days) - min(days):>11} {rate * 365:>12.1f}")
    print("      The rates inside lineages are all over the place,")
    print("      some far above the overall rate and some negative,")
    print("      and none of them means anything. Each lineage spans")
    print("      a few months and ten genomes, so the fit is reading")
    print("      noise: which particular genomes were sequenced")
    print("      matters more than when they were collected.")
    print("      The clock works across the whole sample because")
    print("      four years and a hundred genomes are enough for the")
    print("      signal to beat the noise. A rate is a statement")
    print("      about a range, and quoting one from a short range is")
    print("      how people end up reporting that a virus evolves")
    print("      backwards.")


def exercise_23_6():
    rows = [row for row in read("mutations.tsv")
            if row["gene"] == "S" and row["protein"]
            and row["protein"] != "del"]
    counts = Counter(row["protein"] for row in rows).most_common(10)
    names = [name for name, _ in counts][::-1]
    values = [count for _, count in counts][::-1]
    figure, axes = plt.subplots(figsize=(5.4, 3.2))
    axes.barh(names, values, color="#555555", edgecolor="#222222",
              linewidth=0.5)
    axes.set_xlabel("Genomes carrying the change")
    figure.savefig("ch23-spike-changes.png", dpi=110,
                   bbox_inches="tight")
    plt.close(figure)
    print("\n23.6  The ten commonest spike changes, as a bar chart:")
    for name, count in counts:
        print(f"       {name:8} {count:>4} genomes")
    print("      Horizontal bars, gray fill, no color: the same")
    print("      choices Chapter 5 made, and for the same reason,")
    print("      which is that the book prints in black and white.")
    print("      Saved as ch23-spike-changes.png.")


def exercise_23_7():
    print("\n23.7  Making the pipeline fail on purpose.")
    print("      Deleting reference.gb and running snakemake gives a")
    print("      MissingInputException naming the rule and the file:")
    print("      Snakemake checks every input before it runs anything,")
    print("      so the pipeline stops before mutations.py is called.\n")
    print("      Run mutations.py on its own and the error is uglier,")
    print("      a FileNotFoundError traceback through Biopython. The")
    print("      fix is to say what is missing and why:\n")
    fix = '''ANNOTATION = Path("reference.gb")

if not ANNOTATION.exists():
    raise SystemExit(
        f"{ANNOTATION} is missing; run download.py first")'''
    for line in fix.splitlines():
        print(f"        {line}")
    print("\n      A traceback tells a programmer where the program")
    print("      broke. A message like that one tells anybody what")
    print("      to do next, which is what a reader of your")
    print("      repository needs when a download failed last week")
    print("      and they have forgotten.")


if __name__ == "__main__":
    exercise_23_1()
    exercise_23_2()
    exercise_23_3()
    exercise_23_4()
    exercise_23_5()
    exercise_23_6()
    exercise_23_7()
