"""When each change first appears, and how fast changes accumulate.

Every genome carries a collection date, so the sample is a record of the
virus through time. Two questions follow: how many changes a genome picked
up per year, and when each famous spike change was first collected.

uv run timeline.py         reads mutations.tsv, writes clock.tsv
"""

import csv
from collections import Counter, defaultdict
from datetime import date

from Bio import SeqIO

WATCH = ["D614G", "N501Y", "P681H", "P681R", "L452R", "T478K", "E484K",
         "K417N", "Q498R", "L455S"]
START = date(2019, 12, 24)  # the collection date of the reference genome


def read_days(text):
    year, month, day = (int(part) for part in text.split("-"))
    return (date(year, month, day) - START).days


def main():
    with open("mutations.tsv") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    with open("samples.tsv") as f:
        samples = list(csv.DictReader(f, delimiter="\t"))
    # the genomes the analysis actually ran on, named by the cleaned file,
    # not by samples.tsv, which still lists the ones clean.py dropped
    analyzed = {record.id
                for record in SeqIO.parse("clean.fasta", "fasta")}

    clock = dated_counts(rows, samples, analyzed)
    with open("clock.tsv", "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=clock[0], delimiter="\t")
        writer.writeheader()
        writer.writerows(clock)

    rate, intercept = fit([row["days"] for row in clock],
                          [row["substitutions"] for row in clock])
    days = [row["days"] for row in clock]
    print(f"{len(clock)} genomes, collected over "
          f"{max(days) - min(days)}"
          " days")
    print(f"substitutions per year: {rate * 365:.1f}")
    print("at the reference's own date the fit expects "
          f"{intercept:.1f}\n")

    first_seen(rows)


def dated_counts(rows, samples, analyzed):
    """One row per genome: date, age in days, substitution count."""
    substitutions = Counter()
    for row in rows:
        if row["sample"] != "-":
            substitutions[row["accession"]] += 1
    clock = []
    for sample in samples:
        accession = sample["accession"]
        if accession not in analyzed:
            continue  # clean.py dropped it
        clock.append({
            "accession": accession,
            "lineage": sample["lineage"],
            "date": sample["date"],
            "days": read_days(sample["date"]),
            "substitutions": substitutions.get(accession, 0),
        })
    return clock


def first_seen(rows):
    """The earliest collection date carrying each watched spike change."""
    first = defaultdict(list)
    for row in rows:
        if row["gene"] == "S" and row["protein"] in WATCH:
            first[row["protein"]].append((row["date"], row["lineage"]))
    print(f"{'change':8} {'first seen':12} {'in':10} {'genomes':>8}")
    for change in WATCH:
        seen = sorted(first[change])
        if not seen:
            continue
        print(f"{change:8} {seen[0][0]:12} {seen[0][1]:10} "
              f"{len(seen):>8}")


def fit(x, y):
    """Least squares slope and intercept, written out in full."""
    n = len(x)
    mean_x, mean_y = sum(x) / n, sum(y) / n
    top = sum((a - mean_x) * (b - mean_y) for a, b in zip(x, y))
    bottom = sum((a - mean_x) ** 2 for a in x)
    slope = top / bottom
    return slope, mean_y - slope * mean_x


if __name__ == "__main__":
    main()
