"""Write the project's findings as a report, from the files it produced.

uv run report.py           writes report.md
"""

import csv
import statistics
from collections import Counter, defaultdict

from Bio import SeqIO

WHO = {"B.1.1.7": "Alpha", "B.1.351": "Beta", "P.1": "Gamma",
       "B.1.617.2": "Delta", "BA.1": "Omicron", "BA.2": "Omicron",
       "BA.5.2.1": "Omicron", "XBB.1.5": "Omicron", "JN.1": "Omicron",
       "KP.3": "Omicron"}


def read(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def main():
    mutations = read("mutations.tsv")
    clock = read("clock.tsv")
    per_gene = read("per_gene.tsv")
    analyzed = {record.id
                for record in SeqIO.parse("clean.fasta", "fasta")}

    substitutions = Counter()
    deletions = Counter()
    for row in mutations:
        if row["sample"] == "-":
            deletions[row["accession"]] += 1
        else:
            substitutions[row["accession"]] += 1

    by_lineage = defaultdict(list)
    for row in clock:
        by_lineage[row["lineage"]].append(row["accession"])

    shared_spike = {}
    for lineage, accessions in by_lineage.items():
        counts = Counter(
            row["protein"] for row in mutations
            if row["gene"] == "S" and row["protein"]
            and row["protein"] != "del"
            and row["accession"] in set(accessions))
        threshold = max(1, int(0.8 * len(accessions)))
        shared = [change for change, carriers in counts.items()
                  if carriers >= threshold]
        shared_spike[lineage] = sorted(
            shared, key=lambda change: int(change[1:-1]))

    lines = ["# Tracking a virus as it evolves", "",
             f"{len(analyzed)} SARS-CoV-2 genomes from 13 Pango "
             "lineages,",
             "compared with the Wuhan-Hu-1 reference genome NC_045512.2.",
             "", "## Changes per lineage", "",
             "| Lineage | WHO name | Genomes | Substitutions "
             "| Deleted bases |",
             "|---|---|---|---|---|"]
    for lineage, accessions in by_lineage.items():
        subs = [substitutions[a] for a in accessions]
        dels = [deletions[a] for a in accessions]
        lines.append(
            f"| {lineage} | {WHO.get(lineage, '')} | {len(accessions)} | "
            f"{statistics.median(subs):.0f} | "
            f"{statistics.median(dels):.0f} |")

    lines += ["", "## Where the changes fall", "",
              "| Gene | Bases | Positions changed | Per thousand bases |",
              "|---|---|---|---|"]
    for row in per_gene[:6]:
        lines.append(f"| {row['gene']} | {row['length']} | "
                     f"{row['positions']} | {row['per_kb']} |")

    lines += ["", "## The spike changes each lineage shares", ""]
    for lineage, changes in shared_spike.items():
        if changes:
            lines.append(f"{lineage}: {', '.join(changes)}")
            lines.append("")

    with open("report.md", "w") as out:
        out.write("\n".join(lines) + "\n")
    print(f"report.md written, {len(lines)} lines")


if __name__ == "__main__":
    main()
