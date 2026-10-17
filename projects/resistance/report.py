"""Write the resistance profile as a report, from the files produced.

uv run report.py            writes report.md
"""

import csv
from collections import defaultdict


def read(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def main():
    replicons = read("replicons.tsv")
    profile = read("profile.tsv")
    mutations = read("mutations.tsv")
    version = open("amr_version.txt").read().strip()

    by_class = defaultdict(set)
    for row in profile:
        by_class[row["class"] or "other"].add(row["gene"])
    for row in mutations:
        by_class[row["class"]].add(row["mutation"])

    plasmid_genes = [row for row in profile
                     if "plasmid" in row["replicon"]]
    lines = [
        "# Resistance profile: Klebsiella pneumoniae HS11286", "",
        f"{len(replicons)} replicons, "
        f"{sum(int(r['length']) for r in replicons):,} bases, "
        f"{sum(int(r['proteins']) for r in replicons):,} proteins.",
        f"Searched against NCBI's Reference Gene Catalog, version "
        f"{version}.", "",
        "## What the isolate is resistant to", "",
        "| Class | Genes and mutations |", "|---|---|",
    ]
    for group in sorted(by_class):
        found = ", ".join(sorted(by_class[group]))
        lines.append(f"| {group} | {found} |")

    lines += ["", "## Where the genes sit", "",
              "| Replicon | Bases | GC% | Resistance genes |",
              "|---|---|---|---|"]
    counts = defaultdict(int)
    for row in profile:
        counts[row["accession"]] += 1
    for row in replicons:
        lines.append(f"| {row['kind']} ({row['accession']}) | "
                     f"{int(row['length']):,} | {row['gc']} | "
                     f"{counts[row['accession']]} |")

    lines += ["", "## What can spread", "",
              f"{len(plasmid_genes)} of the {len(profile)} acquired "
              "genes sit on a plasmid, and so can move to another "
              "bacterium:", ""]
    for row in plasmid_genes:
        note = (", beside " + row["mobile_neighbors"].split(";")[0]
                if row["mobile_neighbors"] else "")
        lines.append(f"- {row['gene']} on {row['replicon']}{note}")
    lines += ["",
              "The genes on the chromosome are the ones this species "
              "carries anyway, whatever the isolate has picked up.", ""]

    with open("report.md", "w") as out:
        out.write("\n".join(lines) + "\n")
    print(f"report.md written, {len(lines)} lines")


if __name__ == "__main__":
    main()
