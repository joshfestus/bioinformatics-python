"""Where in the genome the mutations gather.

Counts mutations per gene, and per thousand bases of gene, so that a long
gene is not called a hotspot merely for being long. Then counts how many
of the sampled genomes carry each spike mutation.

uv run locate.py           reads mutations.tsv, writes per_gene.tsv
"""

import csv
from collections import Counter, defaultdict

from Bio import SeqIO


def gene_lengths(record):
    lengths = {}
    for feature in record.features:
        if feature.type != "gene":
            continue
        name = feature.qualifiers["gene"][0]
        lengths[name] = len(feature.location)
    return lengths


def main():
    reference = SeqIO.read("reference.gb", "genbank")
    lengths = gene_lengths(reference)
    with open("mutations.tsv") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    analyzed = sum(1 for _ in SeqIO.parse("clean.fasta", "fasta"))
    # two genomes differ from the reference nowhere at all, so they have
    # no rows in mutations.tsv; Chapter 23 asks why
    genomes = len({row["accession"] for row in rows})

    # one mutation per genome counts once; a position shared by fifty
    # genomes is one place the virus changed, not fifty
    places = defaultdict(set)
    for row in rows:
        places[row["gene"] or "between genes"].add(int(row["position"]))
    events = Counter(row["gene"] or "between genes" for row in rows)

    table = []
    for gene, distinct in places.items():
        length = lengths.get(gene, 0)
        table.append({
            "gene": gene,
            "length": length,
            "positions": len(distinct),
            "per_kb": (round(1000 * len(distinct) / length, 2)
                       if length else ""),
            "changes": events[gene],
        })
    table.sort(key=lambda row: row["per_kb"] or 0, reverse=True)
    with open("per_gene.tsv", "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=table[0], delimiter="\t")
        writer.writeheader()
        writer.writerows(table)

    print(f"{len(rows):,} changes in {genomes} of the "
          f"{analyzed} genomes\n")
    print(f"{'gene':14} {'bases':>6} {'places':>7} {'per kb':>7}")
    for row in table:
        print(f"{row['gene']:14} {row['length'] or '':>6} "
              f"{row['positions']:>7} {row['per_kb'] or '':>7}")

    spike = Counter(row["protein"] for row in rows
                    if row["gene"] == "S" and row["protein"]
                    and row["protein"] != "del")
    print(f"\nthe spike changes most of the {analyzed} genomes share:")
    for change, carriers in spike.most_common(12):
        print(f"  {change:8} {carriers:3} genomes  "
              f"{carriers / analyzed:.0%}")


if __name__ == "__main__":
    main()
