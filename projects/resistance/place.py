"""Ask where each resistance gene sits, and whether it looks mobile.

A resistance gene on the chromosome stays in its own line of descent. One
on a plasmid can move to another bacterium, and one sitting among
transposases can move to another plasmid. The question a hospital asks is
not only what the isolate is resistant to, but what it can pass on.

uv run place.py             reads resistance.tsv, writes profile.tsv
"""

import csv
from collections import defaultdict

from Bio import SeqIO

NEARBY = 5_000  # bases either side counted as the gene's neighborhood
MOBILE = ("transposase", "integrase", "recombinase", "insertion sequence",
          "IS110", "resolvase")


def replicons():
    kinds, neighbors = {}, defaultdict(list)
    for record in SeqIO.parse("genome.gb", "genbank"):
        plasmid = ""
        for feature in record.features:
            if feature.type == "source":
                plasmid = feature.qualifiers.get("plasmid", [""])[0]
            elif feature.type == "CDS":
                neighbors[record.id].append((
                    int(feature.location.start) + 1,
                    feature.qualifiers.get("product", ["unknown"])[0]))
        kinds[record.id] = (f"plasmid {plasmid}" if plasmid
                            else "chromosome")
        neighbors[record.id].sort()
    return kinds, neighbors


def mobile_near(neighbors, start):
    return sorted({product for position, product in neighbors
                   if abs(position - start) <= NEARBY
                   and any(word in product.lower() for word in MOBILE)})


def main():
    kinds, neighbors = replicons()
    with open("resistance.tsv") as f:
        genes = [row for row in csv.DictReader(f, delimiter="\t")
                 if row["subtype"] == "AMR"]

    rows = []
    for gene in genes:
        around = mobile_near(neighbors[gene["replicon"]],
                             int(gene["start"]))
        rows.append({
            "gene": gene["gene"],
            "class": gene["class"],
            "replicon": kinds[gene["replicon"]],
            "accession": gene["replicon"],
            "start": gene["start"],
            "mobile_neighbors": "; ".join(around),
        })
    rows.sort(key=lambda row: (row["replicon"], int(row["start"])))
    with open("profile.tsv", "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=rows[0], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    where = defaultdict(list)
    for row in rows:
        where[row["replicon"]].append(row)
    for replicon, found in where.items():
        print(f"{replicon} ({len(found)} genes)")
        for row in found:
            flag = " next to a mobility gene" if row["mobile_neighbors"] \
                else ""
            print(f"  {row['gene']:14}{row['class'][:16]:18}"
                  f"{int(row['start']):>9,}{flag}")
        print()
    movable = [row for row in rows if "plasmid" in row["replicon"]]
    print(f"{len(movable)} of the {len(rows)} genes sit on a plasmid "
          "and can move to another bacterium;")
    print(f"the other {len(rows) - len(movable)} are chromosomal, ones "
          "this species carries anyway")


if __name__ == "__main__":
    main()
