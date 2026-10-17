"""Summarize the assembly: what the isolate carries, and how big.

uv run describe.py          reads genome.gb, writes replicons.tsv
"""

import csv

from Bio import SeqIO
from Bio.SeqUtils import gc_fraction


def replicons():
    for record in SeqIO.parse("genome.gb", "genbank"):
        plasmid = ""
        for feature in record.features:
            if feature.type == "source":
                plasmid = feature.qualifiers.get("plasmid", [""])[0]
        genes = [f for f in record.features if f.type == "CDS"]
        pseudo = [f for f in genes if "pseudo" in f.qualifiers]
        yield {
            "accession": record.id,
            "kind": f"plasmid {plasmid}" if plasmid else "chromosome",
            "length": len(record),
            "gc": round(100 * gc_fraction(record.seq), 1),
            "proteins": len(genes) - len(pseudo),
            "pseudogenes": len(pseudo),
        }


def main():
    rows = list(replicons())
    with open("replicons.tsv", "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=rows[0], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    print(f"{'accession':14}{'':2}{'kind':18}{'bases':>9} {'GC%':>5} "
          f"{'proteins':>9}")
    for row in rows:
        print(f"{row['accession']:14}{'':2}{row['kind']:18}"
              f"{row['length']:>9,} {row['gc']:>5} {row['proteins']:>9,}")
    total = sum(row["length"] for row in rows)
    plasmid = sum(row["length"] for row in rows
                  if "plasmid" in row["kind"])
    print(f"\n{total:,} bases in {len(rows)} pieces, "
          f"{plasmid:,} of them ({plasmid / total:.1%}) on plasmids")
    print(f"{sum(row['proteins'] for row in rows):,} proteins, "
          f"{sum(row['pseudogenes'] for row in rows)} pseudogenes")


if __name__ == "__main__":
    main()
