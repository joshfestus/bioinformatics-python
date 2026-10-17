"""Keep the hits strong enough to call a gene, and sort them by kind.

A BLAST hit is not a finding, for two reasons. Every bacterium carries
proteins that resemble resistance proteins without conferring resistance,
so a threshold has to be chosen and defended. And the catalog holds two
kinds of entry: genes that confer resistance by being present, and
reference copies of ordinary genes, such as gyrA, that confer resistance
only when they carry a particular mutation. A hit to the second kind says
the genome has a gyrA, which every one of them does.

uv run filter.py            reads hits.tsv, writes resistance.tsv
"""

import argparse
import csv
from collections import Counter

from Bio import SeqIO

FIELDS = ["query", "reference", "identity", "length", "qlen", "slen",
          "qstart", "qend", "sstart", "send", "evalue", "bitscore"]


def read_hits():
    with open("hits.tsv") as f:
        for row in csv.DictReader(f, FIELDS, delimiter="\t"):
            row["identity"] = float(row["identity"])
            row["bitscore"] = float(row["bitscore"])
            # how much of the reference protein the alignment covers
            row["coverage"] = 100 * (int(row["send"]) - int(row["sstart"])
                                     + 1) / int(row["slen"])
            yield row


def catalog():
    """{protein accession: the catalog's row for it}

    An entry is keyed by whichever accessions it has. Most carry a RefSeq
    one, but some are known only by their original GenBank accession, and
    looking up just the first kind leaves those hits unnamed.
    """
    named = {}
    with open("amr_catalog.tsv") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            for column in ("refseq_protein_accession",
                           "genbank_protein_accession"):
                if row.get(column):
                    named.setdefault(row[column], row)
    return named


def where(path="proteins.fasta"):
    """{locus tag: (replicon, start, strand, product)}"""
    place = {}
    for record in SeqIO.parse(path, "fasta"):
        position, product = record.description.split(" ", 2)[1:]
        replicon, start, strand = position.split(":")
        place[record.id] = (replicon, int(start), strand, product)
    return place


def best_hits(min_identity, min_coverage):
    best = {}
    for hit in read_hits():
        if (hit["identity"] < min_identity
                or hit["coverage"] < min_coverage):
            continue
        keeping = best.get(hit["query"])
        if keeping is None or hit["bitscore"] > keeping["bitscore"]:
            best[hit["query"]] = hit
    return best


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--min-identity", type=float, default=90.0)
    parser.add_argument("--min-coverage", type=float, default=90.0)
    args = parser.parse_args()

    named, place = catalog(), where()
    rows = []
    for tag, hit in best_hits(args.min_identity,
                              args.min_coverage).items():
        entry = named.get(hit["reference"].split("|")[0], {})
        replicon, start, strand, annotated = place[tag]
        rows.append({
            "gene": entry.get("allele") or entry.get("gene_family", ""),
            "subtype": entry.get("subtype", ""),
            "locus_tag": tag,
            "replicon": replicon,
            "start": start,
            "strand": strand,
            "identity": round(hit["identity"], 1),
            "coverage": round(hit["coverage"], 1),
            "class": entry.get("class", ""),
            "product": entry.get("product_name") or annotated,
        })
    rows.sort(key=lambda row: (row["subtype"], row["class"], row["gene"]))
    with open("resistance.tsv", "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=rows[0], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    acquired = [row for row in rows if row["subtype"] == "AMR"]
    print(f"{len(rows)} proteins pass {args.min_identity:.0f}% identity "
          f"over {args.min_coverage:.0f}% of the reference")
    for subtype, count in Counter(row["subtype"] for row in rows).items():
        print(f"  {count:>3} {subtype}")

    print("\nresistance genes present, mobile and intrinsic:\n")
    print(f"{'gene':14}{'class':18}{'replicon':14}{'id%':>6}{'cov%':>6}")
    for row in acquired:
        print(f"{row['gene']:14}{row['class'][:17]:18}"
              f"{row['replicon']:14}{row['identity']:>6}"
              f"{row['coverage']:>6}")


if __name__ == "__main__":
    main()
