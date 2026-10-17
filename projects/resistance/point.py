"""Check the point mutations, instead of assuming them from a BLAST hit.

Some resistance comes not from a gene the bacterium acquired but from one
letter of a gene every bacterium already has. The catalog holds the
reference copies of those genes, such as gyrA, so a BLAST hit to one of
them proves only that the isolate has a gyrA. The mutation table says
which position matters, and for which organism: a rule written for
Escherichia says nothing about a Klebsiella.

uv run point.py             reads hits.tsv, writes mutations.tsv
"""

import csv
from collections import defaultdict

from Bio import Align, SeqIO
from Bio.Align import substitution_matrices

ORGANISM = "Klebsiella_pneumoniae"
MIN_IDENTITY, MIN_COVERAGE = 90.0, 90.0
FIELDS = ["query", "reference", "identity", "length", "qlen", "slen",
          "qstart", "qend", "sstart", "send", "evalue", "bitscore"]


def rules():
    """{accession: [(position, symbol, class)]} for one organism."""
    by_protein = defaultdict(list)
    with open("amr_mutations.tsv") as f:
        header = f.readline().lstrip("#").rstrip("\n").split("\t")
        for row in csv.DictReader(f, header, delimiter="\t"):
            if row["taxgroup"] != ORGANISM:
                continue
            by_protein[row["accession_version"]].append((
                int(row["mutation_position"]),
                row["standard_mutation_symbol"],
                row["class"],
            ))
    return by_protein


def aligner():
    pairwise = Align.PairwiseAligner()
    pairwise.substitution_matrix = substitution_matrices.load("BLOSUM62")
    pairwise.open_gap_score = -11
    pairwise.extend_gap_score = -1
    return pairwise


def residue_at(alignment, position):
    """The isolate's residue aligned to `position` of the reference."""
    isolate, reference = alignment.aligned
    for (mine, _), (theirs, ends) in zip(isolate, reference):
        if theirs <= position - 1 < ends:
            return alignment.target[mine + (position - 1 - theirs)]
    return "-"  # the reference position falls in a gap


def candidates(ruled):
    """The isolate's proteins that match a gene the rules are about."""
    best = {}
    with open("hits.tsv") as f:
        for row in csv.DictReader(f, FIELDS, delimiter="\t"):
            accession = row["reference"].split("|")[0]
            if accession not in ruled:
                continue
            coverage = 100 * (int(row["send"]) - int(row["sstart"]) + 1) \
                / int(row["slen"])
            if float(row["identity"]) < MIN_IDENTITY \
                    or coverage < MIN_COVERAGE:
                continue
            key = (row["query"], accession)
            if key not in best or float(row["bitscore"]) > best[key]:
                best[key] = float(row["bitscore"])
    return sorted(best)


def main():
    ruled = rules()
    proteins = SeqIO.to_dict(SeqIO.parse("proteins.fasta", "fasta"))
    references = {}
    for record in SeqIO.parse("amr_proteins.fasta", "fasta"):
        references.setdefault(record.id.split("|")[0], record)

    pairwise = aligner()
    found, checked = [], 0
    for tag, accession in candidates(ruled):
        alignment = pairwise.align(proteins[tag].seq,
                                   references[accession].seq)[0]
        for position, symbol, group in sorted(set(ruled[accession])):
            change = symbol.split("_")[1]
            wild, mutant = change[0], change[-1]
            here = residue_at(alignment, position)
            checked += 1
            if here == mutant and mutant != wild:
                found.append({"mutation": symbol, "locus_tag": tag,
                              "position": position, "residue": here,
                              "class": group})
    with open("mutations.tsv", "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=["mutation", "locus_tag",
                                                 "position", "residue",
                                                 "class"], delimiter="\t")
        writer.writeheader()
        writer.writerows(found)

    print(f"{len(ruled)} reference proteins carry {ORGANISM} rules")
    print(f"{checked} positions checked in this isolate\n")
    for row in found:
        print(f"  {row['mutation']:16} {row['locus_tag']:12} "
              f"{row['class']}")
    print(f"\n{len(found)} resistance mutations found")


if __name__ == "__main__":
    main()
