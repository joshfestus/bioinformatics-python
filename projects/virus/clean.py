"""Keep only genomes that are long enough and have few unknown bases.

The thresholds are arguments, with the defaults the book used, so that the
pipeline can pass the values in config.yaml and you can try others without
editing this file.

uv run clean.py                    reads genomes.fasta, writes clean.fasta
uv run clean.py --max-unknown 0.05  a looser threshold
"""

import argparse
import csv

from Bio import SeqIO


def unknown_share(sequence):
    """The share of positions that are not a plain A, C, G or T.

    N is the commonest, but a consensus may also hold an IUPAC letter
    such as Y for "C or T", and those are not called bases either.
    """
    called = sum(sequence.count(base) for base in "ACGT")
    return 1 - called / len(sequence)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--min-length", type=int, default=29_000,
                        help="shortest genome to keep (the reference is "
                             "29,903 bases)")
    parser.add_argument("--max-unknown", type=float, default=0.01,
                        help="largest share of uncalled bases to allow")
    args = parser.parse_args()

    with open("samples.tsv") as f:
        lineage = {row["accession"]: row["lineage"]
                   for row in csv.DictReader(f, delimiter="\t")}

    kept, dropped = [], []
    for record in SeqIO.parse("genomes.fasta", "fasta"):
        sequence = str(record.seq).upper()
        unknown = unknown_share(sequence)
        if (len(sequence) >= args.min_length
                and unknown <= args.max_unknown):
            record.description = ""
            kept.append(record)
        else:
            dropped.append((record.id, lineage[record.id], len(sequence),
                            unknown))

    SeqIO.write(kept, "clean.fasta", "fasta")
    print(f"kept {len(kept)} of {len(kept) + len(dropped)} genomes, "
          f"at most {args.max_unknown:.0%} unknown")
    for accession, name, length, unknown in dropped:
        print(f"  dropped {accession:11} {name:10} {length:6,} bases, "
              f"{unknown:.1%} unknown")


if __name__ == "__main__":
    main()
