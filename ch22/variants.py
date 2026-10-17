"""Annotate confident variants with the gene and amino acid they change.

Chapter 16's effects.py rebuilt as a command line tool: its inputs and
thresholds are arguments, its progress is logged, its result is a file,
and its logic is in functions that tests can call.

    uv run variants.py calls.vcf.gz genes.gff3 genome.fasta -o out.tsv
"""

import argparse
import csv
import logging

import pysam
from Bio import SeqIO
from Bio.Seq import Seq

log = logging.getLogger("variants")


def read_coding_sequences(gff_path: str) -> list[dict]:
    """CDS features from a GFF3 file, one based closed coordinates."""
    features = []
    with open(gff_path) as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if fields[2] != "CDS":
                continue
            attributes = dict(
                pair.split("=", 1) for pair in fields[8].split(";")
            )
            features.append(
                {
                    "start": int(fields[3]),
                    "end": int(fields[4]),
                    "strand": fields[6],
                    "name": attributes.get("gene")
                    or attributes.get("product", "?")[:24],
                    "frameshift": "frameshift" in fields[8],
                }
            )
    return features


def amino_acid_change(
    gene: Seq, strand: str, start: int, end: int, position: int, alt: str
) -> str:
    """The change a one letter variant makes, such as D127N missense.

    gene is the sequence from start to end on the plus strand; start, end
    and position are one based, as in GFF and VCF.
    """
    offset = position - start
    if strand == "-":
        gene = gene.reverse_complement()
        offset = end - position
        alt = str(Seq(alt).reverse_complement())
    mutated = gene[:offset] + alt + gene[offset + 1 :]
    codon = offset // 3
    before = gene[3 * codon : 3 * codon + 3].translate(table=11)
    after = mutated[3 * codon : 3 * codon + 3].translate(table=11)
    if before == after:
        kind = "synonymous"
    elif after == "*":
        kind = "nonsense"
    else:
        kind = "missense"
    return f"{before}{codon + 1}{after} {kind}"


def annotate(vcf_path, gff_path, fasta_path, min_qual, min_depth):
    """Yield one row per variant that passes the thresholds."""
    genome = SeqIO.read(fasta_path, "fasta").seq
    coding = read_coding_sequences(gff_path)
    log.info("read %d coding sequences", len(coding))
    for record in pysam.VariantFile(vcf_path):
        if record.qual < min_qual or record.info["DP"] < min_depth:
            continue
        alt = record.alts[0]
        hits = [c for c in coding if c["start"] <= record.pos <= c["end"]]
        effect = "between genes"
        if hits:
            cds = hits[0]
            if cds["frameshift"]:
                effect = "frameshifted gene, not translated"
            elif len(record.ref) == 1 and len(alt) == 1:
                gene = genome[cds["start"] - 1 : cds["end"]]
                effect = amino_acid_change(
                    gene,
                    cds["strand"],
                    cds["start"],
                    cds["end"],
                    record.pos,
                    alt,
                )
            else:
                effect = "insertion or deletion"
            effect = f"{cds['name']}: {effect}"
        yield {
            "position": record.pos,
            "ref": record.ref,
            "alt": alt,
            "quality": round(record.qual),
            "effect": effect,
        }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("vcf", help="variant calls, VCF or BCF")
    parser.add_argument("gff", help="gene annotation, GFF3")
    parser.add_argument("fasta", help="reference genome, FASTA")
    parser.add_argument("-o", "--output", required=True, help="TSV file")
    parser.add_argument("--min-qual", type=float, default=30)
    parser.add_argument("--min-depth", type=int, default=10)
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
    )
    log.info(
        "thresholds: QUAL >= %s, depth >= %s",
        args.min_qual,
        args.min_depth,
    )
    rows = list(
        annotate(
            args.vcf, args.gff, args.fasta, args.min_qual, args.min_depth
        )
    )
    with open(args.output, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(rows)
    log.info("wrote %d variants to %s", len(rows), args.output)


if __name__ == "__main__":
    main()
