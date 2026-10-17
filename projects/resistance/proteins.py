"""Write every predicted protein in the genome as one FASTA file.

A complete RefSeq genome arrives annotated: NCBI has already found the
genes and translated them, and each one carries a locus tag, a product
name and a position. That is the annotation step of this project, and
taking it from the record rather than redoing it is the honest choice
when the record is a finished one.

uv run proteins.py          reads genome.gb, writes proteins.fasta
"""

from Bio import SeqIO

written, skipped = 0, 0
with open("proteins.fasta", "w") as out:
    for record in SeqIO.parse("genome.gb", "genbank"):
        for feature in record.features:
            if feature.type != "CDS":
                continue
            translation = feature.qualifiers.get("translation")
            if not translation:
                skipped += 1  # a pseudogene, with no protein to search
                continue
            tag = feature.qualifiers["locus_tag"][0]
            product = feature.qualifiers.get("product", ["unknown"])[0]
            start = int(feature.location.start) + 1
            strand = "+" if feature.location.strand == 1 else "-"
            out.write(f">{tag} {record.id}:{start}:{strand} {product}\n")
            out.write(translation[0] + "\n")
            written += 1

print(f"{written:,} proteins written, {skipped} CDS features skipped")
