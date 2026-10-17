"""List every difference between each genome and the reference.

Reads aligned.fasta, where MAFFT kept the reference's numbering, so a
column of the alignment is a position in NC_045512.2. Writes one row per
mutation per genome, with its amino acid change where it falls in a gene.

uv run mutations.py        writes mutations.tsv
"""

import csv
from collections import Counter

from Bio import SeqIO
from Bio.Seq import Seq

# A consensus sequence is not limited to N for an uncertain position. The
# IUPAC alphabet has a letter for each pair of bases, such as Y for "C or
# T", and a genome that uses them will quietly produce amino acid changes
# that were never called. Anything that is not one plain base is unknown.
BASES = "ACGT"


def coding_positions(record):
    """{gene: [reference position of every base of its CDS, in order]}."""
    genes = {}
    for feature in record.features:
        if feature.type != "CDS":
            continue
        name = feature.qualifiers["gene"][0]
        product = feature.qualifiers["product"][0]
        if name == "ORF1ab" and "ORF1a " in product:
            # ORF1a is the shorter form of the same reading frame
            continue
        positions = []
        for part in feature.location.parts:
            positions += list(range(int(part.start), int(part.end)))
        genes[name] = positions
    return genes


def peptides(record):
    """[(start, end, name)] for the mature proteins cut out of ORF1ab."""
    found, seen = [], set()
    for feature in record.features:
        if feature.type != "mat_peptide":
            continue
        name = feature.qualifiers["product"][0]
        if name in seen:
            continue
        seen.add(name)
        found.append((int(feature.location.start),
                      int(feature.location.end), name))
    return found


def covered(sequence):
    """The first and last position the sequencing actually reached.

    Almost no SARS-CoV-2 genome covers the extreme ends of the molecule,
    and MAFFT writes the missing ends as gaps. Counting those gaps as
    deletions would add hundreds of false mutations to every genome.
    """
    first = len(sequence) - len(sequence.lstrip("-"))
    last = len(sequence.rstrip("-")) - 1
    return first, last


def protein_changes(reference, sample, positions, span):
    """Compare one gene codon by codon. Yields (codon number, change)."""
    for index in range(0, len(positions) - 2, 3):
        codon = positions[index:index + 3]
        if not all(span[0] <= p <= span[1] for p in codon):
            continue
        theirs = "".join(sample[p] for p in codon)
        if any(base not in BASES for base in theirs):
            if theirs == "---":
                yield index // 3 + 1, "del"
            continue
        ours = "".join(reference[p] for p in codon)
        if ours == theirs:
            continue
        before = Seq(ours).translate()
        after = Seq(theirs).translate()
        if before != after:
            yield index // 3 + 1, f"{before}{index // 3 + 1}{after}"


def main():
    annotation = SeqIO.read("reference.gb", "genbank")
    genes = coding_positions(annotation)
    mature = peptides(annotation)
    gene_of = {}
    for name, positions in genes.items():
        for p in positions:
            gene_of.setdefault(p, name)

    records = {r.id: str(r.seq).upper()
               for r in SeqIO.parse("aligned.fasta", "fasta")}
    reference = records.pop(annotation.id)
    with open("samples.tsv") as f:
        samples = {row["accession"]: row
                   for row in csv.DictReader(f, delimiter="\t")}

    rows, counts = [], Counter()
    for accession, sequence in records.items():
        sample = samples[accession]
        span = covered(sequence)
        changes = {}
        for gene, positions in genes.items():
            for codon, change in protein_changes(reference, sequence,
                                                 positions, span):
                changes[gene, codon] = change
        pairs = enumerate(zip(reference, sequence))
        for position, (ours, theirs) in pairs:
            if not span[0] <= position <= span[1]:
                continue
            if ours == theirs or (theirs not in BASES and theirs != "-"):
                continue
            gene = gene_of.get(position, "")
            codon = 0
            if gene:
                codon = genes[gene].index(position) // 3 + 1
            rows.append({
                "accession": accession,
                "lineage": sample["lineage"],
                "date": sample["date"],
                "position": position + 1,  # one based, as papers write it
                "reference": ours,
                "sample": theirs,
                "gene": gene,
                "protein": changes.get((gene, codon), ""),
                "nsp": name_of_peptide(mature, position),
            })
            counts[accession] += 1

    rows.sort(key=lambda r: (r["lineage"], r["accession"], r["position"]))
    with open("mutations.tsv", "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=rows[0], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"{len(rows):,} differences in {len(records)} genomes")
    fewest = min(counts.values())
    most = max(counts.values())
    print(f"from {fewest} in one genome to {most} in another")


def name_of_peptide(mature, position):
    for start, end, name in mature:
        if start <= position < end:
            return name
    return ""


if __name__ == "__main__":
    main()
