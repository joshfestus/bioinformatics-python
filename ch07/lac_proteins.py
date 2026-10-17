# Write the proteins of the lac genes, from their GenBank annotation, as
# a FASTA file.
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

region = SeqIO.read("../ch06/lac_region.gb", "genbank")

proteins = []
for feature in region.features:
    if feature.type != "CDS":
        continue
    gene = feature.qualifiers["gene"][0]
    if not gene.startswith("lac"):
        continue
    protein = SeqRecord(
        Seq(feature.qualifiers["translation"][0]),
        id=feature.qualifiers["protein_id"][0],
        description=f"{gene} {feature.qualifiers['product'][0]}",
    )
    proteins.append(protein)

count = SeqIO.write(proteins, "lac_proteins.fasta", "fasta")
print(f"wrote {count} proteins to lac_proteins.fasta")
