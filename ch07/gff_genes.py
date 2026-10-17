# The genes of E. coli K-12 MG1655 from NCBI's GFF3 annotation.
import pandas as pd

COLUMNS = [
    "seqid",
    "source",
    "type",
    "start",
    "end",
    "score",
    "strand",
    "phase",
    "attributes",
]


def parse_attributes(text: str) -> dict[str, str]:
    """ID=gene-b0344;Name=lacZ;... as a dictionary."""
    return dict(pair.split("=", 1) for pair in text.split(";"))


gff = pd.read_csv(
    "ecoli_k12.gff3.gz", sep="\t", comment="#", names=COLUMNS
)
genes = gff[gff["type"] == "gene"].copy()
attributes = genes["attributes"].map(parse_attributes)
genes["name"] = attributes.map(lambda a: a["Name"])
genes["biotype"] = attributes.map(lambda a: a["gene_biotype"])
genes["length"] = genes["end"] - genes["start"] + 1

print(genes["biotype"].value_counts().to_string())
print()
lac = genes[genes["name"].str.startswith("lac")]
print(lac[["name", "start", "end", "strand", "length"]].to_string())
