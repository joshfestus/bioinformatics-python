# One row per human gene, from NCBI's annotation of GRCh38.p14.
import pandas as pd

features = pd.read_csv(
    "../data/human_feature_table.txt.gz",
    sep="\t",
    dtype={"chromosome": str},
)
features = features.rename(columns={"# feature": "feature"})

# Keep the 24 assembled chromosomes; leave out alternative versions of
# regions and pieces not yet placed on a chromosome.
placed = (features["assembly_unit"] == "Primary Assembly") & (
    features["seq_type"] == "chromosome"
)
features = features[placed]

genes = features[features["feature"] == "gene"]
genes = genes[["symbol", "GeneID", "chromosome", "start", "end", "class"]]
genes = genes.rename(columns={"class": "kind"})
genes["length"] = genes["end"] - genes["start"] + 1

# A gene can make several proteins; keep the length of the longest.
cds = features[features["feature"] == "CDS"]
longest = cds.groupby("GeneID")["product_length"].max()
longest = longest.astype("Int64").rename("protein_length")
genes = genes.merge(longest, on="GeneID", how="left")

genes.to_csv("genes.tsv", sep="\t", index=False)
print(f"{len(genes):,} genes written to genes.tsv")
