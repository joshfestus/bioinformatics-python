# Three ways to make RNA-seq counts comparable between samples, applied to
# the untreated and dexamethasone treated airway cells of GSE52778.
import numpy as np
import pandas as pd

counts = pd.read_csv(
    "../data/GSE52778_counts.tsv.gz", sep="\t", index_col="GeneID"
)
samples = pd.read_csv("samples.tsv", sep="\t")
samples = samples[
    samples["treatment"].isin(["Untreated", "Dexamethasone"])
]
counts = counts[samples["sample"]]
counts = counts[counts.sum(axis=1) > 0]  # drop genes no read reached
print(f"{len(counts):,} genes with any reads in {len(samples)} samples")

# Counts per million: divide by the library size
library = counts.sum()
cpm = counts / library * 1e6

# Transcripts per million: divide by transcript length first
features = pd.read_csv(
    "../data/human_feature_table.txt.gz",
    sep="\t",
    dtype={"chromosome": str},
)
rnas = features[features["# feature"].isin(["mRNA", "ncRNA", "misc_RNA"])]
length = rnas.groupby("GeneID")["product_length"].max() / 1000  # kb
with_length = counts.loc[counts.index.isin(length.index)]
per_kb = with_length.div(length[with_length.index], axis=0)
tpm = per_kb / per_kb.sum() * 1e6

# Median of ratios: compare each sample with a pseudo sample made of every
# gene's geometric mean, using only genes with reads in every sample
expressed = counts[(counts > 0).all(axis=1)]
log_mean = np.log(expressed).mean(axis=1)
ratios = np.log(expressed).sub(log_mean, axis=0)
size_factors = np.exp(ratios.median())

table = pd.DataFrame(
    {
        "treatment": samples.set_index("sample")["treatment"].map(
            {"Untreated": "none", "Dexamethasone": "dex"}
        ),
        "millions": (library / 1e6).round(1),
        "size_factor": size_factors.round(2),
        "FKBP5_cpm": cpm.loc[2289].round(0),
        "FKBP5_tpm": tpm.loc[2289].round(0),
    }
)
print(table.to_string())
cpm.to_csv("cpm.tsv.gz", sep="\t")
