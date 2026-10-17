# Which genes does dexamethasone change in airway smooth muscle cells?
# Run with "paired" (the default) to let each donor's cells be their own
# control, or "unpaired" to ignore which donor each sample came from.
import sys

import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats

counts = pd.read_csv(
    "../data/GSE52778_counts.tsv.gz", sep="\t", index_col="GeneID"
)
samples = pd.read_csv("samples.tsv", sep="\t", index_col="sample")
samples = samples[
    samples["treatment"].isin(["Untreated", "Dexamethasone"])
]

# PyDESeq2 wants samples as rows and genes as columns
counts = counts[samples.index].T
counts = counts.loc[
    :, counts.sum() >= 10
]  # drop genes with almost no reads

# "paired" lets each donor be its own control; "unpaired" ignores donors
paired = sys.argv[1:] != ["unpaired"]
design = "~cell_line + treatment" if paired else "~treatment"
print(f"design {design}")

inference = DefaultInference(n_cpus=2)
dataset = DeseqDataSet(
    counts=counts,
    metadata=samples,
    design=design,
    inference=inference,
    quiet=True,
)
dataset.deseq2()
stats = DeseqStats(
    dataset,
    contrast=["treatment", "Dexamethasone", "Untreated"],
    inference=inference,
    quiet=True,
)
stats.summary()
results = stats.results_df

genes = pd.read_csv(
    "../ch05/genes.tsv", sep="\t", dtype={"chromosome": str}
)
genes = genes.drop_duplicates("GeneID").set_index("GeneID")["symbol"]
results.index = results.index.astype(int)
results.insert(0, "symbol", genes.reindex(results.index))
if paired:
    results.to_csv("results.tsv", sep="\t")

significant = results["padj"] < 0.05
changed = significant & (results["log2FoldChange"].abs() >= 1)
up = changed & (results["log2FoldChange"] > 0)
print(f"{len(results):,} genes tested")
print(f"{significant.sum():,} with adjusted p below 0.05")
print(f"{changed.sum():,} of them changed at least twofold:")
print(f"  {up.sum()} up and {(changed & ~up).sum()} down")
top = results[significant].sort_values("padj").head(10)
table = top[["symbol", "baseMean", "log2FoldChange", "padj"]]
print(
    table.to_string(
        formatters={
            "baseMean": "{:.0f}".format,
            "log2FoldChange": "{:.2f}".format,
            "padj": "{:.1e}".format,
        }
    )
)
