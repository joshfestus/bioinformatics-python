"""Which genes dexamethasone changes, with each donor as its own control.

Chapter 19's program, reading this folder's files and writing the full
table of results for the steps that follow.

uv run de.py                writes results.tsv
"""

import gzip

import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats

KEEP = ["Untreated", "Dexamethasone"]


def symbols():
    """{GeneID: symbol} from NCBI's gene_info."""
    names = {}
    with gzip.open("gene_info.tsv.gz", "rt") as handle:
        header = handle.readline().lstrip("#").rstrip("\n").split("\t")
        for line in handle:
            row = dict(zip(header, line.rstrip("\n").split("\t")))
            names[int(row["GeneID"])] = row["Symbol"]
    return names


def main():
    counts = pd.read_csv("counts.tsv.gz", sep="\t", index_col="GeneID")
    samples = pd.read_csv("samples.tsv", sep="\t", index_col="sample")
    samples = samples[samples["treatment"].isin(KEEP)]
    counts = counts[samples.index].T
    counts = counts.loc[:, counts.sum() >= 10]

    inference = DefaultInference(n_cpus=2)
    dataset = DeseqDataSet(counts=counts, metadata=samples,
                           design="~cell_line + treatment",
                           inference=inference, quiet=True)
    dataset.deseq2()
    stats = DeseqStats(dataset,
                       contrast=["treatment", "Dexamethasone",
                                 "Untreated"],
                       inference=inference, quiet=True)
    stats.summary()
    results = stats.results_df
    results.index = results.index.astype(int)
    results.insert(0, "symbol",
                   pd.Series(symbols()).reindex(results.index))
    results.to_csv("results.tsv", sep="\t")

    significant = results["padj"] < 0.05
    changed = significant & (results["log2FoldChange"].abs() >= 1)
    print(f"{len(results):,} genes tested, "
          f"{significant.sum():,} with adjusted p below 0.05")
    print(f"{changed.sum():,} of those changed at least twofold")


if __name__ == "__main__":
    main()
