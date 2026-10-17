"""Shared loading for Chapter 20: log2 counts per million of all sixteen
GSE52778 samples, with short treatment labels."""

import numpy as np
import pandas as pd

SHORT = {
    "Untreated": "none",
    "Dexamethasone": "dex",
    "Albuterol": "alb",
    "Albuterol_Dexamethasone": "alb+dex",
}


def load() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Samples as rows and expressed genes as columns, plus samples."""
    counts = pd.read_csv(
        "../data/GSE52778_counts.tsv.gz", sep="\t", index_col="GeneID"
    )
    samples = pd.read_csv("samples.tsv", sep="\t", index_col="sample")
    samples["label"] = samples["treatment"].map(SHORT)
    cpm = counts[samples.index] / counts[samples.index].sum() * 1e6
    expressed = cpm[cpm.mean(axis=1) >= 1]
    return np.log2(expressed + 1).T, samples
