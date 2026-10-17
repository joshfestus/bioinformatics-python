"""Check the eight samples before testing anything with them.

Every analysis that follows assumes these libraries are comparable and
that the treated and untreated samples differ. Both are worth confirming,
because a sample sheet can be wrong and a library can fail, and either
produces a result that looks fine.

uv run inspect_samples.py   prints the checks, writes qc.tsv
"""

import csv

import numpy as np
import pandas as pd

KEEP = ["Untreated", "Dexamethasone"]


def load():
    counts = pd.read_csv("counts.tsv.gz", sep="\t", index_col="GeneID")
    samples = pd.read_csv("samples.tsv", sep="\t", index_col="sample")
    samples = samples[samples["treatment"].isin(KEEP)]
    return counts[samples.index], samples


def main():
    counts, samples = load()
    detected = (counts > 0).sum()
    depth = counts.sum()
    samples["reads"] = depth
    samples["genes_detected"] = detected
    samples.to_csv("qc.tsv", sep="\t")

    print(f"{'sample':12}{'donor':10}{'treatment':16}"
          f"{'reads':>12}{'genes':>8}")
    for name, row in samples.iterrows():
        print(f"{name:12}{row['cell_line']:10}{row['treatment']:16}"
              f"{row['reads']:>12,}{row['genes_detected']:>8,}")
    print(f"\nlibrary sizes differ by {depth.max() / depth.min():.1f} "
          "fold, which normalization handles")

    # log counts of the genes seen everywhere, for correlations and PCA
    kept = counts[(counts > 0).all(axis=1)]
    logged = np.log2(kept + 1)
    correlation = logged.corr(method="pearson")
    treated = samples["treatment"] == "Dexamethasone"
    pairs = [(a, b) for i, a in enumerate(samples.index)
             for b in samples.index[i + 1:]]
    same = [correlation.loc[a, b] for a, b in pairs
            if treated[a] == treated[b]]
    across = [correlation.loc[a, b] for a, b in pairs
              if treated[a] != treated[b]]
    print(f"\n{len(kept):,} genes seen in all eight samples")
    print(f"correlation within a treatment {min(same):.3f} to "
          f"{max(same):.3f}")
    print(f"correlation across treatments  {min(across):.3f} to "
          f"{max(across):.3f}")

    centered = (logged - logged.mean(axis=1).values[:, None]).T
    _, variance, components = np.linalg.svd(centered, full_matrices=False)
    share = variance ** 2 / (variance ** 2).sum()
    scores = centered.to_numpy() @ components[:2].T
    print(f"\nPC1 holds {share[0]:.0%} of the variation and PC2 "
          f"{share[1]:.0%}")
    for name, (first, second) in zip(samples.index, scores):
        mark = "dex" if treated[name] else "none"
        print(f"  {name:12}{samples.loc[name, 'cell_line']:10}{mark:6}"
              f"{first:>9.1f}{second:>9.1f}")


if __name__ == "__main__":
    main()
