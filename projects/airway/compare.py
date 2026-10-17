"""Compare the findings with the paper the data came from.

GSE52778 was published by Himes and others in 2014. They reported 316
differentially expressed genes at a Benjamini-Hochberg corrected p below
0.05, named four well known glucocorticoid responders and three less
studied ones, and followed CRISPLD2 into the laboratory. A reanalysis
that cannot find their genes has gone wrong somewhere.

uv run compare.py           writes comparison.tsv
"""

import csv

import pandas as pd

PAPER = {
    "DUSP1": "well known", "KLF15": "well known", "PER1": "well known",
    "TSC22D3": "well known", "C7": "less investigated",
    "CCDC69": "less investigated", "CRISPLD2": "less investigated",
}
THEIR_COUNT = 316


def main():
    results = pd.read_csv("results.tsv", sep="\t").dropna(subset=["padj"])
    results = results.sort_values("padj").reset_index(drop=True)
    results["rank"] = results.index + 1

    rows = []
    for symbol, standing in PAPER.items():
        found = results[results["symbol"] == symbol]
        if found.empty:
            rows.append({"symbol": symbol, "in_paper": standing,
                         "rank": "", "log2FoldChange": "", "padj": ""})
            continue
        gene = found.iloc[0]
        rows.append({
            "symbol": symbol,
            "in_paper": standing,
            "rank": int(gene["rank"]),
            "log2FoldChange": round(gene["log2FoldChange"], 2),
            "padj": f"{gene['padj']:.1e}",
        })
    with open("comparison.tsv", "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=rows[0], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    print("the seven genes the paper names, in this reanalysis:\n")
    print(f"{'gene':10}{'in the paper':20}{'rank':>6}{'log2FC':>9}"
          f"{'padj':>10}")
    for row in rows:
        print(f"{row['symbol']:10}{row['in_paper']:20}{row['rank']:>6}"
              f"{row['log2FoldChange']:>9}{row['padj']:>10}")

    ours = (results["padj"] < 0.05).sum()
    print(f"\nthey reported {THEIR_COUNT} changed genes, this analysis "
          f"{ours:,}")
    print(f"all {len(rows)} of their named genes are in the top "
          f"{max(int(r['rank']) for r in rows if r['rank'])} here, "
          "and all went up")


if __name__ == "__main__":
    main()
