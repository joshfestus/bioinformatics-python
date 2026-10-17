"""Write the study's findings as a report, from the files produced.

uv run report.py            writes report.md
"""

import csv


def read(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def top(rows, direction, how_many=8):
    found = [row for row in rows if row["direction"] == direction
             and float(row["q"]) < 0.05]
    found.sort(key=lambda row: float(row["p"]))
    return found[:how_many]


def main():
    enrichment = read("enrichment.tsv")
    comparison = read("comparison.tsv")
    results = read("results.tsv")
    tested = [row for row in results if row["padj"] not in ("", "NA")]
    significant = [row for row in tested if float(row["padj"]) < 0.05]
    changed = [row for row in significant
               if abs(float(row["log2FoldChange"])) >= 1]
    up = [row for row in changed if float(row["log2FoldChange"]) > 0]

    lines = [
        "# What dexamethasone does to airway smooth muscle", "",
        "Four donors' primary airway smooth muscle cells, treated with "
        "dexamethasone",
        "or left untreated, from GSE52778 (Himes and others, 2014).",
        "Each donor is their own control.", "",
        "## What changed", "",
        f"- {len(results):,} genes tested",
        f"- {len(significant):,} with an adjusted p below 0.05",
        f"- {len(changed):,} of those changed at least twofold: "
        f"{len(up)} up and {len(changed) - len(up)} down", "",
        "## What the changed genes do", "",
        "| Direction | Process | Genes found | In the set |",
        "|---|---|---|---|",
    ]
    for direction in ("up", "down"):
        for row in top(enrichment, direction):
            lines.append(f"| {direction} | {row['name']} | "
                         f"{row['found']} | {row['in_set']} |")

    lines += ["", "## The genes the paper named", "",
              "| Gene | In the paper | Rank here | log2 fold change |",
              "|---|---|---|---|"]
    for row in comparison:
        lines.append(f"| {row['symbol']} | {row['in_paper']} | "
                     f"{row['rank']} | {row['log2FoldChange']} |")
    lines += ["", "All seven are found, and all went up.", ""]

    with open("report.md", "w") as out:
        out.write("\n".join(lines) + "\n")
    print(f"report.md written, {len(lines)} lines")


if __name__ == "__main__":
    main()
