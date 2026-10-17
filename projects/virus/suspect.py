"""Find genomes whose mutations do not fit their collection date.

Quality control is not one step at the start of a project. A genome can
pass every check in clean.py, which only counts letters, and still be
wrong in a way that only the finished analysis reveals.

uv run suspect.py          reads clock.tsv, prints the genomes to look at
"""

import csv

from timeline import fit


def main():
    with open("clock.tsv") as f:
        rows = [{**row,
                 "days": int(row["days"]),
                 "substitutions": int(row["substitutions"])}
                for row in csv.DictReader(f, delimiter="\t")]

    rate, intercept = fit([row["days"] for row in rows],
                          [row["substitutions"] for row in rows])
    for row in rows:
        expected = rate * row["days"] + intercept
        row["expected"] = expected
        row["residual"] = row["substitutions"] - expected

    identical = [row for row in rows if row["substitutions"] == 0]
    print(f"{len(identical)} genome(s) identical to the reference:")
    for row in identical:
        print(f"  {row['accession']:12} {row['lineage']:9} "
              f"{row['date']}  "
              f"{row['days']} days after the reference")

    rows.sort(key=lambda row: row["residual"])
    print("\nfurthest below the number of changes their date predicts:")
    for row in rows[:5]:
        print(f"  {row['accession']:12} {row['lineage']:9} {row['date']} "
              f"{row['substitutions']:4} changes, "
              f"{row['expected']:5.0f} predicted")
    print("\nfurthest above:")
    for row in rows[-3:]:
        print(f"  {row['accession']:12} {row['lineage']:9} {row['date']} "
              f"{row['substitutions']:4} changes, "
              f"{row['expected']:5.0f} predicted")


if __name__ == "__main__":
    main()
