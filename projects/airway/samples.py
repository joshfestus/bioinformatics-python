"""Turn GEO's series matrix into a table of samples.

This is Chapter 18's program with its paths changed, so that the project
folder stands on its own.

uv run samples.py           writes samples.tsv
"""

import gzip

import pandas as pd

fields = {}
with gzip.open("series_matrix.txt.gz", "rt") as handle:
    for line in handle:
        if not line.startswith("!Sample_"):
            continue
        key, *values = line.rstrip("\n").split("\t")
        values = [value.strip('"') for value in values]
        if key == "!Sample_characteristics_ch1":
            key = values[0].split(":")[0]  # "treatment: Untreated"
            values = [value.split(": ", 1)[1] for value in values]
        fields.setdefault(key, values)

samples = pd.DataFrame({
    "sample": fields["!Sample_geo_accession"],
    "cell_line": fields["cell line"],
    "treatment": fields["treatment"],
})
samples.to_csv("samples.tsv", sep="\t", index=False)
print(samples.to_string(index=False))
