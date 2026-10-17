# Turn GEO's series matrix for GSE52778 into a table of samples.
import gzip

import pandas as pd

fields = {}
with gzip.open("../data/GSE52778_series_matrix.txt.gz", "rt") as handle:
    for line in handle:
        if not line.startswith("!Sample_"):
            continue
        key, *values = line.rstrip("\n").split("\t")
        values = [v.strip('"') for v in values]
        if key == "!Sample_characteristics_ch1":
            key = values[0].split(":")[0]  # "treatment: Untreated"
            values = [v.split(": ", 1)[1] for v in values]
        fields.setdefault(key, values)

samples = pd.DataFrame(
    {
        "sample": fields["!Sample_geo_accession"],
        "cell_line": fields["cell line"],
        "treatment": fields["treatment"],
    }
)
samples.to_csv("samples.tsv", sep="\t", index=False)
print(samples.to_string(index=False))
