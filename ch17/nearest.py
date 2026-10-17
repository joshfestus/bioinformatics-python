# For each variant between genes, the nearest gene on each side, and
# whether the variant lies upstream of that gene's start.
from bisect import bisect_right

import pandas as pd

VARIANTS = [281_923, 433_359, 2_333_538, 2_999_330, 3_909_807]  # Ch. 16

columns = [
    "seqid",
    "source",
    "type",
    "start",
    "end",
    "score",
    "strand",
    "phase",
    "attributes",
]
gff = pd.read_csv(
    "../ch15/rel606.gff3", sep="\t", comment="#", names=columns
)
genes = gff[gff["type"] == "gene"].sort_values("start")
genes = genes.assign(
    name=genes["attributes"].str.extract(r"Name=([^;]+)", expand=False)
)
starts = genes["start"].tolist()


def describe(gene, position: int) -> str:
    """Name, distance, and whether the position is before its start."""
    if gene.end < position:
        distance = position - gene.end
    else:
        distance = gene.start - position
    begins = gene.start if gene.strand == "+" else gene.end
    upstream = (begins > position) == (gene.strand == "+")
    side = "upstream" if upstream else "downstream"
    return f"{gene['name']} ({gene.strand}) {distance} bp, {side}"


for position in VARIANTS:
    # bisect finds where the position would go in the sorted starts
    after = bisect_right(starts, position)
    left = genes.iloc[after - 1]
    while left.end >= position:  # skip any gene that contains it
        after -= 1
        left = genes.iloc[after - 1]
    right = genes.iloc[bisect_right(starts, position)]
    print(f"{position:9,}")
    print(f"  left:  {describe(left, position)}")
    print(f"  right: {describe(right, position)}")
