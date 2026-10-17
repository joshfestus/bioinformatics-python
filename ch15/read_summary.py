# What the aligner decided about the first 100,000 reads.
from collections import Counter

import pysam

kinds = Counter()
clipped = 0
with pysam.AlignmentFile("reads.bam") as bam:
    for number, read in enumerate(bam.fetch(until_eof=True)):
        if number == 100_000:
            break
        if read.is_secondary or read.is_supplementary:
            continue
        if read.is_unmapped:
            kinds["unmapped"] += 1
        elif read.mapping_quality < 20:
            kinds["mapped, quality below 20"] += 1
        elif read.is_proper_pair:
            kinds["mapped, properly paired"] += 1
        else:
            kinds["mapped, pair looks unusual"] += 1
        if read.cigartuples and any(
            op == 4 for op, _ in read.cigartuples
        ):
            clipped += 1

for kind, count in kinds.most_common():
    print(f"{kind:28} {count:6,}")
print(f"{'soft clipped at one end':28} {clipped:6,}")
