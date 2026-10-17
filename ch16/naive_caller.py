# The simplest possible variant caller: at each position, report the most
# common base if enough reads cover it and enough of them agree.
import sys

import numpy as np
import pysam
from Bio import SeqIO

MIN_DEPTH = int(sys.argv[1]) if len(sys.argv) > 1 else 10
MIN_SHARE = float(sys.argv[2]) if len(sys.argv) > 2 else 0.8

reference = np.array(
    list(SeqIO.read("../data/rel606.fasta", "fasta").seq)
)
bam = pysam.AlignmentFile("../ch15/reads.bam")
counts = np.array(
    bam.count_coverage(
        bam.references[0], 0, len(reference), quality_threshold=0
    )
)

depth = counts.sum(axis=0)
top_base = np.array(list("ACGT"))[counts.argmax(axis=0)]
share = counts.max(axis=0) / np.maximum(depth, 1)
called = (depth >= MIN_DEPTH) & (share >= MIN_SHARE)
called &= top_base != reference

positions = np.where(called)[0]
print(
    f"depth >= {MIN_DEPTH}, share >= {MIN_SHARE}: {len(positions)} calls"
)
rows = [
    f"{p + 1}\t{reference[p]}\t{top_base[p]}\t{depth[p]}\n"
    for p in positions
]
with open("naive.tsv", "w") as out:
    out.writelines(rows)
