# How deeply the reads cover each base of the REL606 genome, and which
# genes lie in stretches no read covers.
import numpy as np
import pandas as pd
import pysam

bam = pysam.AlignmentFile("reads.bam")
contig, length = bam.references[0], bam.lengths[0]

# count_coverage gives four arrays, the A, C, G and T counts at each base
depth = np.array(
    bam.count_coverage(contig, 0, length, quality_threshold=0)
).sum(axis=0)
print(f"{contig}: {length:,} bases")
print(f"  mean depth {depth.mean():.1f}, median {np.median(depth):.0f}")
print(f"  bases with no reads: {(depth == 0).sum():,}")

# Runs of zero depth: where the depth starts and stops being zero
edges = np.diff(np.concatenate([[0], depth == 0, [0]]).astype(int))
starts, ends = np.where(edges == 1)[0], np.where(edges == -1)[0]
runs = sorted(zip(ends - starts, starts, ends), reverse=True)

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
gff = pd.read_csv("rel606.gff3", sep="\t", comment="#", names=columns)
genes = gff[gff["type"] == "gene"].copy()
genes["name"] = genes["attributes"].str.extract(r"Name=([^;]+)")

print("longest stretches with no reads:")
for size, start, end in runs[:3]:
    inside = genes[(genes["start"] <= end) & (genes["end"] > start)]
    names = ", ".join(inside["name"])
    print(f"  {start + 1:,} to {end:,}, {size:,} bases")
    if names:
        print(f"    genes: {names}")

np.save("depth.npy", depth.astype(np.int32))  # for the figures
