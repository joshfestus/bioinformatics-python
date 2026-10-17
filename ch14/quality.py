# Quality at each position along the reads of a FASTQ file, as a box plot.
import gzip
import sys

import matplotlib
from Bio import SeqIO

matplotlib.use("svg")
import matplotlib.pyplot as plt

path, picture = sys.argv[1], sys.argv[2]

by_position = []  # by_position[i] holds every quality seen at position i
reads = 0
with gzip.open(path, "rt") as handle:
    for record in SeqIO.parse(handle, "fastq"):
        reads += 1
        for i, q in enumerate(record.letter_annotations["phred_quality"]):
            if i == len(by_position):
                by_position.append([])
            by_position[i].append(q)

means = [sum(qs) / len(qs) for qs in by_position]
print(f"{reads:,} reads, up to {len(by_position)} bases long")
for start in (0, 49, 99, 129, 139, 149):
    if start < len(means):
        mean = means[start]
        print(f"  mean quality at position {start + 1:3}: {mean:5.1f}")

fig, ax = plt.subplots(figsize=(5.6, 2.8))
ax.boxplot(
    by_position[::5],
    positions=range(1, len(by_position) + 1, 5),
    widths=3,
    showfliers=False,
    patch_artist=True,
    boxprops={"facecolor": "#BBBBBB", "color": "#444444"},
    medianprops={"color": "#222222"},
)
ax.axhline(20, color="#777777", linestyle="--", linewidth=0.8)
ax.set_xticks(range(0, len(by_position) + 1, 25))
ax.set_xticklabels(range(0, len(by_position) + 1, 25))
ax.set_ylim(0, 42)
ax.set_xlabel("Position in read")
ax.set_ylabel("Phred quality")
fig.savefig(picture, bbox_inches="tight")
print(f"saved {picture}")
