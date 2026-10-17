# The motif around the proximal histidine, the residue that holds the
# heme iron: a position weight matrix and a sequence logo.
import math

import logomaker
import matplotlib
from Bio import AlignIO, motifs

matplotlib.use("svg")
import matplotlib.pyplot as plt

FIRST, LAST = 89, 100  # alignment columns, counted from 1; no gaps here
AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"

alignment = AlignIO.read("globins_aligned.fasta", "fasta")
window = [str(r.seq[FIRST - 1 : LAST]) for r in alignment]
motif = motifs.create(window, alphabet=AMINO_ACIDS)
print(f"{len(window)} instances, consensus {motif.consensus}")

frequencies = motif.counts.normalize()
print("column  most common  share  information (bits)")
for i in range(len(motif)):
    shares = {aa: frequencies[aa][i] for aa in AMINO_ACIDS}
    entropy = -sum(p * math.log2(p) for p in shares.values() if p > 0)
    top = max(shares, key=shares.get)
    information = math.log2(len(AMINO_ACIDS)) - entropy
    share = shares[top]
    print(f"{FIRST + i:6}  {top:>11}  {share:5.2f}  {information:8.2f}")

counts = logomaker.alignment_to_matrix(window)
bits = logomaker.transform_matrix(
    counts, from_type="counts", to_type="information", pseudocount=0
)
fig, ax = plt.subplots(figsize=(5.0, 2.2))
logomaker.Logo(bits, ax=ax, color_scheme="black", font_name="DejaVu Sans")
ax.set_xticks(range(len(window[0])))
ax.set_xticklabels(range(FIRST, LAST + 1))
ax.set_ylim(0, math.log2(len(AMINO_ACIDS)))
ax.set_ylabel("bits")
fig.savefig("proximal-his-logo.svg", bbox_inches="tight")
print("logo saved as proximal-his-logo.svg")
