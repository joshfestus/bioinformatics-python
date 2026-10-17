# Find the membrane crossing helices of human rhodopsin (UniProt P08100)
# from its sequence alone, with a Kyte-Doolittle hydropathy plot.
import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.SeqUtils.ProtParamData import kd

WINDOW = (
    19  # about the number of residues a helix needs to cross a membrane
)
THRESHOLD = 1.6  # Kyte and Doolittle's cut-off for a membrane segment
MAX_GAP = 10  # dips shorter than this do not split a segment

rhodopsin = str(SeqIO.read("rhodopsin.fasta", "fasta").seq)
scores = ProteinAnalysis(rhodopsin).protein_scale(kd, WINDOW)
centers = range(WINDOW // 2 + 1, len(rhodopsin) - WINDOW // 2 + 1)
above = [p for p, score in zip(centers, scores) if score > THRESHOLD]

# Join positions above the threshold into segments, bridging short dips
segments = [[above[0], above[0]]]
for position in above[1:]:
    if position - segments[-1][1] <= MAX_GAP:
        segments[-1][1] = position
    else:
        segments.append([position, position])

known = pd.read_csv("rhodopsin_helices.tsv", sep="\t", comment="#")
print(f"{len(rhodopsin)} residues, {len(segments)} predicted segments")
for first, last in segments:
    overlap = known[(known["start"] <= last) & (known["end"] >= first)]
    helix = overlap["helix"].iloc[0] if len(overlap) else "none"
    print(f"  centers {first:3} to {last:3}   UniProt helix {helix}")

for helix in known.itertuples():
    inside = [
        s
        for p, s in zip(centers, scores)
        if helix.start <= p <= helix.end
    ]
    print(
        f"helix {helix.helix}: {helix.start:3} to {helix.end:3}, "
        f"highest score {max(inside):.2f}"
    )
