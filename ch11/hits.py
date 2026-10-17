# Read BLAST's tabular output (-outfmt 6), keeping the trustworthy hits.
import pandas as pd
from Bio import SeqIO

COLUMNS = [
    "query",
    "subject",
    "identity",
    "length",
    "mismatches",
    "gap_opens",
    "q_start",
    "q_end",
    "s_start",
    "s_end",
    "evalue",
    "bits",
]
CUTOFF = 1e-3

hits = pd.read_csv("hits.tsv", sep="\t", names=COLUMNS)
print(
    f"{len(hits)} hits reported, {(hits['evalue'] < CUTOFF).sum()} "
    f"with an E-value below {CUTOFF}"
)

names = SeqIO.index("human_proteins.fasta", "fasta")
good = hits[hits["evalue"] < CUTOFF].copy()
good["protein"] = [
    names[s].description.split(" ", 1)[1].split(" OS=")[0][:28]
    for s in good["subject"]
]
good["entry"] = good["subject"].str.split("|").str[2]
columns = ["entry", "protein", "identity", "length", "evalue", "bits"]
print(
    good[columns].to_string(
        index=False, formatters={"evalue": "{:.1e}".format}
    )
)
