# BLAST human beta globin against Swiss-Prot on NCBI's servers. The search
# runs once; the result is saved and read from the file afterwards.
from pathlib import Path

from Bio import Blast, SeqIO

Blast.email = "your.name@example.org"  # NCBI asks for a real address
SAVED = Path("hbb_swissprot.xml")

if not SAVED.exists():
    beta = SeqIO.read("hbb.fasta", "fasta").seq
    print("submitting to NCBI; this usually takes a few minutes")
    with Blast.qblast(
        "blastp", "swissprot", beta, hitlist_size=50
    ) as answer:
        SAVED.write_bytes(answer.read())

with SAVED.open("rb") as handle:
    record = Blast.read(handle)

print(f"{len(record)} hits from {record.db}")
for hit in record[:12]:
    best = hit[0]
    species = hit.target.description.split("[")[-1].rstrip("]")
    print(
        f"  {hit.target.id:14} {best.annotations['evalue']:9.1e}"
        f"  {species[:40]}"
    )
