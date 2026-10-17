# Search NCBI for human HBB records and download one, with Biopython's
# Entrez module.
import os
from pathlib import Path

from Bio import Entrez, SeqIO

Entrez.email = "your.name@example.org"  # NCBI asks for a real address
Entrez.api_key = os.environ.get("NCBI_API_KEY")  # optional, see text
Entrez.tool = "bioinformatics-python-book"

query = "HBB[Gene Name] AND Homo sapiens[Organism] AND refseq[filter]"
with Entrez.esearch(db="nuccore", term=query) as handle:
    found = Entrez.read(handle)
print(f"{found['Count']} RefSeq records match the search")

with Entrez.esummary(
    db="nuccore", id=",".join(found["IdList"])
) as handle:
    summaries = Entrez.read(handle)
for summary in summaries:
    print(f"  {summary['AccessionVersion']:15} {summary['Title'][:52]}")

cache = Path("NM_000518.5.gb")
if not cache.exists():
    with Entrez.efetch(
        db="nuccore", id="NM_000518.5", rettype="gb", retmode="text"
    ) as handle:
        cache.write_text(handle.read())
record = SeqIO.read(cache, "genbank")
print(f"\n{record.id}: {len(record)} bp, {record.description}")
