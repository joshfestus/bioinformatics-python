# Download a set of globin proteins from UniProt, once, as FASTA.
import urllib.request
from pathlib import Path

from Bio import SeqIO

# Beta globins from twelve vertebrates, then alpha globins and myoglobin
ACCESSIONS = [
    "P68871",
    "P68873",
    "P02088",
    "P02091",
    "P02070",
    "P02067",
    "P02062",
    "P60524",
    "P02111",
    "P02112",
    "P02132",
    "Q90486",
    "P69905",
    "P01942",
    "P01994",
    "P02144",
]
UNIPROT = "https://rest.uniprot.org/uniprotkb/accessions"
URL = UNIPROT + "?format=fasta&accessions="

target = Path("globins.fasta")
if target.exists():
    print(f"using the copy already in {target}")
else:
    with urllib.request.urlopen(URL + ",".join(ACCESSIONS)) as response:
        target.write_bytes(response.read())
    print(f"downloaded {target}")

for record in SeqIO.parse(target, "fasta"):
    entry = record.id.split("|")[2]
    species = record.description.split("OS=")[1].split(" OX=")[0]
    print(f"  {entry:11} {len(record):4} aa  {species}")
