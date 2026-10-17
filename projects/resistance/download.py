"""Download the genome and the reference catalog of resistance genes.

The genome is Klebsiella pneumoniae HS11286, a complete public assembly
of a hospital isolate: one chromosome and six plasmids, each with its own
RefSeq accession. The catalog is NCBI's Reference Gene Catalog, the same
set of proteins the AMRFinderPlus tool searches.

As in Project One, each file is written under a temporary name and
renamed once it is complete, so an interrupted download cannot leave a
half file behind.

uv run download.py          writes genome.gb and the catalog files
"""

import os
import urllib.parse
import urllib.request

EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
AMR = ("https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance"
       "/AMRFinderPlus/database/latest")

# One chromosome and six plasmids, in the order NCBI lists them.
REPLICONS = [
    ("NC_016845.1", "chromosome"),
    ("NC_016838.1", "pKPHS1"),
    ("NC_016846.1", "pKPHS2"),
    ("NC_016839.1", "pKPHS3"),
    ("NC_016840.1", "pKPHS4"),
    ("NC_016847.1", "pKPHS5"),
    ("NC_016841.1", "pKPHS6"),
]


def fetch(url, name, data=None):
    with urllib.request.urlopen(url, data=data, timeout=900) as reply:
        with open(name + ".part", "wb") as out:
            while chunk := reply.read(1 << 16):
                out.write(chunk)
    os.replace(name + ".part", name)
    print(f"{name}, {os.path.getsize(name) / 1e6:.1f} MB")


def main():
    accessions = [accession for accession, _ in REPLICONS]
    query = urllib.parse.urlencode({
        "db": "nuccore",
        "id": ",".join(accessions),
        "rettype": "gbwithparts",
        "retmode": "text",
    }).encode()
    fetch(EFETCH, "genome.gb", data=query)
    fetch(f"{AMR}/AMRProt.fa", "amr_proteins.fasta")
    fetch(f"{AMR}/ReferenceGeneCatalog.txt", "amr_catalog.tsv")
    fetch(f"{AMR}/AMRProt-mutation.tsv", "amr_mutations.tsv")
    fetch(f"{AMR}/version.txt", "amr_version.txt")


if __name__ == "__main__":
    main()
