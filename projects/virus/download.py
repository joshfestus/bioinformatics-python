"""Download the reference genome and every genome listed in samples.tsv.

Each file is written under a temporary name and renamed once it is
complete, so that a download interrupted halfway cannot leave a truncated
file behind for the next step to read as though it were whole.

uv run download.py          writes reference.fasta, reference.gb and
                            genomes.fasta, about 4 MB, from NCBI
"""

import csv
import os
import urllib.parse
import urllib.request

EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
REFERENCE = "NC_045512.2"  # Wuhan-Hu-1, collected December 2019


def efetch(accessions, rettype):
    data = urllib.parse.urlencode({
        "db": "nuccore",
        "id": ",".join(accessions),
        "rettype": rettype,
        "retmode": "text",
    }).encode()
    with urllib.request.urlopen(EFETCH, data=data, timeout=600) as reply:
        return reply.read().decode()


def main():
    for name, rettype in [("reference.fasta", "fasta"),
                          ("reference.gb", "gb")]:
        with open(name + ".part", "w") as out:
            out.write(efetch([REFERENCE], rettype))
        os.replace(name + ".part", name)

    with open("samples.tsv") as f:
        accessions = [row["accession"]
                      for row in csv.DictReader(f, delimiter="\t")]
    with open("genomes.fasta.part", "w") as out:
        for start in range(0, len(accessions), 40):
            batch = accessions[start:start + 40]
            out.write(efetch(batch, "fasta"))
            print(f"downloaded {start + len(batch)} of {len(accessions)}")
    os.replace("genomes.fasta.part", "genomes.fasta")


if __name__ == "__main__":
    main()
