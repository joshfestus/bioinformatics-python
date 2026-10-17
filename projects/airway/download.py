"""Download the experiment and the gene sets it will be interpreted with.

The counts and sample descriptions are the GSE52778 files of Chapter 18,
and gene_info gives each NCBI GeneID its symbol.
The gene sets are the Gene Ontology: goa_human.gaf.gz says which genes are
annotated to which terms, and go-basic.obo says what each term means and
which terms contain which.

uv run download.py          writes five files, about 50 MB
"""

import os
import urllib.request

GEO = "https://www.ncbi.nlm.nih.gov/geo/download/"
SERIES = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE52nnn/GSE52778"
GENES = ("https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia"
         "/Homo_sapiens.gene_info.gz")
GOA = "https://ftp.ebi.ac.uk/pub/databases/GO/goa/HUMAN/goa_human.gaf.gz"
OBO = "https://current.geneontology.org/ontology/go-basic.obo"

FILES = [
    (f"{GEO}?type=rnaseq_counts&acc=GSE52778&format=file"
     "&file=GSE52778_raw_counts_GRCh38.p13_NCBI.tsv.gz", "counts.tsv.gz"),
    (f"{SERIES}/matrix/GSE52778_series_matrix.txt.gz",
     "series_matrix.txt.gz"),
    (GENES, "gene_info.tsv.gz"),
    (GOA, "goa_human.gaf.gz"),
    (OBO, "go-basic.obo"),
]


def fetch(url, name):
    headers = {"User-Agent": "python"}
    request = urllib.request.Request(url, headers=headers)
    with urllib.request.urlopen(request, timeout=900) as reply:
        with open(name + ".part", "wb") as out:
            while chunk := reply.read(1 << 16):
                out.write(chunk)
    os.replace(name + ".part", name)
    print(f"{name}, {os.path.getsize(name) / 1e6:.1f} MB")


if __name__ == "__main__":
    for url, name in FILES:
        fetch(url, name)
