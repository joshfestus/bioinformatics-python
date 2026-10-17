"""Download the public datasets the book uses, from their sources.

uv run data/fetch.py ch05      the files Chapter 5 needs
uv run data/fetch.py --list    every dataset and where it comes from
"""

import gzip
import sys
import urllib.parse
import urllib.request
from pathlib import Path

HERE = Path(__file__).parent

NCBI_GENOMES = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405"
GRCH38 = "GCF_000001405.40_GRCh38.p14"

UNIPROT_STREAM = "https://rest.uniprot.org/uniprotkb/stream"
HUMAN_REVIEWED = "(reviewed:true) AND (organism_id:9606)"

EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
ENA_READS = "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/003/SRR2584863"
READ_PAIRS = 600_000  # about 39 fold coverage of the E. coli genome

RCSB = "https://files.rcsb.org/download"
ALPHAFOLD = "https://alphafold.ebi.ac.uk/files"

GEO_COUNTS = "https://www.ncbi.nlm.nih.gov/geo/download/"
GEO_SERIES = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE52nnn/GSE52778"

DATASETS = {
    "ch05": [
        (
            "Human genome annotation, GRCh38.p14, NCBI, public domain",
            f"{NCBI_GENOMES}/{GRCH38}/{GRCH38}_feature_table.txt.gz",
            "human_feature_table.txt.gz",
        ),
    ],
    "ch11": [
        (
            "Every reviewed human protein, UniProt Swiss-Prot, CC BY 4.0",
            f"{UNIPROT_STREAM}?compressed=true&format=fasta&query="
            + urllib.parse.quote(HUMAN_REVIEWED),
            "human_proteins.fasta.gz",
        ),
    ],
    "ch15": [
        (
            "E. coli B REL606 genome, NC_012967.1, NCBI, public domain",
            f"{EFETCH}?db=nuccore&id=NC_012967.1&rettype=fasta",
            "rel606.fasta",
        ),
        (
            "SRR2584863 first reads of the first 600,000 pairs",
            f"{ENA_READS}/SRR2584863_1.fastq.gz",
            "SRR2584863_1.fastq.gz",
        ),
        (
            "SRR2584863 second reads of the first 600,000 pairs",
            f"{ENA_READS}/SRR2584863_2.fastq.gz",
            "SRR2584863_2.fastq.gz",
        ),
    ],
    "ch18": [
        (
            "GSE52778 read counts per gene, computed by NCBI, GEO",
            (
                f"{GEO_COUNTS}?type=rnaseq_counts&acc=GSE52778&format=file"
                "&file=GSE52778_raw_counts_GRCh38.p13_NCBI.tsv.gz"
            ),
            "GSE52778_counts.tsv.gz",
        ),
        (
            "GSE52778 sample descriptions, GEO series matrix",
            f"{GEO_SERIES}/matrix/GSE52778_series_matrix.txt.gz",
            "GSE52778_series_matrix.txt.gz",
        ),
    ],
    "ch21": [
        (
            "Human deoxyhemoglobin, PDB 2HHB, RCSB, CC0",
            f"{RCSB}/2HHB.cif",
            "2HHB.cif",
        ),
        (
            "Sickle deoxyhemoglobin, PDB 2HBS, RCSB, CC0",
            f"{RCSB}/2HBS.cif",
            "2HBS.cif",
        ),
        (
            "AlphaFold model of beta globin P68871, CC BY 4.0",
            f"{ALPHAFOLD}/AF-P68871-F1-model_v6.pdb",
            "AF-P68871-F1.pdb",
        ),
    ],
}


def download(url: str, target: Path) -> None:
    """Save a URL to a file; some servers refuse requests with no name."""
    request = urllib.request.Request(
        url, headers={"User-Agent": "bioinformatics-python-book"}
    )
    with urllib.request.urlopen(request) as response:
        target.write_bytes(response.read())


def first_reads(url: str, target: Path, reads: int) -> None:
    """Stream a gzipped FASTQ and keep only its first `reads` records."""
    with (
        urllib.request.urlopen(url) as response,
        gzip.open(response, "rt") as source,
        gzip.open(target, "wt") as output,
    ):
        for number, line in enumerate(source):
            if number == 4 * reads:
                break
            output.write(line)


def fetch(chapter: str) -> None:
    for description, url, name in DATASETS[chapter]:
        target = HERE / name
        if target.exists():
            print(f"already here: {name}")
            continue
        print(f"downloading {description}")
        print(f"  from {url}")
        partial = target.with_suffix(".part")
        if name.startswith("SRR"):
            first_reads(url, partial, READ_PAIRS)
        else:
            download(url, partial)
        partial.rename(target)
        print(f"  saved {name}, {target.stat().st_size:,} bytes")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    if sys.argv[1] == "--list":
        for chapter, files in DATASETS.items():
            for description, url, name in files:
                print(
                    f"{chapter}  {name}\n      {description}\n      {url}"
                )
    else:
        fetch(sys.argv[1])
