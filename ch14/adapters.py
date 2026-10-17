# How many reads contain each of the two common Illumina adapters?
import gzip

from Bio import SeqIO

ADAPTERS = {
    "TruSeq": "AGATCGGAAGAGC",
    "Nextera": "CTGTCTCTTATACACATCT",
}

counts = dict.fromkeys(ADAPTERS, 0)
reads = 0
with gzip.open("reads.fastq.gz", "rt") as handle:
    for record in SeqIO.parse(handle, "fastq"):
        reads += 1
        sequence = str(record.seq)
        for name, adapter in ADAPTERS.items():
            counts[name] += adapter in sequence

for name, count in counts.items():
    print(f"{name:8} {count:6,} of {reads:,} reads ({count / reads:.1%})")
