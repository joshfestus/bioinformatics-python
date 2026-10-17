# Quality summary of the first 2,000 reads of SRR2584863, an Illumina run
# of E. coli from the long term evolution experiment (NCBI SRA).
import gzip

from Bio import SeqIO

reads = bases = high_quality = poor_reads = 0

with gzip.open("reads.fastq.gz", "rt") as handle:
    for record in SeqIO.parse(handle, "fastq"):
        qualities = record.letter_annotations["phred_quality"]
        reads += 1
        bases += len(qualities)
        high_quality += sum(q >= 30 for q in qualities)
        if sum(qualities) / len(qualities) < 20:
            poor_reads += 1

print(f"reads                    {reads:8,}")
print(f"bases                    {bases:8,}")
print(f"bases with quality >= 30 {high_quality / bases:8.1%}")
print(f"reads averaging below 20 {poor_reads:8,}")
