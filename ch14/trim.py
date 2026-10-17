# Remove adapter and low quality ends from reads, and drop reads left too
# short to be useful.
import gzip

from Bio import SeqIO

ADAPTER = "CTGTCTCTTATACACATCT"  # Nextera, the kit this run used
WINDOW, MIN_QUALITY = 4, 20  # trim where 4 bases average below 20
MIN_LENGTH = 50


def trim(record):
    """The record cut at its adapter and its poor quality tail."""
    cut = str(record.seq).find(ADAPTER)
    if cut != -1:
        record = record[:cut]
    qualities = record.letter_annotations["phred_quality"]
    end = len(qualities)
    for start in range(len(qualities) - WINDOW + 1):
        window = qualities[start : start + WINDOW]
        if sum(window) / WINDOW < MIN_QUALITY:
            end = start
            break
    return record[:end]


kept = dropped = adapters = bases_in = bases_out = 0
with (
    gzip.open("reads.fastq.gz", "rt") as source,
    gzip.open("trimmed.fastq.gz", "wt") as target,
):
    for record in SeqIO.parse(source, "fastq"):
        bases_in += len(record)
        adapters += ADAPTER in record.seq
        trimmed = trim(record)
        if len(trimmed) >= MIN_LENGTH:
            SeqIO.write(trimmed, target, "fastq")
            kept += 1
            bases_out += len(trimmed)
        else:
            dropped += 1

print(f"reads with adapter     {adapters:7,}")
print(f"reads kept             {kept:7,}")
print(f"reads dropped          {dropped:7,}")
print(f"bases kept             {bases_out / bases_in:7.1%}")
