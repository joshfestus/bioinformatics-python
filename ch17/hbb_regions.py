# Where in HBB do ClinVar's variants fall, and does that depend on whether
# they cause disease? Uses the canonical transcript, HBB-201.
import json
import urllib.request
from collections import Counter
from pathlib import Path

import pysam

ENSEMBL = "https://rest.ensembl.org/lookup/id/ENST00000335295"
SAVED = Path("hbb_transcript.json")

if not SAVED.exists():
    request = urllib.request.Request(
        ENSEMBL + "?expand=1;utr=1",
        headers={"Accept": "application/json"},
    )
    with urllib.request.urlopen(request) as response:
        SAVED.write_bytes(response.read())
transcript = json.loads(SAVED.read_text())

# Ensembl counts from 1 and includes both ends; convert every interval to
# Python's convention, zero based with the end excluded, once, here.
exons = [(e["start"] - 1, e["end"]) for e in transcript["Exon"]]
coding = (
    transcript["Translation"]["start"] - 1,
    transcript["Translation"]["end"],
)
gene = (transcript["start"] - 1, transcript["end"])


def inside(position: int, interval: tuple[int, int]) -> bool:
    """Is a zero based position inside a half open interval?"""
    start, end = interval
    return start <= position < end


def region(position: int) -> str:
    if not inside(position, gene):
        return "outside the gene"
    if not any(inside(position, exon) for exon in exons):
        return "intron"
    if inside(position, coding):
        return "coding exon"
    return "untranslated exon"


counts = Counter()
for record in pysam.VariantFile("hbb_clinvar.vcf"):
    significance = record.info.get("CLNSIG", ("none",))[0]
    if significance in ("Pathogenic", "Likely_benign"):
        counts[significance, region(record.pos - 1)] += 1

for significance in ("Pathogenic", "Likely_benign"):
    print(significance.replace("_", " "))
    for place in (
        "coding exon",
        "untranslated exon",
        "intron",
        "outside the gene",
    ):
        print(f"  {place:18} {counts[significance, place]:4}")
