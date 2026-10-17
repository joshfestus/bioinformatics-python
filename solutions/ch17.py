"""Solutions to the exercises of Chapter 17.

Run from this folder: uv run ch17.py
Reads Chapter 17's saved transcript and ClinVar file. Exercise 17.3
fetches a wider region from ClinVar and needs the network.
"""

import json
from collections import Counter
from pathlib import Path

import pysam

CH16, CH17 = Path("../ch16"), Path("../ch17")
CLINVAR = ("https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38"
           "/clinvar.vcf.gz")
GENE = (5_225_464, 5_227_071)  # HBB on chromosome 11, one based


def transcript():
    return json.loads((CH17 / "hbb_transcript.json").read_text())


def exons():
    """Exon intervals, one based and inclusive, highest first."""
    data = transcript()
    return sorted(((exon["start"], exon["end"])
                   for exon in data["Exon"]), reverse=True)


def variants(path=CH17 / "hbb_clinvar.vcf"):
    with pysam.VariantFile(path) as vcf:
        for record in vcf:
            info = record.info
            significance = ",".join(info.get("CLNSIG", ()))
            yield record, significance


def exercise_17_1():
    print("17.1  HBB's three exons as BED, converting from the one")
    print("      based inclusive coordinates Ensembl gives:\n")
    print("       chr11\t5226929\t5227071")
    print("       chr11\t5226576\t5226799")
    print("       chr11\t5225463\t5225726\n")
    print(f"       {'exon':6} {'Ensembl':>22} {'BED':>22} {'length':>7}")
    total = 0
    for number, (start, end) in enumerate(exons(), 1):
        length = end - start + 1
        total += length
        print(f"       {number:<6} {start:>10,} to {end:>8,} "
              f"{start - 1:>10,} to {end:>8,} {length:>7}")
    print(f"      The three add to {total:,} bases, which is the")
    print("      length of the messenger RNA, as Chapter 7 found by")
    print("      a different route.")
    print("      Only the start moves. BED counts from zero and")
    print("      excludes the end, Ensembl counts from one and")
    print("      includes it, so subtracting one from the start")
    print("      converts between them and the end stays as it is.")
    print("      Every length is unchanged by the conversion, which")
    print("      is the check worth running: if a length moves, the")
    print("      conversion is wrong.")


def distance_to_exon(position):
    """Bases from a position to the nearest exon edge, 0 if inside."""
    best = None
    for start, end in exons():
        if start <= position <= end:
            return 0
        gap = min(abs(position - start), abs(position - end))
        best = gap if best is None else min(best, gap)
    return best


def exercise_17_2():
    pathogenic, benign = [], []
    for record, significance in variants():
        if distance_to_exon(record.pos) == 0:
            continue
        if not (GENE[0] <= record.pos <= GENE[1]):
            continue
        gap = distance_to_exon(record.pos)
        if "Pathogenic" in significance:
            pathogenic.append(gap)
        elif "Benign" in significance:
            benign.append(gap)
    print("\n17.2  Intron variants, by distance to the nearest exon:")
    for label, found in (("pathogenic", pathogenic),
                         ("likely benign", benign)):
        if not found:
            continue
        close = sum(1 for gap in found if gap <= 10)
        print(f"       {label:14} {len(found):4} variants, "
              f"{close:3} within ten bases ({close / len(found):.0%})")
    print("      The splice sites are the answer. An intron is")
    print("      removed by machinery that recognizes a few bases at")
    print("      each end, so a change there can leave the intron in")
    print("      the message or cut the exon out with it, and either")
    print("      ruins the protein. Ten bases further in, the intron")
    print("      is mostly spacer and a change does nothing.")
    print("      That is why the pathogenic variants crowd the edges")
    print("      and the benign ones are spread through the middle.")
    print("      It is also why a variant caller's output needs")
    print("      intervals: the same change is harmless or disastrous")
    print("      depending on a distance of a few bases.")


def exercise_17_3():
    print("\n17.3  Extending the region 1,000 bases past the gene's")
    print("      high end, which on this minus strand gene is where")
    print("      the promoter sits:")
    target = Path("hbb_upstream.vcf")
    region = (GENE[0], GENE[1] + 1000)
    try:
        with pysam.VariantFile(CLINVAR) as vcf:
            found = list(vcf.fetch("11", region[0] - 1, region[1]))
            header = vcf.header
    except (OSError, ValueError) as problem:
        print(f"      Could not reach ClinVar ({problem}).")
        return
    with pysam.VariantFile(target, "w", header=header) as out:
        for record in found:
            out.write(record)
    upstream = [r for r in found if r.pos > GENE[1]]
    pathogenic = [r for r in upstream
                  if "Pathogenic" in ",".join(r.info.get("CLNSIG", ()))]
    print(f"       {len(found):,} records in the wider region, "
          f"{len(upstream)} of them upstream.")
    print(f"       {len(pathogenic)} of those are pathogenic:")
    for record in pathogenic[:8]:
        print(f"         {record.pos:,}, {record.pos - GENE[1]} bases "
              f"past the gene, {record.ref}>{record.alts[0]}")
    print("      These are promoter variants. The gene's start is at")
    print("      its high coordinate because HBB is on the minus")
    print("      strand, so the sequence that switches it on lies")
    print("      above 5,227,071 and outside every region Chapter 16")
    print("      asked for. A variant there makes too little beta")
    print("      globin without changing a single codon.")
    print("      The lesson is about the region, not the gene: ask")
    print("      only for the gene and you will never see the")
    print("      variants that control it, and nothing will warn you")
    print("      that a class of disease is missing from your table.")
    target.unlink()


def exercise_17_4():
    print("\n17.4  Testing a deletion by its whole interval.")
    print("      Chapter 17 found six pathogenic variants that fall")
    print("      outside HBB when each is judged by its POS alone.")
    outside = []
    for record, significance in variants():
        if "Pathogenic" not in significance:
            continue
        if GENE[0] <= record.pos <= GENE[1]:
            continue
        span = (record.pos, record.pos + len(record.ref) - 1)
        outside.append((record.pos, len(record.ref), span))
    print(f"      There are {len(outside)} of them, and every one is a")
    print("      deletion that starts below the gene and reaches up")
    print("      into it:\n")
    print(f"       {'POS':>11} {'length':>8}  reaches")
    for position, length, span in sorted(outside):
        touched = []
        for number, (start, end) in enumerate(exons(), 1):
            if span[0] <= end and span[1] >= start:
                touched.append(f"exon {number}")
        if span[1] > GENE[1]:
            touched.append("past the promoter")
        print(f"       {position:>11,} {length:>8,}  "
              f"{', '.join(touched) or 'no exon'}")
    print("\n      Judged by POS they are outside the gene and")
    print("      therefore nothing to do with it. Judged by their")
    print("      whole interval, every one of them removes at least")
    print("      the third exon, and two remove all three exons and")
    print("      reach past the promoter as well. These are the")
    print("      deletions that cause the severest thalassemias.")
    print("      Note that the longest, 7,722 bases, is not one of")
    print("      those two: it starts further down and runs out")
    print("      before the first exon. Length alone does not tell")
    print("      you what a deletion costs; where it lands does.")
    print("      A single coordinate cannot describe a thing that")
    print("      has length. For anything but a substitution, test")
    print("      the interval, and the only cost is remembering that")
    print("      the interval runs from POS to POS plus the length of")
    print("      REF, less one.")


def exercise_17_5():
    print("\n17.5  All 25 confident E. coli variants, inside or out:")
    import pandas as pd
    columns = ["seqid", "source", "type", "start", "end", "score",
               "strand", "phase", "attributes"]
    gff = pd.read_csv("../ch15/rel606.gff3", sep="\t", comment="#",
                      names=columns)
    cds = gff[gff["type"] == "CDS"]
    spans = list(zip(cds["start"], cds["end"], cds["attributes"]))
    inside = outside = 0
    with pysam.VariantFile(CH16 / "haploid.vcf.gz") as vcf:
        calls = [r for r in vcf
                 if r.qual and r.qual >= 30 and r.info.get("DP", 0) >= 10]
    for record in calls:
        hit = next((attributes for start, end, attributes in spans
                    if start <= record.pos <= end), None)
        if hit:
            inside += 1
        else:
            outside += 1
    print(f"       {len(calls)} confident calls: {inside} inside a")
    print(f"       coding sequence, {outside} between genes.")
    print("      Chapter 16's effects.py reported 20 inside a gene,")
    print("      13 of them missense, and 5 between genes, which is")
    print("      the same split: the counts match because both")
    print("      programs ask the same question of the same")
    print("      annotation.")
    print("      That agreement is the point of the exercise. Two")
    print("      programs written at different times, one walking")
    print("      the CDS lines and one using bisect on sorted gene")
    print("      starts, should give the same answer, and when they")
    print("      do not it is nearly always an off by one at an")
    print("      interval edge.")


if __name__ == "__main__":
    exercise_17_1()
    exercise_17_2()
    exercise_17_3()
    exercise_17_4()
    exercise_17_5()
