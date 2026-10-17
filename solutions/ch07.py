"""Solutions to the exercises of Chapter 7.

Run from this folder: uv run ch07.py
Figures are saved here as PNG files rather than shown on screen.
"""

import gzip
from pathlib import Path

import matplotlib
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

CH07 = Path("../ch07")
CH06 = Path("../ch06")
COLUMNS = ["seqid", "source", "type", "start", "end", "score", "strand",
           "phase", "attributes"]


def reads():
    with gzip.open(CH07 / "reads.fastq.gz", "rt") as handle:
        yield from SeqIO.parse(handle, "fastq")


def exercise_7_1():
    totals, counts = [], []
    for record in reads():
        qualities = record.letter_annotations["phred_quality"]
        while len(totals) < len(qualities):
            totals.append(0)
            counts.append(0)
        for position, quality in enumerate(qualities):
            totals[position] += quality
            counts[position] += 1
    means = [total / count for total, count in zip(totals, counts)]
    figure, axes = plt.subplots(figsize=(6, 3.2))
    axes.plot(range(1, len(means) + 1), means, color="#333333")
    axes.set_xlabel("Position in the read")
    axes.set_ylabel("Mean quality")
    axes.set_ylim(0, 42)
    figure.savefig("ch07-quality.png", dpi=110, bbox_inches="tight")
    plt.close(figure)
    best = max(range(len(means)), key=lambda i: means[i])
    falls = next(i for i in range(best, len(means)) if means[i] < 30)
    print(f"7.1  {len(means)} positions, mean quality {means[0]:.1f} at "
          f"the first base,")
    print(f"     highest {means[best]:.1f} at position {best + 1}, "
          f"{means[-1]:.1f} at the last.")
    print(f"     It first drops below 30 at position {falls + 1}.")
    print("     Quality is low for the first few bases, climbs quickly,")
    print("     then falls away towards the end of the read. That shape")
    print("     is what sequencing by synthesis does: the chemistry")
    print("     falls out of step as the cycles go on. Saved as")
    print("     ch07-quality.png.")


def exercise_7_2():
    shares = []
    for record in reads():
        sequence = str(record.seq).upper()
        called = sum(sequence.count(base) for base in "ACGT")
        if called:
            shares.append(
                (sequence.count("G") + sequence.count("C")) / called)
    series = pd.Series(shares)
    figure, axes = plt.subplots(figsize=(6, 3.2))
    axes.hist(series, bins=40, color="#777777", edgecolor="#222222",
              linewidth=0.4)
    axes.axvline(0.51, color="#222222", lw=1, ls=(0, (4, 3)))
    axes.set_xlabel("GC content of a read")
    axes.set_ylabel("Reads")
    figure.savefig("ch07-gc.png", dpi=110, bbox_inches="tight")
    plt.close(figure)
    print(f"\n7.2  {len(series):,} reads, mean GC {series.mean():.1%}, "
          f"median {series.median():.1%}")
    print(f"     The middle half lies between {series.quantile(.25):.1%} "
          f"and {series.quantile(.75):.1%}.")
    print("     The histogram is centered close to the 51 percent the")
    print("     genome has, which is the answer the exercise wants, but")
    print("     it is wide: a 150 base window of a real genome varies a")
    print("     lot, and the spread here is that variation rather than")
    print("     a fault. A second peak away from the genome's value is")
    print("     what contamination looks like. Saved as ch07-gc.png.")


def exercise_7_3():
    record = SeqIO.read(CH07 / "hbb.gb", "genbank")
    exons = [f for f in record.features if f.type == "exon"]
    lengths = [len(f.location) for f in exons]
    print(f"\n7.3  {len(exons)} exons in {record.id}:")
    for number, feature in enumerate(exons, 1):
        start = int(feature.location.start) + 1
        print(f"       exon {number}: {start:4} to "
              f"{int(feature.location.end):4}, {len(feature.location):3} "
              "bases")
    print(f"     They add to {sum(lengths)} bases, and the record itself")
    print(f"     is {len(record)} bases.")
    print("     They match, because this record is the messenger RNA,")
    print("     which is what is left after the introns are cut out. A")
    print("     genomic record of the same gene would be 1,600 bases")
    print("     with the introns still in it.")


def exercise_7_4():
    region = SeqIO.read(CH06 / "lac_region.gb", "genbank")
    # headers read ">NP_414876.1 lacA galactoside O-acetyltransferase",
    # so the gene name is the second word
    expected = {}
    for record in SeqIO.parse(CH07 / "lac_proteins.fasta", "fasta"):
        words = record.description.split()
        if len(words) > 1:
            expected[words[1]] = str(record.seq)
    print("\n7.4  Cutting each gene out of the region and translating:")
    for feature in region.features:
        name = feature.qualifiers.get("gene", [""])[0]
        if feature.type != "CDS" or name not in expected:
            continue
        dna = feature.extract(region.seq)
        plain = str(dna.translate(table=11))[:-1]
        as_cds = str(dna.translate(table=11, cds=True))
        print(f"       {name}: {len(dna)} bases, {len(as_cds)} amino "
              f"acids, begins {plain[0]} translated plainly and "
              f"{as_cds[0]} as a CDS")
        print(f"            against the real protein: "
              f"{'matches' if as_cds == expected[name] else 'DIFFERS'}")
    print("     All four match, but only when the translation is told")
    print("     it is looking at a coding sequence. lacA begins TTG and")
    print("     lacI begins GTG, and a plain translation renders those")
    print("     as leucine and valine, which is what those codons mean")
    print("     in the middle of a gene. The cell puts methionine first")
    print("     whichever start codon it reads, and `cds=True` is how")
    print("     you tell Biopython that the first codon is a start.")
    print("     Two of the four proteins are wrong in their first")
    print("     letter without it, and nothing warns you.")
    print("     Biopython's extract also does the strand work: a gene")
    print("     on the minus strand comes back already reverse")
    print("     complemented and reading the right way. By hand you")
    print("     would subtract the region's offset, slice, then reverse")
    print("     complement only the minus strand genes.")


def exercise_7_5():
    gff = pd.read_csv(CH07 / "ecoli_k12.gff3.gz", sep="\t", comment="#",
                      names=COLUMNS)
    genes = gff[gff["type"] == "gene"].copy()
    genes["length"] = genes["end"] - genes["start"] + 1
    genes["biotype"] = genes["attributes"].str.extract(
        r"gene_biotype=([^;]+)")
    genes["name"] = genes["attributes"].str.extract(r"Name=([^;]+)")
    coding = genes[genes["biotype"] == "protein_coding"]
    cds = gff[gff["type"] == "CDS"].copy()
    cds["parent"] = cds["attributes"].str.extract(r"Parent=([^;]+)")
    cds["product"] = cds["attributes"].str.extract(r"product=([^;]+)")
    products = dict(zip(cds["parent"], cds["product"]))
    genes["id"] = genes["attributes"].str.extract(r"ID=([^;]+)")
    longest = coding.loc[coding["length"].idxmax()]
    shortest = coding.loc[coding["length"].idxmin()]
    print(f"\n7.5  {len(coding):,} protein coding genes.")
    for label, gene in (("longest", longest), ("shortest", shortest)):
        row = genes.loc[gene.name]
        product = products.get(row["id"], "no product line")
        print(f"       {label}: {row['name']}, {row['length']:,} bases")
        print(f"            {product.replace('%2C', ',')[:58]}")
    biggest = longest["length"] // 3 - 1
    smallest = shortest["length"] // 3 - 1
    print(f"     That is {biggest:,} amino acids against {smallest}.")
    print("     The longest is a surface protein that anchors in the")
    print("     outer membrane, and the shortest is a peptide of eight")
    print("     amino acids. Both are real annotations, and the short")
    print("     one is exactly the gene an ORF finder with a 100 amino")
    print("     acid threshold can never see, which is Exercise 6.3's")
    print("     cost made concrete.")


if __name__ == "__main__":
    exercise_7_1()
    exercise_7_2()
    exercise_7_3()
    exercise_7_4()
    exercise_7_5()
