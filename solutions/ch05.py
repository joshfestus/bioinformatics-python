"""Solutions to the exercises of Chapter 5.

Run from this folder: uv run ch05.py
Needs ../ch05/genes.tsv, which Chapter 5's gene_table.py writes.
Figures are saved here as PNG files rather than shown on screen.
"""

from pathlib import Path

import matplotlib
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

TABLE = Path("../ch05/genes.tsv")
LENGTHS = {"1": 248_956_422, "19": 58_617_616}


def genes():
    if not TABLE.exists():
        raise SystemExit("run ../ch05/gene_table.py first")
    return pd.read_csv(TABLE, sep="\t", dtype={"chromosome": str})


def exercise_5_1(table):
    coding = table[table["kind"] == "protein_coding"]
    print("5.1  Protein coding genes per million bases:")
    density = {}
    for name, bases in LENGTHS.items():
        count = (coding["chromosome"] == name).sum()
        density[name] = count / (bases / 1e6)
        print(f"       chromosome {name:2}: {count:5} genes in "
              f"{bases / 1e6:6.1f} Mb, {density[name]:.1f} per Mb")
    times = density["19"] / density["1"]
    print(f"     Chromosome 19 is {times:.1f} times denser.")
    print("     Gene density varies several fold across the genome, so")
    print("     a count of genes means little without the length.")


def exercise_5_2(table):
    kinds = ["protein_coding", "lncRNA", "pseudogene"]
    shades = ["#222222", "#777777", "#BBBBBB"]
    figure, axes = plt.subplots(figsize=(6, 4))
    print("\n5.2  Median length of each kind, in bases:")
    for kind, shade in zip(kinds, shades):
        lengths = table.loc[table["kind"] == kind, "length"].dropna()
        if lengths.empty:
            print(f"       {kind}: no rows with this exact name")
            continue
        axes.hist(lengths, bins=60, range=(0, 200_000), color=shade,
                  histtype="step", linewidth=1.2, label=kind)
        print(f"       {kind:15} {len(lengths):6,} genes, "
              f"median {lengths.median():9,.0f}")
    axes.set_xlabel("Gene length in bases")
    axes.set_ylabel("Genes")
    axes.legend(frameon=False)
    figure.savefig("ch05-kinds.png", dpi=110, bbox_inches="tight")
    plt.close(figure)
    print("     Pseudogenes are the shortest: most are decayed copies of")
    print("     a messenger RNA, so they carry no introns. Saved as")
    print("     ch05-kinds.png.")


def exercise_5_3(table):
    print("\n5.3  gene_table.py keeps `strand` by adding it to the list")
    print("     of columns it selects:\n")
    print('       genes = genes[["symbol", "GeneID", "chromosome",')
    print('                      "start", "end", "class", "strand"]]\n')
    if "strand" not in table.columns:
        print("     The saved table here has no strand column, so rerun")
        print("     gene_table.py with the change to try it. In the full")
        print("     annotation the split is close to even on every")
        print("     chromosome, a few percent either way, because there")
        print("     is no reason for a gene to prefer a strand.")
        return
    coding = table[table["kind"] == "protein_coding"]
    counts = coding.groupby(["chromosome", "strand"]).size().unstack()
    share = counts["+"] / counts.sum(axis=1)
    print(f"     Plus strand share runs {share.min():.1%} to "
          f"{share.max():.1%} across chromosomes.")


def exercise_5_4(table):
    coding = table[(table["kind"] == "protein_coding")
                   & table["protein_length"].notna()].copy()
    coding["coding_share"] = ((coding["protein_length"] + 1) * 3
                              / coding["length"])
    top = coding.nlargest(20, "coding_share")
    impossible = coding[coding["coding_share"] > 1]
    print("\n5.4  The 20 protein coding genes with the highest share of")
    print("     their length spent coding:")
    for start in range(0, 20, 5):
        names = top["symbol"].astype(str).iloc[start:start + 5]
        print(f"       {', '.join(names)}")
    print(f"     Their share is 99% or more, against a median of "
          f"{coding['coding_share'].median():.1%}")
    print("     for all protein coding genes. Olfactory receptors, the")
    print("     OR symbols, take most of the list, with the late")
    print("     cornified envelope genes LCE3A and LCE3B and a number of")
    print("     uncharacterized LOC genes, with one thing in common:")
    print("     they are single exon genes, with no introns to dilute")
    print("     the coding part.")
    print(f"     {len(impossible)} gene comes out above 100%, which is")
    print("     impossible and so is worth looking at:")
    for _, row in impossible.iterrows():
        where = row["chromosome"]
        print(f"       {row['symbol']} on chromosome {where}: "
              f"{row['length']:.0f} bases")
        print(f"            but a protein of {row['protein_length']:.0f} "
              "amino acids needs 168 with its stop codon")
    print("     The formula assumes the gene record spans the whole")
    print("     coding sequence. For this one annotation it does not, by")
    print("     three bases. One row in twenty thousand, and the only")
    print("     reason you see it is that the number went somewhere a")
    print("     number cannot go.")


def exercise_5_5(table):
    hb = table[table["symbol"].astype(str).str.startswith("HB")]
    hb = hb[hb["kind"] == "protein_coding"]
    print(f"\n5.5  {len(hb)} protein coding genes have a symbol starting")
    print("     HB, on these chromosomes:")
    for chromosome, group in hb.groupby("chromosome"):
        symbols = ", ".join(sorted(group["symbol"].astype(str)))
        span = group["end"].max() - group["start"].min()
        print(f"       chromosome {chromosome:2}: {symbols}")
        print(f"            spread over {span / 1000:,.0f} kb")
    print("     The beta like genes sit together in a cluster on")
    print("     chromosome 11 and the alpha like genes on chromosome 16.")
    print("     Each cluster holds copies of one ancestral gene,")
    print("     which is what Chapter 13 builds a tree from. Not every")
    print("     HB symbol is a globin, though: check before assuming.")


def exercise_5_6(table):
    counts = table["kind"].value_counts().head(10).sort_values()
    figure, axes = plt.subplots(figsize=(6, 4))
    axes.barh(counts.index, counts.values, color="#555555")
    axes.set_xlabel("Genes")
    figure.savefig("ch05-kinds-bar.png", dpi=110, bbox_inches="tight")
    plt.close(figure)
    print("\n5.6  The ten commonest kinds, most to least:")
    for kind, count in counts.sort_values(ascending=False).items():
        print(f"       {kind:22} {count:7,}")
    print("     barh puts the names along the vertical axis, where long")
    print("     names fit. Saved as ch05-kinds-bar.png.")


if __name__ == "__main__":
    data = genes()
    exercise_5_1(data)
    exercise_5_2(data)
    exercise_5_3(data)
    exercise_5_4(data)
    exercise_5_5(data)
    exercise_5_6(data)
