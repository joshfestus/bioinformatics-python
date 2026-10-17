# Three views of the human gene table made by gene_table.py.
import matplotlib

matplotlib.use("svg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# The book's style: gray scale, no top or right frame, small type
plt.rcParams.update(
    {
        "font.size": 9,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.edgecolor": "#444444",
        "xtick.color": "#444444",
        "ytick.color": "#444444",
    }
)
INK, MID, LIGHT = "#222222", "#777777", "#BBBBBB"
CHROMOSOMES = [str(n) for n in range(1, 23)] + ["X", "Y"]

genes = pd.read_csv("genes.tsv", sep="\t", dtype={"chromosome": str})
coding = genes[genes["kind"] == "protein_coding"]


def length_histogram(folder: str = ".") -> None:
    fig, ax = plt.subplots(figsize=(5.4, 2.8))
    bins = np.logspace(2, 7, 61)
    ax.hist(coding["length"], bins=bins, color=MID)
    median = coding["length"].median()
    ax.axvline(median, color=INK, linestyle="--", linewidth=0.9)
    ax.text(
        median / 1.15,
        ax.get_ylim()[1] * 0.95,
        f"median {median:,.0f} bp",
        ha="right",
        color=INK,
        fontsize=8,
    )
    ax.set_xscale("log")
    ax.set_xlabel("Gene length on the genome (base pairs, log scale)")
    ax.set_ylabel("Number of genes")
    fig.savefig(f"{folder}/gene-lengths.svg", bbox_inches="tight")


def per_chromosome(folder: str = ".") -> None:
    counts = coding["chromosome"].value_counts().reindex(CHROMOSOMES)
    fig, ax = plt.subplots(figsize=(5.4, 2.6))
    ax.bar(counts.index, counts.to_numpy(), color=MID, width=0.7)
    ax.set_xlabel("Chromosome")
    ax.set_ylabel("Protein coding genes")
    ax.tick_params(axis="x", labelsize=7.5)
    fig.savefig(f"{folder}/genes-per-chromosome.svg", bbox_inches="tight")


def protein_against_gene(folder: str = ".") -> None:
    fig, ax = plt.subplots(figsize=(5.0, 3.6))
    ax.scatter(
        coding["length"],
        coding["protein_length"],
        s=2,
        color=LIGHT,
        linewidths=0,
    )
    for symbol in ("HBB", "TTN", "RBFOX1"):
        row = coding[coding["symbol"] == symbol].iloc[0]
        ax.scatter(row["length"], row["protein_length"], s=14, color=INK)
        ax.annotate(
            symbol,
            (row["length"], row["protein_length"]),
            xytext=(5, -10),
            textcoords="offset points",
            fontsize=8,
            color=INK,
        )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Gene length on the genome (base pairs)")
    ax.set_ylabel("Longest protein (amino acids)")
    fig.savefig(f"{folder}/protein-vs-gene.svg", bbox_inches="tight")


if __name__ == "__main__":
    length_histogram()
    per_chromosome()
    protein_against_gene()
    print("saved gene-lengths.svg, genes-per-chromosome.svg,")
    print("      protein-vs-gene.svg")
