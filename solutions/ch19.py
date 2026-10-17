"""Solutions to the exercises of Chapter 19.

Run from this folder: uv run ch19.py
Exercise 19.4 runs PyDESeq2 again and takes about a minute; pass
--quick to skip it.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control

CH05, CH18, CH19, DATA = (Path("../ch05"), Path("../ch18"),
                          Path("../ch19"), Path("../data"))


def results():
    table = pd.read_csv(CH19 / "results.tsv", sep="\t", index_col=0)
    return table.dropna(subset=["padj"])


def exercise_19_1():
    table = results()
    changed = table[table["padj"] < 0.05]
    print("19.1  The ten genes with the largest fall, as asked:")
    print(f"       {'gene':14} {'log2':>7} {'fold':>7} {'baseMean':>10}")
    for gene_id, row in changed.nsmallest(10,
                                          "log2FoldChange").iterrows():
        print(f"       {str(row['symbol']):14} "
              f"{row['log2FoldChange']:>7.2f} "
              f"{2 ** abs(row['log2FoldChange']):>6.0f}x "
              f"{row['baseMean']:>10,.0f}")
    print("      Every one of them has a baseMean under thirty, and")
    print("      more than half have no name beyond a LOC or LINC")
    print("      identifier. That is what sorting by fold change")
    print("      alone does: a gene averaging eight reads needs very")
    print("      little movement to show a fiftyfold ratio, and the")
    print("      top of the list fills with genes nobody can follow")
    print("      up. They passed the test, so they are probably")
    print("      real; they are just not where the biology is.")
    print("\n      The same question of genes with real expression")
    print("      behind them, a baseMean of 1,000 or more:\n")
    jobs = {"CXCL12": "a chemokine, which calls immune cells in",
            "BDKRB2": "the bradykinin receptor, inflammation",
            "SOX4": "a transcription factor",
            "SCD": "makes unsaturated fatty acids",
            "FST": "follistatin, a growth factor binder",
            "KCTD12": "a potassium channel regulator",
            "ADAMTS14": "a protease that trims collagen",
            "LRRC17": "secreted, acts on bone turnover",
            "CCKAR": "a hormone receptor",
            "MEST": "an imprinted gene of development"}
    solid = changed[changed["baseMean"] >= 1000]
    print(f"       {'gene':14} {'log2':>7} {'baseMean':>10}  what")
    for gene_id, row in solid.nsmallest(10, "log2FoldChange").iterrows():
        symbol = str(row["symbol"])
        print(f"       {symbol:14} {row['log2FoldChange']:>7.2f} "
              f"{row['baseMean']:>10,.0f}  {jobs.get(symbol, '')}")
    print("      This list can be read. CXCL12 is a chemokine, one")
    print("      of the signals that summons immune cells, and it is")
    print("      running at 36,000 counts before the drug and a")
    print("      quarter of that after. BDKRB2 receives bradykinin,")
    print("      which inflames and constricts the airway. A steroid")
    print("      turning both down is the whole reason it is given")
    print("      for asthma.")
    print("      The lesson is not that fold change is wrong. It is")
    print("      that fold change and expression answer different")
    print("      questions, and a list sorted by one should be read")
    print("      with the other in view.")


def benjamini_hochberg(pvalues):
    """The step-up procedure, written out as Chapter 19 describes it."""
    order = np.argsort(pvalues)
    ranked = np.asarray(pvalues)[order]
    total = len(pvalues)
    adjusted = ranked * total / np.arange(1, total + 1)
    # walk back from the largest, keeping each value no bigger than the
    # one after it, so the adjustment never falls as the p-value rises
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    out = np.empty(total)
    out[order] = np.clip(adjusted, 0, 1)
    return out


def exercise_19_2():
    table = results()
    mine = benjamini_hochberg(table["pvalue"].to_numpy())
    theirs = false_discovery_control(table["pvalue"].to_numpy())
    difference = np.abs(mine - theirs).max()
    print(f"\n19.2  {len(table):,} p-values adjusted two ways.")
    print(f"       largest disagreement: {difference:.2e}")
    print(f"       mine says {int((mine < 0.05).sum()):,} genes pass "
          f"0.05, SciPy says {int((theirs < 0.05).sum()):,}")
    print(f"       against the padj column: "
          f"{int((table['padj'] < 0.05).sum()):,}")
    print("      The whole procedure is three lines: multiply each")
    print("      p-value by the number of tests and divide by its")
    print("      rank, then walk back from the largest making sure")
    print("      the list never rises as you go down it.")
    print("      That last step is the one people leave out. Without")
    print("      it a gene can end up with a smaller adjusted value")
    print("      than a gene with a smaller p-value, which is")
    print("      nonsense, and it is why the procedure is called")
    print("      step-up rather than just an arithmetic correction.")


def exercise_19_3():
    table = results()
    strict = table[(table["padj"] < 0.01)
                   & (table["log2FoldChange"].abs() >= 2)]
    up = (strict["log2FoldChange"] > 0).sum()
    print(f"\n19.3  {len(strict):,} genes pass an adjusted p below 0.01")
    print(f"       and a fourfold change: {up} up and "
          f"{len(strict) - up} down.")
    print("       The three with the most reads behind them:")
    for gene_id, row in strict.nlargest(3, "baseMean").iterrows():
        print(f"         {str(row['symbol']):10} baseMean "
              f"{row['baseMean']:>9,.0f}, "
              f"{2 ** row['log2FoldChange']:>5.1f} fold, padj "
              f"{row['padj']:.1e}")
    print("      baseMean is the average normalized count, so these")
    print("      are the genes that are both highly expressed and")
    print("      strongly changed. They are the ones worth a")
    print("      laboratory's time: a gene that doubles from two")
    print("      reads to four may be real and will be hard to")
    print("      confirm at the bench.")


def exercise_19_4():
    print("\n19.4  Albuterol against untreated, all four donors:")
    if "--quick" in sys.argv:
        print("      Skipped with --quick.")
        return
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.default_inference import DefaultInference
    from pydeseq2.ds import DeseqStats
    counts = pd.read_csv(DATA / "GSE52778_counts.tsv.gz", sep="\t",
                         index_col="GeneID")
    samples = pd.read_csv(CH18 / "samples.tsv", sep="\t",
                          index_col="sample")
    samples = samples[samples["treatment"].isin(["Untreated",
                                                 "Albuterol"])]
    table = counts[samples.index].T
    table = table.loc[:, table.sum() >= 10]
    inference = DefaultInference(n_cpus=2)
    dataset = DeseqDataSet(counts=table, metadata=samples,
                           design="~cell_line + treatment",
                           inference=inference, quiet=True)
    dataset.deseq2()
    stats = DeseqStats(dataset,
                       contrast=["treatment", "Albuterol", "Untreated"],
                       inference=inference, quiet=True)
    stats.summary()
    found = stats.results_df.dropna(subset=["padj"])
    significant = (found["padj"] < 0.05).sum()
    twofold = ((found["padj"] < 0.05)
               & (found["log2FoldChange"].abs() >= 1)).sum()
    dex = results()
    print(f"       albuterol:     {significant:>6,} genes at padj<0.05, "
          f"{twofold:>5,} of them twofold")
    dex_significant = (dex["padj"] < 0.05).sum()
    dex_twofold = ((dex["padj"] < 0.05)
                   & (dex["log2FoldChange"].abs() >= 1)).sum()
    print(f"       dexamethasone: {dex_significant:>6,} genes at "
          f"padj<0.05, {dex_twofold:>5,} of them twofold")
    print("      Albuterol changes far less. It relaxes the muscle")
    print("      through a receptor on the cell surface, and that is")
    print("      mostly a fast signal rather than a change in which")
    print("      genes are read. A steroid works by walking into the")
    print("      nucleus and binding DNA, so changing transcription")
    print("      is the whole of what it does.")
    print("      Remember Chapter 20's warning: two donors' albuterol")
    print("      labels look swapped with their albuterol plus")
    print("      dexamethasone ones. A comparison using those samples")
    print("      is built on labels this book does not fully trust,")
    print("      and the honest thing is to say so beside the number.")


def exercise_19_5():
    table = results()
    tests = len(table)
    bonferroni = np.clip(table["pvalue"] * tests, 0, 1)
    fdr = table["padj"]
    both = ((bonferroni < 0.05) & (fdr < 0.05)).sum()
    print(f"\n19.5  {tests:,} tests, so Bonferroni multiplies every")
    print("      p-value by that number.")
    strict = int((bonferroni < 0.05).sum())
    loose = int((fdr < 0.05).sum())
    print(f"       Bonferroni below 0.05: {strict:,} genes")
    print(f"       false discovery rate:  {loose:,} genes")
    print(f"       in both lists:         {int(both):,}")
    print("      Bonferroni asks for a near certainty that not one")
    print("      gene on the list is a false positive. The false")
    print("      discovery rate asks that at most one in twenty is,")
    print("      which for a screen is the more useful promise: you")
    print("      are going to check the hits at the bench anyway.")
    print("      Every Bonferroni gene is in the FDR list, because")
    print("      the stricter test cannot pass something the looser")
    print("      one rejects. Which to use is a question about what")
    print("      the list is for, not about which is correct.")


if __name__ == "__main__":
    exercise_19_1()
    exercise_19_2()
    exercise_19_3()
    exercise_19_4()
    exercise_19_5()
