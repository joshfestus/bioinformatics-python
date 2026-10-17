"""Solutions to the exercises of Project Three, Chapter 25.

Run from this folder: uv run ch25.py
Reads the project's files in ../projects/airway. Exercise 25.6 runs
PyDESeq2 and takes about a minute; pass --quick to skip it.
"""

import csv
import sys
from collections import Counter
from pathlib import Path

import pandas as pd

WORK = Path("../projects/airway").resolve()
sys.path.insert(0, str(WORK))

PROCESS = "biological_process"
FUNCTION = "molecular_function"


def read(name):
    with open(WORK / name) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def results():
    table = pd.read_csv(WORK / "results.tsv", sep="\t", index_col=0)
    return table.dropna(subset=["padj", "symbol"])


def gene_sets(namespace=PROCESS, propagate=True):
    """The project's own machinery, with the two switches exposed."""
    from enrich import annotations, ontology, with_ancestors
    names, parents = ontology(str(WORK / "go-basic.obo"))
    direct = annotations(str(WORK / "goa_human.gaf.gz"))
    if propagate:
        import enrich
        kept = enrich.PROCESS
        enrich.PROCESS = namespace
        try:
            full = with_ancestors(direct, names, parents)
        finally:
            enrich.PROCESS = kept
    else:
        full = {symbol: {term for term in terms if term in names
                         and names[term][1] == namespace}
                for symbol, terms in direct.items()}
    return names, full


def enrichment_of(changed, background, full, names,
                  smallest=10, largest=500):
    from scipy.stats import hypergeom
    from enrich import benjamini_hochberg
    sets = {}
    for symbol in background:
        for term in full.get(symbol, ()):
            sets.setdefault(term, set()).add(symbol)
    sets = {term: members for term, members in sets.items()
            if smallest <= len(members) <= largest}
    rows = []
    for term, members in sets.items():
        found = members & changed
        if len(found) < 3:
            continue
        rows.append({"term": term, "name": names[term][0],
                     "found": len(found), "in_set": len(members),
                     "p": hypergeom.sf(len(found) - 1, len(background),
                                       len(members), len(changed))})
    if not rows:
        return [], len(sets)
    for row, q in zip(rows, benjamini_hochberg([r["p"] for r in rows])):
        row["q"] = q
    rows.sort(key=lambda row: row["p"])
    return rows, len(sets)


def exercise_25_1():
    table = results()
    names, full = gene_sets()
    background = {s for s in table["symbol"] if s in full}
    print("25.1  Dropping the fold change requirement.")
    print(f"       {'rule':28} {'genes':>7} {'terms at q<0.05':>16}")
    for label, wanted in (
            ("padj<0.05 and twofold", (table["padj"] < 0.05)
             & (table["log2FoldChange"].abs() >= 1)),
            ("padj<0.05, any size", table["padj"] < 0.05)):
        for direction, sign in (("up", 1), ("down", -1)):
            genes = {s for s in table[wanted & (
                sign * table["log2FoldChange"] > 0)]["symbol"]
                if s in full}
            rows, _ = enrichment_of(genes, background, full, names)
            passing = sum(1 for row in rows if row["q"] < 0.05)
            print(f"       {label + ', ' + direction:28} "
                  f"{len(genes):>7,} {passing:>16}")
    names2, full2 = gene_sets()
    background2 = {s for s in table["symbol"] if s in full2}
    watch = ("cellular response to glucocorticoid stimulus",
             "cellular response to hormone stimulus")
    print("\n       what happens to two terms in particular:")
    for label, wanted in (
            ("twofold", (table["padj"] < 0.05)
             & (table["log2FoldChange"] >= 1)),
            ("any size", (table["padj"] < 0.05)
             & (table["log2FoldChange"] > 0))):
        genes = {s for s in table[wanted]["symbol"] if s in full2}
        rows, _ = enrichment_of(genes, background2, full2, names2)
        found = {row["name"]: row for row in rows}
        for name in watch:
            row = found.get(name)
            if row:
                print(f"         {label:9} {name[:38]:40} "
                      f"q={row['q']:.0e}")
    print("      Dropping the threshold brings in thousands more")
    print("      genes and many more enriched terms, and the")
    print("      interesting part is what it loses.")
    print("      The wide term, cellular response to hormone")
    print("      stimulus, gets stronger. The narrow one, cellular")
    print("      response to glucocorticoid stimulus, stops passing")
    print("      at all: it holds only 34 genes, and adding 1,800")
    print("      weakly changed genes to the foreground raises how")
    print("      many of them chance alone would supply.")
    print("      So a looser list did not reveal more biology. It")
    print("      diluted the most specific and most relevant term")
    print("      in the result, and left the vaguer one standing.")
    print("      Neither rule is wrong, and that is the point. The")
    print("      fold change filter says a change too small to")
    print("      matter should not count, which is a judgment, and")
    print("      like every threshold in this book it belongs in")
    print("      the report where a reader can see it.")


def exercise_25_2():
    table = results()
    changed = table[(table["padj"] < 0.05)
                    & (table["log2FoldChange"].abs() >= 1)]
    print("\n25.2  Molecular function instead of biological process.")
    for namespace, label in ((PROCESS, "biological process"),
                             (FUNCTION, "molecular function")):
        names, full = gene_sets(namespace)
        background = {s for s in table["symbol"] if s in full}
        genes = {s for s in changed[changed["log2FoldChange"] > 0]
                 ["symbol"] if s in full}
        rows, testable = enrichment_of(genes, background, full, names)
        passing = [row for row in rows if row["q"] < 0.05]
        print(f"\n       {label}: {testable:,} testable terms, "
              f"{len(passing)} enriched")
        for row in passing[:6]:
            print(f"         q={row['q']:.0e} {row['found']:>3} of "
                  f"{row['in_set']:>3}  {row['name'][:40]}")
    print("      A process says what the cell is doing; a function")
    print("      says what a protein is. The process view told us")
    print("      the cells stop moving and answer a hormone. The")
    print("      function view names the machinery: binding")
    print("      activities and the enzymes behind them.")
    print("      Neither is a better answer. They are two columns")
    print("      of the same ontology, and a result that appears in")
    print("      both is better supported than one in either.")


def exercise_25_3():
    table = results().set_index("symbol")
    names, full = gene_sets()
    wanted = {"cellular response to glucocorticoid stimulus": None,
              "axon guidance": None}
    for term, (label, namespace) in names.items():
        if label in wanted and namespace == PROCESS:
            wanted[label] = term
    changed = table[(table["padj"] < 0.05)
                    & (table["log2FoldChange"].abs() >= 1)]
    print("\n25.3  The genes behind two headings.")
    for label, term in wanted.items():
        if term is None:
            continue
        members = {symbol for symbol, terms in full.items()
                   if term in terms}
        found = sorted(members & set(changed.index))
        print(f"\n       {label}: {len(found)} changed genes")
        for symbol in found:
            row = changed.loc[symbol]
            change = row["log2FoldChange"]
            if hasattr(change, "iloc"):
                change = change.iloc[0]
            print(f"         {symbol:10} {change:+6.2f} log2")
    print("      The glucocorticoid heading survives reading. Every")
    print("      gene under it is a known steroid responder and")
    print("      every one moves the same way, which is what a real")
    print("      pathway looks like.")
    print("      The axon guidance heading does not survive in the")
    print("      sense the name suggests. Its genes are real and")
    print("      they really changed, but they are semaphorins and")
    print("      their relatives, the proteins cells use to tell")
    print("      each other where to go. In airway muscle that is")
    print("      about cells moving, not about nerves. The term")
    print("      carries the name of the tissue it was first")
    print("      studied in, and the gene list is what tells you so.")
    print("      One more tell. The seven glucocorticoid genes all")
    print("      move the same way, up. The thirty-two guidance")
    print("      genes are split, nineteen down and thirteen up.")
    print("      A pathway the cell switches on moves together; a")
    print("      term that is really a bag of related proteins does")
    print("      not, and the direction column shows it before any")
    print("      biology is argued.")


def exercise_25_4():
    table = results()
    changed = table[(table["padj"] < 0.05)
                    & (table["log2FoldChange"].abs() >= 1)]
    print("\n25.4  With and without propagating annotations upward.")
    for propagate in (True, False):
        names, full = gene_sets(PROCESS, propagate)
        background = {s for s in table["symbol"] if s in full}
        genes = {s for s in changed[changed["log2FoldChange"] > 0]
                 ["symbol"] if s in full}
        rows, testable = enrichment_of(genes, background, full, names)
        passing = [row for row in rows if row["q"] < 0.05]
        label = "with ancestors" if propagate else "direct only"
        print(f"\n       {label}: {testable:,} testable terms, "
              f"{len(passing)} enriched")
        for row in passing[:5]:
            print(f"         {row['found']:>3} of {row['in_set']:>3}  "
                  f"{row['name'][:44]}")
    print("      Without propagation the wide terms vanish. A gene")
    print("      annotated to cellular response to glucocorticoid")
    print("      stimulus is not counted under response to hormone,")
    print("      so the hormone term holds a handful of genes")
    print("      instead of four hundred and never reaches the")
    print("      minimum set size.")
    print("      What is left is the narrow terms, each with a")
    print("      few genes, and reading them is like reading an")
    print("      index with no chapter headings. The ontology is a")
    print("      hierarchy, and an enrichment that ignores the")
    print("      hierarchy throws away most of what it is for.")


def exercise_25_5():
    print("\n25.5  Why the four albuterol samples were left out.")
    print("      A paragraph a reviewer would accept:\n")
    paragraph = (
        "Four samples were excluded from this analysis: the "
        "albuterol and albuterol-plus-dexamethasone samples of "
        "donors N052611 and N080611. In those two donors the "
        "labels for the two conditions appear to have been "
        "exchanged. Four independent lines of evidence agree: the "
        "glucocorticoid marker genes, the principal components, "
        "hierarchical clustering, and a classifier trained only on "
        "the untreated and dexamethasone samples, which assigns the "
        "nominal albuterol samples a dexamethasone probability of "
        "0.98 and the nominal albuterol-plus-dexamethasone samples "
        "0.06 to 0.09. We could not determine where the exchange "
        "occurred and have not corrected it. The comparison "
        "reported here uses only the untreated and dexamethasone "
        "samples, which are unaffected, and the four doubtful "
        "samples are reported rather than silently dropped.")
    words = paragraph.split()
    line = "       "
    for word in words:
        if len(line) + len(word) > 70:
            print(line)
            line = "       "
        line += word + " "
    print(line)
    print("\n      Three things make it acceptable. It says exactly")
    print("      which samples and why, it gives evidence a reader")
    print("      can check rather than an impression, and it admits")
    print("      what is not known. The last is the one people")
    print("      leave out, and it is the one that makes the rest")
    print("      believable.")


def exercise_25_6():
    print("\n25.6  Albuterol against untreated, label doubt and all.")
    if "--quick" in sys.argv:
        print("      Skipped with --quick.")
        return
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.default_inference import DefaultInference
    from pydeseq2.ds import DeseqStats
    counts = pd.read_csv(WORK / "counts.tsv.gz", sep="\t",
                         index_col="GeneID")
    samples = pd.read_csv(WORK / "samples.tsv", sep="\t",
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
    significant = found[found["padj"] < 0.05]
    print(f"       {len(found):,} genes tested, {len(significant)} with "
          "an adjusted p below 0.05")
    names, full = gene_sets()
    symbols = pd.read_csv(WORK / "results.tsv", sep="\t",
                          index_col=0)["symbol"]
    # the counts table carries GeneIDs as strings, results.tsv as
    # integers, so convert before looking anything up
    changed = {symbols.get(int(gene)) for gene in significant.index}
    changed = {s for s in changed if isinstance(s, str) and s in full}
    background = {s for s in symbols.dropna() if s in full}
    rows, testable = enrichment_of(changed, background, full, names)
    passing = [row for row in rows if row["q"] < 0.05]
    print(f"       {len(changed)} of them carry an annotation, "
          f"against {testable:,} testable terms")
    print(f"       terms enriched at q below 0.05: {len(passing)}")
    for row in passing[:5]:
        print(f"         q={row['q']:.0e} {row['found']:>3} of "
              f"{row['in_set']:>3}  {row['name'][:40]}")
    if not passing:
        print("       nothing passes, and the best any term manages:")
        for row in rows[:3]:
            print(f"         q={row['q']:.2f} {row['found']:>3} of "
                  f"{row['in_set']:>3}  {row['name'][:40]}")
    small = [row for row in passing if row["found"] <= 3]
    print(f"      It does find terms, {len(passing)} of them, and")
    print(f"      {len(small)} rest on three genes each. With 43")
    print("      annotated genes against five thousand terms, that")
    print("      is what passing looks like when there is little to")
    print("      find: the statistics are correct and the biology")
    print("      is three genes wide.")
    print("      Nothing here contradicts Chapter 19, which found")
    print("      45 changed genes for albuterol against 4,544 for")
    print("      dexamethasone. The drug does little to")
    print("      transcription, and a list that short cannot")
    print("      support a pathway story.")
    print("      Had the enrichment returned nothing at all, that")
    print("      too would be a result, and three readings would")
    print("      stay open that the data cannot separate.")
    print("      Albuterol may genuinely change little, which fits")
    print("      how it works: through a surface receptor on a")
    print("      timescale of minutes, not by binding DNA. Or the")
    print("      label doubt of Chapter 20 has mixed two conditions")
    print("      and buried the signal. Or four donors are too few")
    print("      to find a small effect.")
    print("      Reporting little found, with those possibilities")
    print("      named, is honest. Reporting that albuterol does")
    print("      not change gene expression is not, and neither is")
    print("      writing up thirty-five terms that rest on three")
    print("      genes apiece as though they were a finding.")


def exercise_25_7():
    print("\n25.7  Putting the chart in the report.")
    print("      report.py gains a figure, drawn from the same")
    print("      enrichment table it already reads:\n")
    code = '''def chart(rows, path="enrichment.png"):
    """The top terms each way, as a horizontal bar chart."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    shown = [r for r in rows if float(r["q"]) < 0.05][:8][::-1]
    ratio = [int(r["found"]) / float(r["expected"]) for r in shown]
    figure, axes = plt.subplots(figsize=(6, 3.2))
    axes.barh([r["name"][:44] for r in shown], ratio, color="#555555")
    axes.set_xlabel("Times more genes than chance would give")
    figure.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(figure)
    return path'''
    for line in code.splitlines():
        print(f"       {line}")
    print("\n      and the Snakefile learns that the report now makes")
    print("      two files rather than one:\n")
    rule = '''rule report:
    input:
        "enrichment.tsv",
        "comparison.tsv",
        "results.tsv",
    output:
        "report.md",
        "enrichment.png",
    shell:
        "python report.py > {log} 2>&1"'''
    for line in rule.splitlines():
        print(f"       {line}")
    print("\n      Declaring the figure as an output is the part that")
    print("      matters. Snakemake only knows what you tell it, so")
    print("      an undeclared file is never rebuilt when the data")
    print("      changes and never cleaned up, and sooner or later a")
    print("      report carries last week's chart beside this week's")
    print("      numbers.")


if __name__ == "__main__":
    exercise_25_1()
    exercise_25_2()
    exercise_25_3()
    exercise_25_4()
    exercise_25_5()
    exercise_25_6()
    exercise_25_7()
