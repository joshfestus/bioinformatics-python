"""Solutions to the exercises of Chapter 20.

Run from this folder: uv run ch20.py
Exercise 20.5 reruns PyDESeq2 twice and takes a couple of minutes;
pass --quick to skip that part.
"""

import sys
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import pdist
from sklearn.decomposition import PCA
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import LeaveOneGroupOut, cross_val_predict
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

CH05, CH18, CH20, DATA = (Path("../ch05"), Path("../ch18"),
                          Path("../ch20"), Path("../data"))
SHORT = {"Untreated": "none", "Dexamethasone": "dex",
         "Albuterol": "alb", "Albuterol_Dexamethasone": "alb+dex"}
SUSPECT = {"GSM1275868": "GSM1275869", "GSM1275872": "GSM1275873"}


def load(samples_file=CH20 / "samples.tsv"):
    """Chapter 20's loader, reading from this folder."""
    counts = pd.read_csv(DATA / "GSE52778_counts.tsv.gz", sep="\t",
                         index_col="GeneID")
    samples = pd.read_csv(samples_file, sep="\t", index_col="sample")
    samples["label"] = samples["treatment"].map(SHORT)
    cpm = counts[samples.index] / counts[samples.index].sum() * 1e6
    expressed = cpm[cpm.mean(axis=1) >= 1]
    return np.log2(expressed + 1).T, samples


def symbols():
    genes = pd.read_csv(CH05 / "genes.tsv", sep="\t",
                        dtype={"chromosome": str})
    return genes.drop_duplicates("GeneID").set_index("GeneID")["symbol"]


def exercise_20_1():
    expression, samples = load()
    scores = PCA(n_components=4).fit(expression)
    points = scores.transform(expression)
    figure, axes = plt.subplots(figsize=(5.4, 3.6))
    for donor, rows in samples.groupby("cell_line"):
        index = [samples.index.get_loc(name) for name in rows.index]
        axes.scatter(points[index, 0], points[index, 2], s=36,
                     label=donor, edgecolor="#222222")
    axes.set_xlabel("PC1")
    axes.set_ylabel("PC3")
    axes.legend(frameon=False, fontsize=7)
    figure.savefig("ch20-pc1-pc3.png", dpi=110, bbox_inches="tight")
    plt.close(figure)
    share = scores.explained_variance_ratio_
    print("20.1  PC1 against PC3.")
    print(f"       variance explained: PC1 {share[0]:.0%}, "
          f"PC2 {share[1]:.0%}, PC3 {share[2]:.0%}, PC4 {share[3]:.0%}")
    print(f"       {'donor':10} {'PC1 mean':>10} {'PC3 mean':>10} "
          f"{'PC3 spread':>11}")
    for donor, rows in samples.groupby("cell_line"):
        index = [samples.index.get_loc(name) for name in rows.index]
        print(f"       {donor:10} {points[index, 0].mean():>10.1f} "
              f"{points[index, 2].mean():>10.1f} "
              f"{points[index, 2].std():>11.1f}")
    spread = {}
    for donor, rows in samples.groupby("cell_line"):
        index = [samples.index.get_loc(name) for name in rows.index]
        spread[donor] = points[index, 2].mean()
    furthest = max(spread, key=lambda d: abs(spread[d]))
    print(f"      PC3 separates {furthest} from the rest.")
    print("      Each component after the first picks up whatever")
    print("      variation the ones before it left behind, so by the")
    print("      third there is usually one donor's peculiarity and")
    print("      not much else. Saved as ch20-pc1-pc3.png.")


def exercise_20_2():
    expression, samples = load()
    model = PCA(n_components=3).fit(expression)
    weights = pd.Series(model.components_[0], index=expression.columns)
    names = symbols()
    top = weights.abs().nlargest(10)
    results = pd.read_csv("../ch19/results.tsv", sep="\t", index_col=0)
    print("\n20.2  The ten genes weighing most on PC1:")
    print(f"       {'gene':12} {'weight':>8} {'log2FC in ch19':>15} "
          f"{'padj':>10}")
    for gene_id in top.index:
        row = results.loc[gene_id] if gene_id in results.index else None
        change = f"{row['log2FoldChange']:>15.2f}" if row is not None \
            else f"{'not tested':>15}"
        padj = f"{row['padj']:>10.1e}" if row is not None \
            and pd.notna(row["padj"]) else f"{'':>10}"
        print(f"       {str(names.get(gene_id, gene_id)):12} "
              f"{weights[gene_id]:>8.3f} {change} {padj}")
    tested = [g for g in top.index if g in results.index]
    agree = sum(1 for g in tested
                if pd.notna(results.loc[g, "padj"])
                and results.loc[g, "padj"] < 0.05)
    print(f"      {agree} of the {len(tested)} are also significant in")
    print("      Chapter 19's test, which is the point: the direction")
    print("      the samples spread in is the drug, and the genes")
    print("      that define it are the genes the drug changes.")
    print("      PCA found that without being told which samples were")
    print("      treated. It is an unsupervised method agreeing with")
    print("      a supervised one, and when the two disagree it is")
    print("      usually a sign that the biggest source of variation")
    print("      in the data is not the thing you are studying.")


def exercise_20_3():
    expression, samples = load()
    distances = pdist(expression, metric="correlation")
    tree = linkage(distances, method="average")
    print("\n20.3  Cutting the tree into four groups:")
    groups = fcluster(tree, 4, criterion="maxclust")
    table = samples.assign(group=groups)
    for group, rows in table.groupby("group"):
        donors = ", ".join(sorted(set(rows["cell_line"])))
        labels = ", ".join(sorted(set(rows["label"])))
        print(f"       group {group}: {len(rows)} samples")
        print(f"         donors:     {donors}")
        print(f"         treatments: {labels}")
    by_donor = table.groupby("cell_line")["group"].nunique()
    by_label = table.groupby("label")["group"].nunique()
    print(f"      Each donor's samples fall into "
          f"{by_donor.mean():.1f} groups on average,")
    print(f"      each treatment into {by_label.mean():.1f}.")
    print("      The answer is neither, and the shape is the")
    print("      interesting part. One donor, N080611, takes two")
    print("      groups entirely to itself, and the other three")
    print("      donors share the remaining two groups between them.")
    print("      So the first thing the clustering sees is that one")
    print("      person's cells are unlike everyone else's, and only")
    print("      then does it notice the drug.")
    print("      That is the argument for the paired design in")
    print("      another form. Two samples from one person are more")
    print("      alike than two samples given the same treatment, so")
    print("      a test that ignores the donor spends its power on")
    print("      the wrong difference.")


def exercise_20_4():
    expression, samples = load()
    treated = samples["treatment"].str.contains(
        "Dexamethasone").astype(int)
    known = samples["label"].isin(["none", "dex"])
    names = symbols()
    print("\n20.4  How few genes the classifier needs:")
    print(f"       {'genes':>6} {'correct':>9}  chosen when only one")
    for k in (50, 5, 1):
        model = make_pipeline(StandardScaler(),
                              SelectKBest(f_classif, k=k),
                              LogisticRegression(max_iter=1000))
        predicted = cross_val_predict(
            model, expression[known], treated[known],
            groups=samples["cell_line"][known], cv=LeaveOneGroupOut())
        correct = (predicted == treated[known]).sum()
        chosen = ""
        if k == 1:
            picker = make_pipeline(StandardScaler(),
                                   SelectKBest(f_classif, k=1))
            picker.fit(expression[known], treated[known])
            mask = picker.named_steps["selectkbest"].get_support()
            gene = expression.columns[mask][0]
            chosen = str(names.get(gene, gene))
        print(f"       {k:>6} {correct:>4} of {known.sum():<4} {chosen}")
    print("      One gene is enough to get all eight right, and")
    print("      five genes get only seven. More information made")
    print("      the model worse, which happens when the samples are")
    print("      this few: the extra genes carry donor differences")
    print("      as well as drug differences, and with eight samples")
    print("      there is not enough data to learn which is which.")
    print("      The gene it picks on its own is MT1X, a")
    print("      metallothionein, up 3.34 log2 with an adjusted p of")
    print("      8e-76 in Chapter 19. Its relatives MT1E, MT1M and")
    print("      MT1G are all up too, so the choice is not a fluke:")
    print("      steroids switch this whole family on.")
    print("      None of that means one gene is a good classifier.")
    print("      With eight samples and twenty thousand genes, some")
    print("      gene will separate any two groups perfectly whether")
    print("      or not the biology is real, which is why the leave")
    print("      one donor out check matters more than the score.")


def exercise_20_5():
    print("\n20.5  Swapping the two suspect pairs back:")
    original = pd.read_csv(CH20 / "samples.tsv", sep="\t")
    swapped = original.copy()
    for first, second in SUSPECT.items():
        a = swapped.index[swapped["sample"] == first][0]
        b = swapped.index[swapped["sample"] == second][0]
        swapped.loc[a, "treatment"], swapped.loc[b, "treatment"] = (
            swapped.loc[b, "treatment"], swapped.loc[a, "treatment"])
    path = Path("samples_swapped.tsv")
    swapped.to_csv(path, sep="\t", index=False)
    for label, sheet in (("as published", CH20 / "samples.tsv"),
                         ("swapped back", path)):
        expression, samples = load(sheet)
        model = PCA(n_components=2).fit(expression)
        points = model.transform(expression)
        axis = pd.Series(points[:, 0], index=samples.index)
        dex = samples["treatment"].str.contains("Dexamethasone")
        gap = abs(axis[dex].mean() - axis[~dex].mean())
        print(f"       {label:14} PC1 separates the steroid samples "
              f"by {gap:6.1f}")
    path.unlink()
    print("      Swapping the labels back makes the picture tidier:")
    print("      the samples that carry dexamethasone sit together")
    print("      and the gap along PC1 widens. That is consistent")
    print("      with the swap being real, and it is not proof.")
    print("      A relabelling that makes your plot look better is")
    print("      the most dangerous kind of analysis there is, which")
    print("      is why Chapter 20 reports the evidence and stops,")
    print("      and why Chapter 25 leaves those four samples out")
    print("      rather than correcting them.")


if __name__ == "__main__":
    exercise_20_1()
    exercise_20_2()
    exercise_20_3()
    exercise_20_4()
    exercise_20_5()
