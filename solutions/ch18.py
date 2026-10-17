"""Solutions to the exercises of Chapter 18.

Run from this folder: uv run ch18.py
Reads the GSE52778 counts in ../data and Chapter 5's genes.tsv.
"""

from pathlib import Path

import numpy as np
import pandas as pd

CH05, CH18, DATA = Path("../ch05"), Path("../ch18"), Path("../data")
COUNTS = DATA / "GSE52778_counts.tsv.gz"
DEX = ["Untreated", "Dexamethasone"]
FKBP5, CRISPLD2 = 2289, 83716
HOUSEKEEPING = {2597: "GAPDH", 60: "ACTB", 6175: "RPLP0"}


def counts(treatments=None):
    table = pd.read_csv(COUNTS, sep="\t", index_col="GeneID")
    samples = pd.read_csv(CH18 / "samples.tsv", sep="\t")
    if treatments:
        samples = samples[samples["treatment"].isin(treatments)]
    table = table[samples["sample"]]
    return table[table.sum(axis=1) > 0], samples


def symbols():
    genes = pd.read_csv(CH05 / "genes.tsv", sep="\t",
                        dtype={"chromosome": str})
    return genes.drop_duplicates("GeneID").set_index("GeneID")


def exercise_18_1():
    table, samples = counts(DEX)
    library = table.sum()
    print("18.1  FKBP5 by hand: its count, the library size, and the")
    print("      count divided by the library and multiplied by a")
    print("      million.\n")
    print(f"       {'sample':12} {'treatment':14} {'count':>7} "
          f"{'library':>12} {'CPM by hand':>12}")
    for _, row in samples.iterrows():
        name = row["sample"]
        count = table.loc[FKBP5, name]
        by_hand = count / library[name] * 1e6
        print(f"       {name:12} {row['treatment']:14} {count:>7,} "
              f"{library[name]:>12,} {by_hand:>12.1f}")
    program = table.loc[FKBP5] / library * 1e6
    agree = np.allclose(program.to_numpy(),
                        [table.loc[FKBP5, s] / library[s] * 1e6
                         for s in samples["sample"]])
    print(f"      The program's numbers agree: {agree}.")
    print("      CPM is the simplest correction there is, and doing")
    print("      it by hand once shows what it does and does not fix.")
    print("      It makes libraries of different sizes comparable. It")
    print("      does nothing about a handful of very high genes")
    print("      taking up most of the reads, which is what the size")
    print("      factors in the rest of the chapter are for.")


def exercise_18_2():
    table, samples = counts()
    library = table.sum()
    cpm = table.loc[FKBP5] / library * 1e6
    print("\n18.2  FKBP5 in all sixteen samples:")
    print(f"       {'donor':10} {'treatment':24} {'CPM':>9}")
    joined = samples.assign(cpm=cpm[samples["sample"]].to_numpy())
    ordered = joined.sort_values(["cell_line", "treatment"])
    for _, row in ordered.iterrows():
        print(f"       {row['cell_line']:10} {row['treatment']:24} "
              f"{row['cpm']:>9.1f}")
    by_treatment = joined.groupby("treatment")["cpm"].mean()
    print("\n       mean by treatment:")
    for treatment, value in by_treatment.sort_values().items():
        print(f"         {treatment:24} {value:8.1f}")
    print("      Albuterol alone does not switch FKBP5 on. Its mean")
    print("      sits with the untreated samples, while everything")
    print("      with dexamethasone in it is many times higher.")
    print("      FKBP5 is a glucocorticoid responsive gene, and")
    print("      albuterol is not a glucocorticoid: it opens the")
    print("      airways through a different receptor. The gene")
    print("      answers the steroid and ignores the other drug,")
    print("      which is what a marker should do.")
    print("      Chapter 20 found reason to doubt the albuterol")
    print("      labels in two donors, so read these four rows with")
    print("      that in mind.")


def exercise_18_3():
    table, samples = counts(DEX)
    logged = np.log(table.replace(0, np.nan))
    reference = logged.mean(axis=1)
    ratios = logged.sub(reference, axis=0)
    all_genes = np.exp(ratios.median())
    three = np.exp(ratios.loc[list(HOUSEKEEPING)].median())
    print("\n18.3  Size factors from every gene, and from three:")
    print(f"       {'sample':12} {'all genes':>10} {'three':>8} "
          f"{'difference':>11}")
    worst, gap = None, 0
    for name in table.columns:
        difference = three[name] - all_genes[name]
        if abs(difference) > gap:
            worst, gap = name, abs(difference)
        print(f"       {name:12} {all_genes[name]:>10.3f} "
              f"{three[name]:>8.3f} {difference:>+11.3f}")
    print(f"      {worst} disagrees most, by {gap:.3f}.")
    print("      Three genes give almost the same answer as twenty")
    print("      thousand, which is the assumption behind every")
    print("      housekeeping gene normalization ever done, and here")
    print("      it holds.")
    print("      It holds because these three are high, steady and")
    print("      unaffected by the drug. Pick three genes the")
    print("      treatment happens to change and the correction")
    print("      would remove the very effect you are measuring,")
    print("      with nothing in the output to show it. The median")
    print("      over all genes needs no such luck.")


def exercise_18_4():
    table, samples = counts(DEX)
    feature = symbols()
    first = samples["sample"].iloc[0]
    lengths = feature["length"].reindex(table.index)
    kinds = feature["kind"].reindex(table.index)
    per_kb = table[first] / (lengths / 1000)
    tpm = per_kb / per_kb.sum() * 1e6
    print(f"\n18.4  The ten highest TPM genes in {first}:")
    print(f"       {'gene':12} {'TPM':>9} {'kind':16} {'bases':>8}")
    highest = tpm.sort_values(ascending=False).head(10)
    for gene_id, value in highest.items():
        print(f"       {str(feature['symbol'].get(gene_id)):12} "
              f"{value:>9,.0f} {str(kinds.get(gene_id)):16} "
              f"{lengths.get(gene_id):>8,.0f}")
    print("      Half of that list is microRNA and snoRNA genes of")
    print("      fifty to a hundred bases, and they are there because")
    print("      of the arithmetic rather than the biology. TPM")
    print("      divides by gene length, so a gene a twentieth the")
    print("      size of its neighbor needs a twentieth of the reads")
    print("      to reach the same value. On genes this short the")
    print("      division amplifies a handful of reads into a number")
    print("      that leads the table.")
    print("      Ask the question the exercise meant, of the protein")
    print("      coding genes only:\n")
    coding = tpm[kinds == "protein_coding"]
    print(f"       {'gene':12} {'TPM':>9}  what it makes")
    jobs = {"FTL": "ferritin, which stores iron",
            "EEF1A1": "elongation factor, protein synthesis",
            "FTH1": "ferritin, the other chain",
            "S100A6": "a calcium binding protein",
            "ACTB": "beta actin, the cytoskeleton",
            "ACTG1": "gamma actin, the cytoskeleton",
            "GAPDH": "glycolysis",
            "TMSB4X": "thymosin, which binds actin",
            "TMSB10": "thymosin, the other one",
            "RPLP0": "a ribosomal protein"}
    for gene_id, value in coding.sort_values(
            ascending=False).head(10).items():
        symbol = str(feature["symbol"].get(gene_id))
        print(f"       {symbol:12} {value:>9,.0f}  "
              f"{jobs.get(symbol, '')}")
    print("      Now it reads as a cell. Two actins and two thymosins")
    print("      that bind actin, which is the machinery a muscle")
    print("      cell contracts with; the ribosome and the elongation")
    print("      factor, which build it; GAPDH, which powers it; and")
    print("      both chains of ferritin, because these cells store a")
    print("      lot of iron.")
    print("      Notice ACTB, GAPDH and RPLP0 in that list. They are")
    print("      the three housekeeping genes of Exercise 18.3, and")
    print("      they are there for the reason they are used: high,")
    print("      and doing a job every cell needs doing.")


def exercise_18_5():
    table, samples = counts()
    library = table.sum()
    cpm = table.loc[CRISPLD2] / library * 1e6
    joined = samples.assign(cpm=cpm[samples["sample"]].to_numpy())
    print("\n18.5  CRISPLD2, GeneID 83716, in all sixteen samples:")
    print(f"       {'donor':10} {'treatment':24} {'CPM':>9}")
    for _, row in joined.sort_values(["cell_line",
                                      "treatment"]).iterrows():
        print(f"       {row['cell_line']:10} {row['treatment']:24} "
              f"{row['cpm']:>9.1f}")
    paired = joined[joined["treatment"].isin(DEX)]
    means = paired.groupby("treatment")["cpm"].mean()
    fold = means["Dexamethasone"] / means["Untreated"]
    print(f"\n      Untreated {means['Untreated']:.1f} CPM against "
          f"{means['Dexamethasone']:.1f} with the steroid,")
    print(f"      a {fold:.1f} fold increase, and it goes the same way")
    print("      in every donor.")
    print("      That is a smaller change than FKBP5's, and it is the")
    print("      gene the paper is named after. Size is not")
    print("      importance: CRISPLD2 was worth a paper because of")
    print("      what it does and because variants in it are")
    print("      associated with how patients respond to steroids,")
    print("      not because it moved the furthest.")


if __name__ == "__main__":
    exercise_18_1()
    exercise_18_2()
    exercise_18_3()
    exercise_18_4()
    exercise_18_5()
