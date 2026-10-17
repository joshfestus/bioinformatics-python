"""Solutions to the exercises of Chapter 9.

Run from this folder: uv run ch09.py
Exercise 9.2 downloads one protein from UniProt; the rest are local.
"""

import urllib.error
import urllib.request
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.SeqUtils.ProtParamData import kd

CH08, CH09 = Path("../ch08"), Path("../ch09")
WINDOW, THRESHOLD, MAX_GAP = 19, 1.6, 10
UNIPROT = "https://rest.uniprot.org/uniprotkb"


def segments_above(sequence, threshold=THRESHOLD, window=WINDOW):
    """The stretches whose hydropathy stays above a threshold."""
    scores = ProteinAnalysis(sequence).protein_scale(kd, window)
    centers = range(window // 2 + 1, len(sequence) - window // 2 + 1)
    above = [p for p, score in zip(centers, scores) if score > threshold]
    if not above:
        return []
    found = [[above[0], above[0]]]
    for position in above[1:]:
        if position - found[-1][1] <= MAX_GAP:
            found[-1][1] = position
        else:
            found.append([position, position])
    return found


def exercise_9_1():
    rows = []
    for record in SeqIO.parse(CH08 / "globins.fasta", "fasta"):
        protein = ProteinAnalysis(str(record.seq))
        rows.append({
            "name": record.id.split("|")[2],
            "length": len(record.seq),
            "weight": round(protein.molecular_weight() / 1000, 1),
            "pI": round(protein.isoelectric_point(), 2),
            "gravy": round(protein.gravy(), 3),
        })
    table = pd.DataFrame(rows).sort_values("pI").reset_index(drop=True)
    print("9.1  The sixteen globins, by isoelectric point:\n")
    print(table.to_string(index=False))
    betas = table[table["name"].str.startswith("HBB")]
    print(f"\n     The beta globins run pI {betas['pI'].min()} to "
          f"{betas['pI'].max()} and {betas['weight'].min()} to "
          f"{betas['weight'].max()} kDa.")
    print("     They are alike in size and weight, as proteins doing")
    print("     the same job in different species should be, and their")
    print("     isoelectric points differ by more than you might")
    print("     expect: a few charged residues swapped between species")
    print("     move the pI without changing what the protein does.")


def exercise_9_2():
    print("\n9.2  The beta-2 adrenergic receptor, P07550:")
    path = Path("P07550.fasta")
    if not path.exists():
        try:
            request = urllib.request.Request(
                f"{UNIPROT}/P07550.fasta",
                headers={"User-Agent": "python"})
            with urllib.request.urlopen(request, timeout=120) as reply:
                path.write_text(reply.read().decode())
        except (urllib.error.URLError, TimeoutError) as problem:
            print(f"     Could not download it ({problem}).")
            return
    sequence = str(SeqIO.read(path, "fasta").seq)
    found = segments_above(sequence)
    print(f"     {len(sequence)} residues, {len(found)} segments above "
          f"{THRESHOLD}:")
    for first, last in found:
        print(f"       {first:4} to {last:4}, {last - first + 1:3} "
              "residues")
    print("     This receptor crosses the membrane seven times, like")
    print("     rhodopsin, and the same method finds a similar number")
    print("     of them. It is the receptor that the albuterol of")
    print("     Chapters 18 to 20 acts on, so the same molecule turns")
    print("     up in the structural half of the book and the")
    print("     expression half.")
    path.unlink()


def exercise_9_3():
    globin = str(SeqIO.read(CH09 / "hbb_human.fasta", "fasta").seq)
    found = segments_above(globin)
    print(f"\n9.3  Human beta globin, {len(globin)} residues: "
          f"{len(found)} segments above {THRESHOLD}.")
    scores = ProteinAnalysis(globin).protein_scale(kd, WINDOW)
    print(f"     Its highest hydropathy anywhere is {max(scores):.2f}, "
          f"against a threshold of {THRESHOLD}.")
    print("     It should find none, and it finds none. Beta globin is")
    print("     a soluble protein that floats inside a red blood cell,")
    print("     so it has no membrane crossing helix to find. Running")
    print("     a method on a case where you know the answer is a")
    print("     negative control, and a method that found helices here")
    print("     would be one you could not trust on rhodopsin.")


def exercise_9_4():
    known = pd.read_csv(CH09 / "rhodopsin_helices.tsv", sep="\t",
                        comment="#")
    rhodopsin = str(SeqIO.read(CH09 / "rhodopsin.fasta", "fasta").seq)
    globin = str(SeqIO.read(CH09 / "hbb_human.fasta", "fasta").seq)
    print("\n9.4  Lowering the threshold, for rhodopsin's seven helices:")
    print(f"       {'threshold':>9} {'segments':>9} {'helices found':>14}"
          f" {'in globin':>10}")
    for threshold in (1.6, 1.4, 1.2, 1.0, 0.8):
        found = segments_above(rhodopsin, threshold)
        hit = sum(any(first <= row["end"] and last >= row["start"]
                      for first, last in found)
                  for _, row in known.iterrows())
        false_positives = len(segments_above(globin, threshold))
        print(f"       {threshold:>9} {len(found):>9} {hit:>14} "
              f"{false_positives:>10}")
    print("     The seventh helix needs the threshold down near 1.2,")
    print("     because it holds the lysine that binds retinal and a")
    print("     charged residue drags the average down. At 1.2 beta")
    print("     globin still gives nothing, so there is a window where")
    print("     all seven are found and the soluble protein stays")
    print("     quiet. Go on down to 0.8 and the method starts calling")
    print("     segments in a protein that has none. The margin is")
    print("     thin, and that is why real predictors use a model")
    print("     trained on known structures rather than one number.")


def exercise_9_5():
    beta = str(SeqIO.read(CH09 / "hbb_human.fasta", "fasta").seq)
    # the chapter's sickle_charge.py works on the mature protein, with
    # the initial methionine removed, so position 6 is position 6
    mature = beta[1:]
    sickle = mature[:5] + "V" + mature[6:]
    hemoglobin_c = mature[:5] + "K" + mature[6:]
    print("\n9.5  On the mature protein, as sickle_charge.py does, "
          "position 6 is")
    print(f"     {mature[5]} normally, {sickle[5]} in sickle "
          f"hemoglobin and {hemoglobin_c[5]} in hemoglobin C.\n")
    print(f"       {'variant':14} {'charge at pH 7.4':>17} {'pI':>7} "
          f"{'weight':>9}")
    for name, sequence in (("normal", mature), ("sickle", sickle),
                           ("hemoglobin C", hemoglobin_c)):
        protein = ProteinAnalysis(sequence)
        print(f"       {name:14} {protein.charge_at_pH(7.4):>17.2f} "
              f"{protein.isoelectric_point():>7.2f} "
              f"{protein.molecular_weight() / 1000:>8.2f} kDa")
    print("     Glutamic acid is negative, valine is neutral and lysine")
    print("     is positive, so each swap moves the charge the same way")
    print("     and hemoglobin C moves twice as far as sickle does.")
    print("     In Pauling's experiment the molecules are pulled")
    print("     through a field and separate by charge, so the three")
    print("     run at three speeds, with hemoglobin C the slowest")
    print("     towards the positive electrode. That separation is how")
    print("     these variants were found, and it is what made sickle")
    print("     cell disease the first molecular disease.")


if __name__ == "__main__":
    exercise_9_1()
    exercise_9_2()
    exercise_9_3()
    exercise_9_4()
    exercise_9_5()
