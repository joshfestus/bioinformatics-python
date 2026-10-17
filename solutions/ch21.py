"""Solutions to the exercises of Chapter 21.

Run from this folder: uv run ch21.py
Reads the structures in ../data. Exercise 21.5 downloads one AlphaFold
model and needs the network.
"""

import urllib.error
import urllib.request
from pathlib import Path

import matplotlib
import numpy as np
from Bio.PDB import MMCIFParser, PDBParser
from Bio.PDB.DSSP import residue_max_acc
from Bio.PDB.SASA import ShrakeRupley

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

DATA = Path("../data")
LARGEST = residue_max_acc["Sander"]
ALPHAFOLD = "https://alphafold.ebi.ac.uk/files"
CHARGED = {"ARG", "LYS", "ASP", "GLU", "HIS"}


def model_of(code):
    path = DATA / f"{code}.cif"
    return MMCIFParser(QUIET=True).get_structure(code, path)[0]


def exercise_21_1():
    model = model_of("2HHB")
    print("21.1  The iron's two histidines in every chain of 2HHB:")
    print(f"       {'chain':7} {'kind':7} {'proximal':>12} "
          f"{'distal':>10}")
    for chain_id, proximal, distal, kind in (("A", 87, 58, "alpha"),
                                             ("C", 87, 58, "alpha"),
                                             ("B", 92, 63, "beta"),
                                             ("D", 92, 63, "beta")):
        chain = model[chain_id]
        heme = next(r for r in chain if r.get_resname() == "HEM")
        iron = heme["FE"]
        near = iron - chain[proximal]["NE2"]
        far = iron - chain[distal]["NE2"]
        print(f"       {chain_id:7} {kind:7} {near:>9.2f} A "
              f"{far:>8.2f} A")
    print("      The alpha chains hold their iron the same way the")
    print("      beta chains do: the proximal histidine bonded at")
    print("      about two angstroms and the distal one four or more")
    print("      away, leaving the gap where oxygen binds.")
    print("      The residue numbers differ, 87 and 58 against 92 and")
    print("      63, because the two chains are different lengths.")
    print("      The arrangement does not. Chapter 13's tree said")
    print("      these two genes are cousins from one ancestor, and")
    print("      this is that ancestry seen in three dimensions: the")
    print("      sequences drifted, the geometry around the iron did")
    print("      not, because a globin that lost it would not work.")


def exercise_21_2():
    model = model_of("2HHB")
    chain = model["B"]
    residues = [r for r in chain if r.id[0] == " " and "CA" in r]
    coordinates = np.array([r["CA"].coord for r in residues])
    numbers = [r.id[1] for r in residues]
    close = []
    for index in range(len(residues) - 4):
        here = coordinates[index]
        gaps = [np.linalg.norm(here - coordinates[index + k])
                for k in (3, 4)]
        close.append(max(gaps) < 6.5)  # a helix keeps i to i+4 near
    runs, start = [], None
    for index, helical in enumerate(close):
        if helical and start is None:
            start = index
        elif not helical and start is not None:
            if index - start >= 4:
                runs.append((numbers[start], numbers[index + 3]))
            start = None
    print("\n21.2  Helices found from the contact map alone, by")
    print("      looking for residues close to the one four along:")
    for first, last in runs:
        print(f"       {first:>4} to {last:>4}  ({last - first + 1} "
              "residues)")
    print(f"      {len(runs)} stretches, against the eight helices")
    print("      UniProt lists for P68871 and names A to H. The")
    print("      agreement is close but not exact, and both kinds of")
    print("      error are worth seeing: a short helix can fall below")
    print("      the rule and vanish, and a kinked one can come back")
    print("      as two. The boundaries drift by a residue or two as")
    print("      well, because a helix frays at its ends.")
    print("      The stretch around 60 is helix E, which carries the")
    print("      distal histidine of Exercise 21.1.")
    print("      Nothing here knows what a helix is. The pattern is")
    print("      simply that a residue sits near the one four along,")
    print("      which is what winding into a spiral does, and the")
    print("      contact map makes it visible as a thick band beside")
    print("      the diagonal.")


def exercise_21_3():
    model = model_of("2HHB")
    for chain in list(model):
        if chain.id not in "ABCD":
            model.detach_child(chain.id)
    ShrakeRupley().compute(model, level="R")
    exposure = []
    for residue in model["B"]:
        if residue.id[0] != " ":
            continue
        name = residue.get_resname()
        if name not in LARGEST:
            continue
        share = residue.sasa / LARGEST[name]
        exposure.append((share, residue.id[1], name))
    exposure.sort(reverse=True)
    top = exposure[:20]
    charged = sum(1 for _, _, name in top if name in CHARGED)
    print(f"\n21.3  Chain B has {len(exposure)} residues with an area.")
    print(f"       most exposed: {top[0][2]}{top[0][1]} at "
          f"{top[0][0]:.0%}")
    print(f"       least exposed: {exposure[-1][2]}{exposure[-1][1]} at "
          f"{exposure[-1][0]:.0%}")
    print(f"       median exposure: "
          f"{np.median([s for s, _, _ in exposure]):.0%}")
    print(f"      Of the twenty most exposed, {charged} are charged:")
    print("       " + ", ".join(f"{name}{number}"
                                for _, number, name in top
                                if name in CHARGED))
    print("      A protein that dissolves in water keeps its charges")
    print("      outside, where they can meet the water, and its oily")
    print("      residues inside. That is most of what folding is,")
    print("      and it is why Chapter 9's hydropathy plot could find")
    print("      membrane helices: a stretch with nowhere to hide its")
    print("      oily residues is a stretch that lives in a membrane.")
    print("      It is also why the sickle mutation matters. Chapter")
    print("      21 found valine 6 is 59 percent exposed, an oily")
    print("      residue sitting on a wet surface, which is exactly")
    print("      the kind of mistake this rule says should not")
    print("      happen, and it sticks to the first oily pocket it")
    print("      meets.")


def exercise_21_4():
    print("\n21.4  What the other three valines touch.")
    print("      2HBS holds two whole hemoglobin molecules, chains A")
    print("      to D and E to H. The four sickle valines are the")
    print("      residue 6 of the four beta chains, B, D, F and H.")
    model = model_of("2HBS")
    chains = [c.id for c in model]
    print(f"       chains present: {', '.join(chains)}")
    valines = [(c.id, c[6].get_resname()) for c in model
               if 6 in [r.id[1] for r in c if r.id[0] == " "]]
    print("       residue 6 in each chain: "
          + ", ".join(f"{c}:{n}" for c, n in valines))
    atoms = {c.id: [a for r in c if r.id[0] == " " for a in r]
             for c in model}
    print("\n       what each valine 6 reaches, within 4.5 angstroms:")
    for chain_id, name in valines:
        if name != "VAL":
            continue
        mine = [a for a in model[chain_id][6]]
        partners = set()
        for other, theirs in atoms.items():
            if other == chain_id:
                continue
            for atom in mine:
                for neighbor in theirs:
                    if (atom - neighbor) <= 4.5:
                        partners.add(other)
                        break
        found = ", ".join(sorted(partners)) if partners else "nothing"
        print(f"         chain {chain_id}: {found}")
    print("      Only one of the four, chain H, touches another")
    print("      chain inside this file, and that is the contact")
    print("      Chapter 21 measured against chain B. The other")
    print("      three reach nothing here at all.")
    print("      They are not free in the fiber. Their partners are")
    print("      molecules in the neighboring unit cells of the")
    print("      crystal, which the file records only as the")
    print("      symmetry operations needed to place them.")
    print("      Two ways to see it. Ask the PDB for the biological")
    print("      assembly, which applies those operations and hands")
    print("      back the molecules already placed. Or apply them")
    print("      yourself from the cell parameters in the file, which")
    print("      Biopython does not do for you and which is the")
    print("      reason the first route exists.")
    print("      The crystallographic lesson is worth keeping: a PDB")
    print("      file is one asymmetric unit, not one molecule and")
    print("      not the biological thing. Chapter 21 checked that")
    print("      2HBS held two tetramers before trusting any contact")
    print("      it found, and that check is the habit.")


def exercise_21_5():
    print("\n21.5  pLDDT along a predicted structure.")
    accession = "Q8N4S0"  # a human protein of unknown function
    path = Path(f"AF-{accession}-F1.pdb")
    if not path.exists():
        try:
            # the version suffix moves; v6 is current, and the
            # service answers 403 without a User-Agent header
            url = f"{ALPHAFOLD}/AF-{accession}-F1-model_v6.pdb"
            request = urllib.request.Request(
                url, headers={"User-Agent": "python"})
            with urllib.request.urlopen(request, timeout=120) as reply:
                path.write_bytes(reply.read())
        except (urllib.error.URLError, TimeoutError) as problem:
            print(f"      Could not download {accession} ({problem}).")
            return
    known = PDBParser(QUIET=True).get_structure(
        "hbb", DATA / "AF-P68871-F1.pdb")[0]["A"]
    unknown = PDBParser(QUIET=True).get_structure(
        accession, path)[0]["A"]
    figure, axes = plt.subplots(figsize=(6, 3.2))
    for label, chain, shade in (("beta globin", known, "#222222"),
                                (accession, unknown, "#888888")):
        scores = [r["CA"].bfactor for r in chain if "CA" in r]
        axes.plot(range(1, len(scores) + 1), scores, color=shade,
                  lw=1, label=label)
        low = sum(1 for s in scores if s < 70)
        print(f"       {label:12} {len(scores):4} residues, mean pLDDT "
              f"{np.mean(scores):5.1f}, {low:4} below 70")
    axes.axhline(70, color="#bbbbbb", lw=0.8, ls=(0, (4, 3)))
    axes.set_xlabel("Residue")
    axes.set_ylabel("pLDDT")
    axes.legend(frameon=False, fontsize=8)
    figure.savefig("ch21-plddt.png", dpi=110, bbox_inches="tight")
    plt.close(figure)
    path.unlink()
    print("      Beta globin is a protein with thousands of relatives")
    print("      and dozens of solved structures, and AlphaFold is")
    print("      confident about nearly every residue of it.")
    print("      The comparison protein is not, and the plot shows")
    print("      where a model gets unsure: the ends, and any stretch")
    print("      with no fixed shape to predict. A low pLDDT is not")
    print("      always a failure. Some proteins genuinely have")
    print("      disordered regions, and the model saying so is the")
    print("      model being right.")
    print("      Saved as ch21-plddt.png.")


if __name__ == "__main__":
    exercise_21_1()
    exercise_21_2()
    exercise_21_3()
    exercise_21_4()
    exercise_21_5()
