"""Solutions to the exercises of Chapter 13.

Run from this folder: uv run ch13.py
The bootstrap in 13.3 resamples 100 times and takes a minute or two.
Exercise 13.5 needs MAFFT and the network.
"""

import random
import shutil
import subprocess
import tempfile
import urllib.error
import urllib.request
from io import StringIO
from pathlib import Path

from Bio import AlignIO, Phylo, SeqIO
from Bio.Phylo.Consensus import bootstrap_trees
from Bio.Phylo.TreeConstruction import (
    DistanceCalculator,
    DistanceTreeConstructor,
)

CH13 = Path("../ch13")
UNIPROT = "https://rest.uniprot.org/uniprotkb"
# the alpha and beta type globins of Chapter 11's BLAST hits
FAMILY = {"P68871": "HBB", "P02042": "HBD", "P69891": "HBG1",
          "P69892": "HBG2", "P02100": "HBE1", "P69905": "HBA1",
          "P02008": "HBZ", "P09105": "HBQ1", "P02144": "MYG"}


def alignment():
    found = AlignIO.read(CH13 / "globins_aligned.fasta", "fasta")
    for record in found:
        record.id = record.id.split("|")[2]
    return found


def build(method, measure="blosum62"):
    constructor = DistanceTreeConstructor(DistanceCalculator(measure))
    data = alignment()
    tree = getattr(constructor, method)(
        DistanceCalculator(measure).get_distance(data))
    tree.root_with_outgroup("MYG_HUMAN")
    tree.ladderize()
    return tree


def leaves_under(tree, names):
    """Whether these names sit together on a branch of their own."""
    try:
        clades = [tree.find_any(name) for name in names]
    except Exception:
        return False
    if any(clade is None for clade in clades):
        return False
    return tree.is_monophyletic(clades) is not False


def drawn(tree):
    out = StringIO()
    Phylo.draw_ascii(tree, file=out, column_width=62)
    return out.getvalue()


def exercise_13_1():
    neighbor = build("nj")
    upgma = build("upgma")
    print("13.1  Neighbor joining against UPGMA, same distances:\n")
    print(drawn(upgma))
    groups = {
        "human with chimpanzee": ["HBB_HUMAN", "HBB_PANTR"],
        "mouse with rat": ["HBB1_MOUSE", "HBB1_RAT"],
        "dog with the primates": ["HBB_CANLF", "HBB_HUMAN",
                                  "HBB_PANTR"],
        "the three alpha globins": ["HBA_HUMAN", "HBA_MOUSE",
                                    "HBA_CHICK"],
    }
    print(f"       {'group':26} {'neighbor joining':>17} {'UPGMA':>8}")
    for label, names in groups.items():
        print(f"       {label:26} "
              f"{str(leaves_under(neighbor, names)):>17} "
              f"{str(leaves_under(upgma, names)):>8}")
    print("      UPGMA assumes every lineage changes at the same rate,")
    print("      so it places a leaf by its average distance to")
    print("      everything else. Where a lineage has changed faster")
    print("      than its relatives, UPGMA reads that extra distance")
    print("      as extra time and pushes the leaf towards the root.")
    print("      Neighbor joining does not assume a clock, which is")
    print("      why it is the method the chapter uses and the one")
    print("      used in practice.")


def exercise_13_2():
    print("\n13.2  The identity distance instead of BLOSUM62:")
    for measure in ("blosum62", "identity"):
        tree = build("nj", measure)
        with_primates = leaves_under(
            tree, ["HBB_CANLF", "HBB_HUMAN", "HBB_PANTR"])
        with_laurasiatheria = leaves_under(
            tree, ["HBB_CANLF", "HBB_PIG", "HBB_BOVIN", "HBB_HORSE"])
        print(f"       {measure:9}: dog with the primates "
              f"{with_primates}, dog with")
        print(f"                  pig, cow and horse "
              f"{with_laurasiatheria}")
    print("      The identity distance counts differences and nothing")
    print("      else. BLOSUM62 weighs them, so a conservative swap")
    print("      costs less than a radical one. They give the same")
    print("      answer here, which is reassuring: a result that")
    print("      changed with the distance measure would be a result")
    print("      about the measure rather than about the globins.")
    print("      The dog sitting with the primates is wrong either")
    print("      way, and Chapter 13 says why: one gene of 147")
    print("      residues does not hold enough signal to resolve")
    print("      branches that are close together in time.")


def exercise_13_3():
    print("\n13.3  More groups, 100 bootstrap replicates each:")
    random.seed(1)
    data = alignment()
    constructor = DistanceTreeConstructor(
        DistanceCalculator("blosum62"), "nj")
    trees = list(bootstrap_trees(data, 100, constructor))
    for tree in trees:
        tree.root_with_outgroup("MYG_HUMAN")
    groups = {
        "human and chimpanzee": ["HBB_HUMAN", "HBB_PANTR"],
        "dog with the primates": ["HBB_CANLF", "HBB_HUMAN",
                                  "HBB_PANTR"],
        "placental mammals": ["HBB_HUMAN", "HBB_PANTR", "HBB1_MOUSE",
                              "HBB1_RAT", "HBB_BOVIN", "HBB_PIG",
                              "HBB_HORSE", "HBB_CANLF"],
        "all mammals, with platypus": ["HBB_HUMAN", "HBB_PANTR",
                                       "HBB1_MOUSE", "HBB1_RAT",
                                       "HBB_BOVIN", "HBB_PIG",
                                       "HBB_HORSE", "HBB_CANLF",
                                       "HBB_ORNAN"],
        "amniotes with the frog": ["HBB_HUMAN", "HBB_PANTR",
                                   "HBB1_MOUSE", "HBB1_RAT",
                                   "HBB_BOVIN", "HBB_PIG",
                                   "HBB_HORSE", "HBB_CANLF",
                                   "HBB_ORNAN", "HBB_CHICK",
                                   "HBB1_XENLA"],
    }
    for label, names in groups.items():
        count = sum(leaves_under(tree, names) for tree in trees)
        print(f"       {label:28} {count:3} of 100")
    print("      A group in nearly every replicate is one the data")
    print("      supports. The mammals, with and without the platypus,")
    print("      come out above ninety, and they are true groupings.")
    print("      Dog with the primates comes out at 82 and is wrong,")
    print("      which is the warning in this exercise: bootstrap")
    print("      support measures how consistently the data says a")
    print("      thing, not whether the thing is true. A systematic")
    print("      error repeats in every replicate.")
    print("      The amniotes with the frog almost never appear, and")
    print("      that is the same mistake seen from the other side.")
    print("      In these trees the frog branches outside the")
    print("      zebrafish, so a group holding the frog and the land")
    print("      vertebrates but not the fish cannot form. One gene of")
    print("      147 residues is not enough to place a branch that old.")


def exercise_13_4():
    tree = Phylo.read(CH13 / "globins.nwk", "newick")
    before = [tip.name for tip in tree.get_terminals()]
    # pick a node with two branches that actually differ, not the
    # outgroup leaf, or the reversal changes nothing to look at
    node = next(clade for clade in tree.get_nonterminals()
                if len(clade.clades) == 2
                and all(len(child.get_terminals()) > 1
                        for child in clade.clades))
    node.clades.reverse()
    after = [tip.name for tip in tree.get_terminals()]
    moved = sum(a != b for a, b in zip(before, after))
    print(f"\n13.4  Reversing one node moves {moved} of {len(before)} "
          "leaves in the drawing.")
    print(f"       before: {', '.join(before[:4])} ...")
    print(f"       after:  {', '.join(after[:4])} ...")
    print("      Nothing about the history changed. A tree records")
    print("      which leaves share an ancestor and how far apart they")
    print("      are, and swapping the two children of a node changes")
    print("      neither. The drawing looks different because a")
    print("      drawing has to put one branch above the other, and")
    print("      that choice carries no information.")
    print("      This matters when reading a published tree: two trees")
    print("      that look unalike may be the same tree, and the test")
    print("      is the groups, not the order down the page.")


def exercise_13_5():
    print("\n13.5  A tree of the human globin family:")
    if not shutil.which("mafft"):
        print("      MAFFT is not on the PATH; see Appendix A.")
        return
    records = []
    for accession, name in FAMILY.items():
        try:
            request = urllib.request.Request(
                f"{UNIPROT}/{accession}.fasta",
                headers={"User-Agent": "python"})
            with urllib.request.urlopen(request, timeout=120) as reply:
                record = SeqIO.read(StringIO(reply.read().decode()),
                                    "fasta")
        except (urllib.error.URLError, TimeoutError) as problem:
            print(f"      Could not download {accession} ({problem}).")
            return
        record.id = name
        record.description = ""
        records.append(record)
    with tempfile.TemporaryDirectory() as work:
        source = Path(work) / "family.fasta"
        SeqIO.write(records, source, "fasta")
        done = subprocess.run(["mafft", "--quiet", "--auto", str(source)],
                              capture_output=True, text=True, check=True)
        aligned = AlignIO.read(StringIO(done.stdout), "fasta")
    calculator = DistanceCalculator("blosum62")
    tree = DistanceTreeConstructor(calculator, "nj").build_tree(aligned)
    tree.root_with_outgroup("MYG")
    tree.ladderize()
    print(drawn(tree))
    print("      The two clusters are the two gene families. HBA1,")
    print("      HBZ and HBQ1 are the alpha type globins on")
    print("      chromosome 16; HBB, HBD, HBG1, HBG2 and HBE1 are the")
    print("      beta type globins on chromosome 11, in a cluster")
    print("      Chapter 5 found spread over 44 kb.")
    print("      Inside the beta family the fetal pair HBG1 and HBG2")
    print("      sit together, as two copies of the same gene should,")
    print("      and the embryonic HBE1 and HBZ branch away from their")
    print("      adult relatives. The tree of the gene family tracks")
    print("      the order the duplications happened in, which is a")
    print("      different history from the tree of species in")
    print("      Exercise 13.3 even though both are built the same way.")


if __name__ == "__main__":
    exercise_13_1()
    exercise_13_2()
    exercise_13_3()
    exercise_13_4()
    exercise_13_5()
