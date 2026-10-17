"""Solutions to the exercises of Chapter 11.

Run from this folder: uv run ch11.py
Needs BLAST+ on the PATH and the database Chapter 11 built in ../ch11.
Exercise 11.3 downloads one protein from UniProt.
"""

import shutil
import subprocess
import urllib.error
import urllib.request
from pathlib import Path

from Bio import SeqIO
from Bio.Align import substitution_matrices

CH04, CH11 = Path("../ch04"), Path("../ch11")
BLOSUM62 = substitution_matrices.load("BLOSUM62")
UNIPROT = "https://rest.uniprot.org/uniprotkb"
FORMAT = "6 qseqid sseqid pident length evalue bitscore"


def have(tool):
    return shutil.which(tool) is not None


def blast(program, query, database, extra=()):
    """Run a BLAST program and return its tabular lines."""
    command = [program, "-query", str(query), "-db", str(database),
               "-outfmt", FORMAT, *extra]
    done = subprocess.run(command, capture_output=True, text=True,
                          cwd=CH11)
    if done.returncode != 0:
        raise RuntimeError(done.stderr.strip()[:200])
    return [line.split("\t") for line in done.stdout.splitlines()]


def seeds(word, threshold):
    """Chapter 11's seed finder, with the word size as an argument."""
    beta = str(SeqIO.read(CH11 / "hbb.fasta", "fasta").seq)
    alpha = str(SeqIO.read(CH11 / "hba_human.fasta", "fasta").seq)
    found = []
    for i in range(len(beta) - word + 1):
        for j in range(len(alpha) - word + 1):
            score = sum(BLOSUM62[x][y]
                        for x, y in zip(beta[i:i + word],
                                        alpha[j:j + word]))
            if score >= threshold:
                found.append((i, j, score))
    return found


def exercise_11_1():
    print("11.1  Seeds at two word sizes:")
    for word, threshold in ((3, 11), (2, 8)):
        found = seeds(word, threshold)
        near = [s for s in found if abs(s[0] - s[1]) <= 10]
        share = len(near) / len(found) if found else 0
        print(f"       word {word}, threshold {threshold}: "
              f"{len(found):5} seeds, {len(near):4} near the diagonal "
              f"({share:.0%})")
    print("      A shorter word finds far more seeds and a smaller")
    print("      share of them mean anything, because two residues")
    print("      match by chance far more often than three do. Every")
    print("      seed costs an extension attempt, so BLAST would spend")
    print("      most of its time chasing noise. Three is the length")
    print("      where a protein word is rare enough to be worth")
    print("      following and short enough to be found at all.")


def exercise_11_2():
    print("\n11.2  tblastn: a protein query against translated DNA.")
    if not have("tblastn") or not have("makeblastdb"):
        print("      BLAST+ is not on the PATH; see Appendix A.")
        return
    subprocess.run(["makeblastdb", "-in", str(CH04 / "globin_cds.fasta"),
                    "-dbtype", "nucl", "-out", "globin_cds"],
                   capture_output=True, text=True, cwd=CH11, check=True)
    command = ["tblastn", "-query", "hbb.fasta", "-db", "globin_cds",
               "-outfmt", "6 sseqid pident length evalue sframe"]
    done = subprocess.run(command, capture_output=True, text=True,
                          cwd=CH11, check=True)
    print(f"       {'subject':28} {'ident':>6} {'len':>4} "
          f"{'E value':>10} {'frame':>6}")
    for line in done.stdout.splitlines():
        subject, pident, length, evalue, frame = line.split("\t")
        print(f"       {subject[:28]:28} {pident:>6} {length:>4} "
              f"{evalue:>10} {frame:>6}")
    for name in CH11.glob("globin_cds.n*"):
        name.unlink()
    print("      Every hit is in frame 1, because these records are")
    print("      coding sequences that begin at their first base. A")
    print("      search against a genome would find hits in all six")
    print("      frames, and the frame is how you know which strand")
    print("      and which reading the protein came from.")


def exercise_11_3():
    print("\n11.3  Human beta galactosidase, P16278, against E. coli:")
    if not have("blastp") or not have("makeblastdb"):
        print("      BLAST+ is not on the PATH; see Appendix A.")
        return
    query = CH11 / "P16278.fasta"
    if not query.exists():
        try:
            request = urllib.request.Request(
                f"{UNIPROT}/P16278.fasta",
                headers={"User-Agent": "python"})
            with urllib.request.urlopen(request, timeout=120) as reply:
                query.write_text(reply.read().decode())
        except (urllib.error.URLError, TimeoutError) as problem:
            print(f"      Could not download it ({problem}).")
            return
    subprocess.run(["makeblastdb", "-in", "../ch07/lac_proteins.fasta",
                    "-dbtype", "prot", "-out", "lac"],
                   capture_output=True, text=True, cwd=CH11, check=True)
    hits = blast("blastp", "P16278.fasta", "lac", ["-evalue", "10"])
    names = {"NP_414876.1": "lacA, the acetyltransferase",
             "NP_414877.1": "lacY, the membrane permease",
             "NP_414878.1": "lacZ, beta-galactosidase",
             "NP_414879.3": "lacI, the repressor"}
    if not hits:
        print("       Nothing at all below an E value of 10.")
    for hit in hits:
        print(f"       {names.get(hit[1], hit[1]):28} identity "
              f"{float(hit[2]):5.1f}%, E {float(hit[4]):8.2g}")
    lacz = [h for h in hits if h[1] == "NP_414878.1"]
    best_lacz = min((float(h[4]) for h in lacz), default=None)
    print("      The enzyme we were looking for is lacZ, and its best")
    if best_lacz is None:
        print("      hit does not appear at all below an E value of 10.")
    else:
        print(f"      hit here is E = {best_lacz:.2g}, which is no")
        print("      better than chance would give.")
    print("      Human GLB1 and E. coli LacZ cut the same bond on the")
    print("      same sugar, and BLAST cannot see any relationship")
    print("      between them: they belong to different families of")
    print("      glycoside hydrolase that arrived at the job")
    print("      separately. Doing the same chemistry is not evidence")
    print("      of common ancestry.")
    print("      The one hit below 1 is to lacY, a membrane transport")
    print("      protein, which shares nothing with a soluble enzyme")
    print("      except stretches of hydrophobic residues. That is")
    print("      what a chance hit looks like, and the reason to read")
    print("      what a subject actually is before believing a hit.")
    query.unlink()
    for name in CH11.glob("lac.p*"):
        name.unlink()


def exercise_11_4():
    print("\n11.4  The same search at four thresholds:")
    if not have("blastp"):
        print("      BLAST+ is not on the PATH; see Appendix A.")
        return
    for threshold in ("0.001", "0.1", "10", "100"):
        hits = blast("blastp", "hbb.fasta", "human",
                     ["-evalue", threshold])
        names = [hit[1] for hit in hits]
        has_ngb = any("NGB" in name or "NEUG" in name for name in names)
        print(f"       E <= {threshold:>5}: {len(hits):3} hits"
              f"{'  including neuroglobin' if has_ngb else ''}")
    print("      Neuroglobin is the twilight zone case of Chapter 11:")
    print("      23 percent identical to beta globin, which is real")
    print("      homology that the statistics can barely distinguish")
    print("      from chance. The chapter found it only with a")
    print("      threshold of 100, because the first pass scores")
    print("      without the composition correction that later rescues")
    print("      it. A threshold is not a measure of truth; it decides")
    print("      how much noise you will read through.")


def exercise_11_5():
    print("\n11.5  Myoglobin against the human proteins:")
    if not have("blastp"):
        print("      BLAST+ is not on the PATH; see Appendix A.")
        return
    hits = blast("blastp", "myg.fasta", "human", ["-evalue", "10"])
    for hit in hits[:12]:
        print(f"       {hit[1][:34]:34} identity {hit[2]:>6}%, "
              f"E {hit[4]:>9}")
    names = " ".join(hit[1] for hit in hits)
    wide = blast("blastp", "myg.fasta", "human", ["-evalue", "1000"])
    anywhere = any("HBB" in hit[1] for hit in wide)
    present = "is" if "HBB" in names else "is not"
    print(f"      {len(hits)} hits below an E value of 10. Beta globin")
    print(f"      {present} among them,")
    ever = "it appears" if anywhere else "it does not appear"
    print(f"      and {ever} even at an E value of 1000.")
    print("      Look at what it did find: cytoglobin, and the alpha,")
    print("      zeta, theta, mu, epsilon and both gamma chains. The")
    print("      gamma and epsilon chains are beta family globins,")
    print("      more than seventy percent identical to beta globin")
    print("      itself, so myoglobin finds beta globin's close")
    print("      relatives and never beta globin.")
    print("      That is what a detection floor looks like from")
    print("      underneath. These proteins are all about a quarter")
    print("      identical to myoglobin, which is where the statistics")
    print("      stop being able to separate homology from chance, and")
    print("      on which side of the line a given pair lands is close")
    print("      to arbitrary. The structures in Chapter 21 settle what")
    print("      the sequences cannot.")


if __name__ == "__main__":
    exercise_11_1()
    exercise_11_2()
    exercise_11_3()
    exercise_11_4()
    exercise_11_5()
