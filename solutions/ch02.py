"""Solutions to the exercises of Chapter 2.

Run from this folder: uv run ch02.py
These exercises are about the workstation rather than about biology, so
most of the answers are commands to type. This program prints them, and
checks the ones it can check on your own machine.
"""

import shutil
import subprocess
from pathlib import Path

REPO = Path("..").resolve()


def exercise_2_1():
    """Moving around with pwd, cd and ls."""
    print("2.1  The three commands, from your home folder:\n")
    print("       cd bioinformatics-python/ch10")
    print("       pwd")
    print("       ls")
    print("       head -n 3 globins.fasta\n")
    listing = sorted(p.name for p in (REPO / "ch10").iterdir()
                     if p.is_file())
    print(f"     ch10 holds {len(listing)} files: {', '.join(listing)}")
    first = (REPO / "ch10" / "globins.fasta").read_text().splitlines()[:3]
    for line in first:
        print(f"       {line[:62]}")
    print("     `pwd` prints where you are, which is the question to ask")
    print("     first whenever a program cannot find its data.")


def exercise_2_2():
    """Adding pysam to the setup check."""
    print("\n2.2  Add the name to the list and rerun:\n")
    print('       PACKAGES = ["Bio", "numpy", "pandas", "matplotlib",')
    print('                   "pysam"]\n')
    print("     At Chapter 2 it reports pysam as missing, and that is")
    print("     correct: the book installs each package in the chapter")
    print("     that first needs one, and pysam is not needed until")
    print("     Chapter 15. A missing package is only a problem when a")
    print("     program you are running asks for it.")


def exercise_2_3():
    """One cell, then four, run out of order."""
    print("\n2.3  Split into four cells, the order you run them in")
    print("     decides what works. The base counts cell fails with")
    print("     NameError if the reading cell has not been run, because")
    print("     `gene` does not exist yet. Running the sickle change")
    print("     twice is harmless here, since it rebuilds `sickle` from")
    print("     `gene` each time, but a cell that changed `gene` itself")
    print("     would give a different answer the second time.")
    print("     A notebook keeps its variables between cells, which is")
    print("     what makes it convenient and what makes a notebook run")
    print("     out of order hard to trust. Restart and run all before")
    print("     believing a notebook's output.")


def exercise_2_4():
    """A project folder of your own."""
    print("\n2.4  The folder and the one change the program needs:\n")
    print("       mkdir -p ~/globin/data ~/globin/scripts")
    print("       mkdir -p ~/globin/results")
    print("       touch ~/globin/notes.md")
    print("       cp hbb.txt ~/globin/data/")
    print("       cp gene.py ~/globin/scripts/\n")
    print("     Then in scripts/gene.py the first line becomes:\n")
    print('       gene = open("../data/hbb.txt").read().strip()\n')
    print("     The path is relative to where you run the program from,")
    print("     not to where the program is saved, which is the thing")
    print("     that catches everyone once.")


def exercise_2_5():
    """Reading the lock file."""
    print("\n2.5  uv.lock records every package and what it depends on.")
    lock = REPO / "uv.lock"
    if not lock.exists():
        print("     No uv.lock here yet; run `uv sync` first.")
        return
    text = lock.read_text()
    blocks = text.split("[[package]]")
    for block in blocks:
        # the first `name =` in a block is the package the block is about;
        # the others are the packages it depends on
        first = next((line for line in block.splitlines()
                      if line.startswith("name = ")), "")
        if first == 'name = "biopython"':
            needs = [line.split('"')[1]
                     for line in block.splitlines()
                     if line.strip().startswith("{ name =")]
            print(f"     Biopython depends on {len(needs)}: "
                  f"{', '.join(needs) if needs else 'nothing'}")
    depend_on_numpy = sum(1 for block in blocks
                          if '{ name = "numpy"' in block)
    print(f"     {depend_on_numpy} packages list numpy as a dependency.")
    print("     Biopython needs almost nothing, which is why it installs")
    print("     in seconds; numpy sits underneath most of the rest.")


def exercise_2_check():
    """What this machine actually has, for comparison."""
    print("\n     On this machine, for comparison:")
    for tool in ("python3", "git", "uv"):
        where = shutil.which(tool)
        version = "not found"
        if where:
            got = subprocess.run([where, "--version"],
                                 capture_output=True, text=True)
            version = got.stdout.strip() or got.stderr.strip()
        print(f"       {tool:8} {version}")


if __name__ == "__main__":
    exercise_2_1()
    exercise_2_2()
    exercise_2_3()
    exercise_2_4()
    exercise_2_5()
    exercise_2_check()
