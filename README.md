# Bioinformatics with Python for Beginners

The code for *Bioinformatics with Python for Beginners: A Hands-On Guide to
DNA and Protein Sequence Analysis, Genomics and Biological Data Science with
Biopython* by Josh Festus. Every program printed in the book is here, with a
Jupyter notebook for each chapter, a script that downloads every dataset, the
three projects, and worked solutions to all 134 exercises.

Every listing in the book was run when the book was typeset, and the output
printed under it is what that run produced. If you type a program in and get
something different, compare it with the file here.

## Getting started

Install uv and the command line tools as Chapter 2 and Appendix A show, then:

```
git clone https://github.com/joshfestus/bioinformatics-python
cd bioinformatics-python
uv sync
```

`uv sync` reads `pyproject.toml` and `uv.lock` and installs the exact versions
the book was built with. After that, run any program from its own folder:

```
cd ch01
uv run gene.py
```

## What is here

| Folder | What it holds |
|---|---|
| `ch01` to `ch22` | The programs for each chapter that has them, named as the book names them, with the chapter's notebook beside them |
| `data/fetch.py` | Downloads the larger datasets, chapter by chapter |
| `projects/virus` | Project One, Chapter 23: tracking SARS-CoV-2 through its variants |
| `projects/resistance` | Project Two, Chapter 24: resistance genes in a Klebsiella genome |
| `projects/airway` | Project Three, Chapter 25: an RNA-seq study from counts to biology |
| `solutions` | One program per chapter that prints each exercise's answer and the reasoning behind it |

Chapter 3, the molecular biology chapter, and Chapter 26 have no programs of
their own, so they have no folder.

## Getting the data

Small datasets whose terms allow it sit in the chapter folders. Everything
larger is downloaded from where it was published, never committed here:

```
uv run data/fetch.py --list
uv run data/fetch.py ch15
```

The first command lists every dataset and its source, and the second fetches
what one chapter needs. Each project downloads its own data with the
`download.py` in its folder. Appendix B of the book lists every dataset with
its source and its terms of use.

## Solutions

Try each exercise before you look. Then run the chapter's solution from the
`solutions` folder:

```
cd solutions
uv run ch05.py
```

Where an exercise asks you to work something out by hand, the program gives
the hand answer first and then checks it with code.

## Versions

The book's output was produced with Python 3.14.7, BLAST+ 2.16.0, MAFFT
7.525, FastQC 0.12.1, minimap2 2.27, samtools 1.22.1 and bcftools 1.22, with
the Python packages pinned in `pyproject.toml` and `uv.lock`. When a later
release changes an example, the code here is updated.

## License

The code is MIT licensed, as the book's copyright page says: use it, change it
and share it, even commercially, as long as the license goes with it. See
`LICENSE`. The datasets keep their own terms, listed in Appendix B.
