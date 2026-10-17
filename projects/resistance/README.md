# Which antibiotics still work, and what can spread

A resistance profile for a hospital isolate of *Klebsiella pneumoniae*,
built from its public genome and NCBI's Reference Gene Catalog.

This is Project Two of *Bioinformatics with Python for Beginners*
(Chapter 24). A full run takes about five and a half minutes on two cores.

## What it finds

The isolate is HS11286: one chromosome of 5,333,942 bases and six
plasmids, 5,779 proteins in all.

- 19 resistance genes in 8 classes, including blaKPC-2, the carbapenemase
  that defeats the drugs held in reserve for this species.
- 2 resistance mutations, gyrA S83I and parC S80I, found by reading the
  residue at the position the catalog names rather than by trusting a
  BLAST hit. Of 41 hits to point-mutation references, only these 2 carry
  the mutation.
- 15 of the 19 genes sit on plasmids and can move to another bacterium.
  The 4 on the chromosome are ones this species carries anyway.
- 11 genes sit together in a resistance island on plasmid pKPHS3, nearly
  all beside a transposase or integrase, built around a class 1 integron
  whose cassettes are dfrA12 and aadA2.

## Running it

```
uv run download.py              # genome and catalog, about 19 MB
snakemake --cores 2 all         # every step, to report.md
```

## The steps

| File | What it does | Writes |
|---|---|---|
| `download.py` | Fetches the 7 replicons and the catalog | `genome.gb` |
| `describe.py` | Sizes, GC content and protein counts | `replicons.tsv` |
| `proteins.py` | Every annotated protein as FASTA | `proteins.fasta` |
| `search.sh` | blastp against the catalog | `hits.tsv` |
| `filter.py` | Applies the thresholds, sorts by kind | `resistance.tsv` |
| `point.py` | Checks the point mutations by alignment | `mutations.tsv` |
| `place.py` | Which replicon, and what lies nearby | `profile.tsv` |
| `report.py` | The profile as Markdown | `report.md` |

## Thresholds

A hit counts when it is at least 90 percent identical to the reference
protein over at least 90 percent of that protein's length. Both values
live in `config.yaml` and are printed in the report with the catalog
version, because they decide the answer. Point mutation rules are applied
only for *Klebsiella pneumoniae*: a rule written for *Escherichia* is not
evidence about this organism.

## Data

The genome is NCBI RefSeq accessions NC_016845.1 (chromosome) and
NC_016838.1, NC_016846.1, NC_016839.1, NC_016840.1, NC_016847.1 and
NC_016841.1 (plasmids pKPHS1 to pKPHS6), public domain. The catalog is
NCBI's Reference Gene Catalog, downloaded by `download.py`, which also
records its version. Nothing is redistributed here.

## What it cannot say

This project found genes and mutations in DNA. It did not measure
resistance. A gene can be present and silent, or broken in a way the
search did not notice. The laboratory's susceptibility test is the
measurement; a genome explains and predicts it, and it is normal for the
two to disagree.
