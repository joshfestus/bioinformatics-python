# Tracking a virus as it evolves

How SARS-CoV-2 changed between the first genome sequenced in Wuhan and the
variants that followed, and where in its genome the changes gathered.

This is Project One of *Bioinformatics with Python for Beginners*
(Chapter 23). It runs on a laptop in about a minute once the genomes are
downloaded.

## What it finds

114 public genomes from 13 Pango lineages, compared with the reference
genome NC_045512.2:

- Substitutions from the reference rise from a median of 5 in lineage A to
  122 in KP.3, collected in late 2024, along a molecular clock of about 28
  substitutions a year.
- The spike gene and ORF8 change most per thousand bases, at about 56 and
  60 changed positions; the replication machinery in ORF1ab changes at 19,
  and the envelope protein E least of all at 18.
- Each named variant's defining spike mutations come out of the data
  without being told: Alpha's A570D, D614G, P681H, T716I, S982A, D1118H and
  N501Y, Gamma's twelve, Omicron's thirty and the rest.
- The neighbor joining tree groups 10 of the 13 lineages into exactly one
  branch each. The three it does not are right not to: B and B.1 are the
  ancestral stock the others came from, and KP.3 sits inside JN.1 because
  it descends from it.
- Two genomes in the sample are identical to the reference at all 29,903
  positions although they were collected 14 months later, which is what a
  consensus pipeline that fills uncovered positions with the reference
  base produces.

## Running it

```
uv run download.py              # the reference and 130 genomes, about 4 MB
snakemake --cores 2 all         # every step, to report.md
```

`choose.py` is the step that picks the sample from NCBI Virus. It is not
part of the pipeline, because NCBI adds genomes every week and the sample
would change under you; `samples.tsv` holds the accessions this project
used. Run it to choose your own, into `my_samples.tsv`.

## The steps

| File | What it does | Writes |
|---|---|---|
| `choose.py` | Picks genomes per lineage from NCBI Virus | `my_samples.tsv` |
| `download.py` | Fetches the reference and every chosen genome | `genomes.fasta` |
| `clean.py` | Drops genomes too short or too uncertain | `clean.fasta` |
| `align.sh` | MAFFT, keeping the reference's numbering | `aligned.fasta` |
| `mutations.py` | Every difference, with its gene and amino acid change | `mutations.tsv` |
| `locate.py` | Changes per gene and per thousand bases | `per_gene.tsv` |
| `timeline.py` | The molecular clock and first appearances | `clock.tsv` |
| `tree.py` | Distances and a neighbor joining tree | `genomes.nwk` |
| `suspect.py` | Genomes whose mutations do not fit their date | prints |
| `report.py` | The findings as Markdown | `report.md` |

## Data

Every genome comes from NCBI GenBank through NCBI Virus and the Datasets
service, and is in the public domain. Nothing is redistributed here:
`download.py` fetches each record from NCBI by accession, and `.gitignore`
keeps the downloads out of the repository. The reference genome is
NC_045512.2, Wuhan-Hu-1, collected in December 2019.

## What it cannot say

The sample is ten genomes per lineage, taken from what NCBI returned
first. It is not a random sample of infections, so it cannot say how
common a lineage was, when a mutation first arose anywhere, or where it
spread. It can say what the lineages differ by, and that is the question
the project asks.
