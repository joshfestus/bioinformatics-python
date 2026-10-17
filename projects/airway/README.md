# What a steroid does to airway smooth muscle

One published RNA-seq experiment carried from a table of counts to a
statement about biology, and checked against the paper it came from.

This is Project Three of *Bioinformatics with Python for Beginners*
(Chapter 25). A full run takes under a minute once the data is downloaded.

## What it finds

Four donors' primary airway smooth muscle cells, treated with
dexamethasone or left untreated (GSE52778, Himes and others, 2014):

- 21,833 genes tested, 4,544 changed at an adjusted p below 0.05, and
  1,171 of those by at least twofold, 590 up and 581 down.
- Among the genes that went up: muscle processes, response to hormone, and
  the suppression of cell movement. "Cellular response to glucocorticoid
  stimulus" is enriched sevenfold, which is the analysis recovering the
  drug from the biology without being told it.
- Among the genes that went down: cell shape, cell projections and G
  protein coupled receptor signaling.
- All seven genes the paper names, DUSP1, KLF15, PER1, TSC22D3, C7,
  CCDC69 and CRISPLD2, are found here, all increased, four of them in the
  top thirteen genes of twenty-one thousand.

## Running it

```
uv run download.py              # counts, sample sheet and gene sets, 50 MB
snakemake --cores 2 all         # every step, to report.md
```

## The steps

| File | What it does | Writes |
|---|---|---|
| `download.py` | Fetches the experiment and the Gene Ontology | five files |
| `samples.py` | GEO's series matrix into a sample table | `samples.tsv` |
| `inspect_samples.py` | Library sizes, correlations, components | `qc.tsv` |
| `de.py` | Differential expression, donors paired | `results.tsv` |
| `enrich.py` | Which processes those genes belong to | `enrichment.tsv` |
| `compare.py` | The genes the original paper named | `comparison.tsv` |
| `report.py` | The findings as Markdown | `report.md` |

## Choices worth knowing about

The albuterol samples are left out. Chapter 20 found evidence that in two
donors those labels are swapped, and nothing here needs them.

The design is `~cell_line + treatment`, so each donor is their own
control. It matters: the largest axis of variation in these eight samples
is the donor, not the drug.

The enrichment background is the tested genes that carry an annotation,
14,053 of them, not every gene in the genome. Annotations are propagated
up the ontology, so a gene annotated to a narrow term counts for the wider
terms that contain it.

## What it cannot say

Enrichment names hypotheses, not mechanisms. A GO term records what a gene
was seen to do somewhere, often in another cell type, and some annotations
are computational predictions. The headings are a place to start reading
the gene lists, not a result in themselves.
