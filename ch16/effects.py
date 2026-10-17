# What each confident variant does: which gene it falls in and, for
# single letter changes in a coding sequence, the amino acid change.
import pandas as pd
import pysam
from Bio import SeqIO
from Bio.Seq import Seq

reference = SeqIO.read("../data/rel606.fasta", "fasta").seq
columns = [
    "seqid",
    "source",
    "type",
    "start",
    "end",
    "score",
    "strand",
    "phase",
    "attributes",
]
gff = pd.read_csv(
    "../ch15/rel606.gff3", sep="\t", comment="#", names=columns
)
cds = gff[gff["type"] == "CDS"].drop_duplicates(subset=["start", "end"])
cds = cds.assign(
    gene=cds["attributes"].str.extract(r"gene=([^;]+)", expand=False),
    product=cds["attributes"].str.extract(
        r"product=([^;]+)", expand=False
    ),
)


def effect(row, position: int, alt: str) -> str:
    """The amino acid change a single letter change makes in one CDS."""
    start, end = row.start - 1, row.end  # zero based, end excluded
    gene = reference[start:end]
    offset = position - 1 - start
    if row.strand == "-":
        gene = gene.reverse_complement()
        offset = end - position
        alt = str(Seq(alt).reverse_complement())
    mutated = gene[:offset] + alt + gene[offset + 1 :]
    codon = offset // 3
    before = gene[3 * codon : 3 * codon + 3].translate(table=11)
    after = mutated[3 * codon : 3 * codon + 3].translate(table=11)
    kind = (
        "synonymous"
        if before == after
        else ("nonsense" if after == "*" else "missense")
    )
    return f"{before}{codon + 1}{after} {kind}"


for record in pysam.VariantFile("haploid.vcf.gz"):
    if record.qual < 30 or record.info["DP"] < 10:
        continue
    alt = record.alts[0]
    hits = cds[(cds["start"] <= record.pos) & (cds["end"] >= record.pos)]
    if hits.empty:
        where = "between genes"
    else:
        row = hits.iloc[0]
        name = (
            row.gene if isinstance(row.gene, str) else row["product"][:24]
        )
        if "IS1" in row["product"]:
            name = "IS1 element"
        if "frameshift" in row["attributes"]:
            where = f"{name}: frameshifted gene, not translated"
        elif len(record.ref) == 1 and len(alt) == 1:
            where = f"{name}: {effect(row, record.pos, alt)}"
        else:
            where = f"{name}: insertion or deletion"
    change = f"{record.ref[:4]}>{alt[:4]}"
    print(f"{record.pos:9,}  {change:10} {where}")
