# Compare the naive caller's single letter calls with bcftools'.
import pysam

with open("naive.tsv") as handle:
    naive = {int(line.split("\t")[0]) for line in handle}
bcftools = {
    record.pos
    for record in pysam.VariantFile("haploid.vcf.gz")
    if record.qual >= 30
    and record.info["DP"] >= 10
    and len(record.ref) == 1
    and len(record.alts[0]) == 1
}
print(
    f"naive {len(naive)}, bcftools {len(bcftools)}, "
    f"both {len(naive & bcftools)}"
)
print(f"only bcftools: {sorted(bcftools - naive)}")
extra = sorted(naive - bcftools)
if extra:
    print(f"only naive: {len(extra)}, the first at {extra[:5]}")
else:
    print("only naive: none")
