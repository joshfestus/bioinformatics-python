# What ClinVar records about variants in HBB, fetched by region from
# NCBI's indexed ClinVar file without downloading all of it.
from collections import Counter
from pathlib import Path

import pysam

CLINVAR = (
    "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz"
)
HBB = ("11", 5_225_463, 5_227_071)  # GRCh38, zero based start
SAVED = Path("hbb_clinvar.vcf")

if not SAVED.exists():
    with (
        pysam.VariantFile(CLINVAR) as remote,
        pysam.VariantFile(SAVED, "w", header=remote.header) as out,
    ):
        for record in remote.fetch(*HBB):
            out.write(record)


def readable(values: tuple) -> str:
    """ClinVar's ("criteria_provided", "_multiple_submitters") as text."""
    return ",".join(values).replace("_", " ").strip()


records = list(pysam.VariantFile(SAVED))
print(f"{len(records)} ClinVar records in HBB")

significance = Counter(
    r.info.get("CLNSIG", ("none",))[0] for r in records
)
for label, count in significance.most_common(6):
    print(f"  {readable((label,)):44} {count:5}")

sickle = next(r for r in records if r.info.get("RS") == ("334",))
change = f"{sickle.ref}>{sickle.alts[0]}"
print(f"\nrs334 at chromosome {sickle.chrom}, {sickle.pos:,}, {change}")
print(f"  significance  {readable(sickle.info['CLNSIG'])}")
print(f"  review status {readable(sickle.info['CLNREVSTAT'])}")
print(f"  description   {sickle.info['CLNHGVS'][0]}")
print(f"  consequence   {sickle.info['MC'][0]}")
