"""Choose a sample of SARS-CoV-2 genomes, a few from each major lineage.

Asks NCBI Virus, through the NCBI Datasets API, for complete genomes of
each lineage and keeps ten with a full collection date. NCBI adds genomes
every week, so a run today chooses different ones from the book's; the
book's choice is in samples.tsv, which download.py reads.

uv run choose.py            writes my_samples.tsv
"""

import csv
import json
import random
import time
import urllib.error
import urllib.parse
import urllib.request

API = "https://api.ncbi.nlm.nih.gov/datasets/v2/virus/taxon/2697049"

# Pango lineages in the order they spread, with the name the World Health
# Organization gave the variants of concern.
LINEAGES = {
    "A": "",
    "B": "",
    "B.1": "",
    "B.1.1.7": "Alpha",
    "B.1.351": "Beta",
    "P.1": "Gamma",
    "B.1.617.2": "Delta",
    "BA.1": "Omicron",
    "BA.2": "Omicron",
    "BA.5.2.1": "Omicron",
    "XBB.1.5": "Omicron",
    "JN.1": "Omicron",
    "KP.3": "Omicron",
}
PER_LINEAGE = 10

# The service runs out of time on Alpha, which has hundreds of thousands
# of genomes, unless the question is narrowed to one country.
NARROW = {"B.1.1.7": "Nigeria"}
ATTEMPTS = 4  # the service fails now and then, and works on a second ask


def reports(lineage, how_many=200):
    filters = {
        "filter.pangolin_classification": lineage,
        "filter.complete_only": "true",
        "page_size": how_many,
    }
    if lineage in NARROW:
        filters["filter.geo_location"] = NARROW[lineage]
    query = urllib.parse.urlencode(filters)
    for attempt in range(1, ATTEMPTS + 1):
        try:
            with urllib.request.urlopen(f"{API}/dataset_report?{query}",
                                        timeout=300) as response:
                return json.load(response).get("reports", [])
        except urllib.error.HTTPError as problem:
            if problem.code < 500 or attempt == ATTEMPTS:
                raise
            print(f"  {lineage}: {problem.code} from the service, "
                  f"asking again ({attempt} of {ATTEMPTS - 1})")
            time.sleep(5 * attempt)


def main():
    random.seed(1)
    rows = []
    for lineage, who in LINEAGES.items():
        usable = []
        for r in reports(lineage):
            date = r.get("isolate", {}).get("collection_date", "")
            human = r.get("host", {}).get("tax_id") == 9606
            if human and len(date) == 10 and r.get("length", 0) >= 29_000:
                usable.append({
                    "accession": r["accession"],
                    "lineage": lineage,
                    "who": who,
                    "date": date,
                    "country": r.get("location", {})
                    .get("geographic_location", "").split(":")[0],
                })
        chosen = random.sample(usable, min(PER_LINEAGE, len(usable)))
        print(f"{lineage:10} {len(usable):4} usable, "
              f"{len(chosen)} chosen")
        rows += sorted(chosen, key=lambda row: row["date"])
    with open("my_samples.tsv", "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=rows[0], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
