# The human globin genes in Ensembl, joined to NCBI's gene table from
# Chapter 5 on their shared NCBI GeneID.
import json
import time
import urllib.request

import pandas as pd

ENSEMBL = "https://rest.ensembl.org"
GENES = ["HBB", "HBD", "HBG1", "HBG2", "HBE1", "HBA1", "HBA2", "HBZ"]


def get_json(path: str) -> dict | list:
    request = urllib.request.Request(
        ENSEMBL + path, headers={"Content-Type": "application/json"}
    )
    with urllib.request.urlopen(request) as response:
        return json.load(response)


rows = []
for symbol in GENES:
    gene = get_json(f"/lookup/symbol/homo_sapiens/{symbol}")
    xrefs = get_json(f"/xrefs/id/{gene['id']}?all_levels=1")
    by_database = {x["dbname"]: x["primary_id"] for x in xrefs}
    rows.append(
        {
            "symbol": symbol,
            "ensembl": gene["id"],
            "GeneID": int(by_database["EntrezGene"]),
            "uniprot": by_database.get("Uniprot/SWISSPROT"),
            "ens_start": gene["start"],
            "ens_end": gene["end"],
        }
    )
    time.sleep(0.1)  # Ensembl allows 15 requests a second

ensembl = pd.DataFrame(rows)
ncbi = pd.read_csv(
    "../ch05/genes.tsv", sep="\t", dtype={"chromosome": str}
)
table = ensembl.merge(
    ncbi[["GeneID", "chromosome", "start", "end"]], on="GeneID"
)
table = table.rename(columns={"start": "ncbi_start", "end": "ncbi_end"})
print(table[["symbol", "ensembl", "GeneID", "uniprot"]].to_string())
print()
print(
    table[
        [
            "symbol",
            "chromosome",
            "ncbi_start",
            "ens_start",
            "ncbi_end",
            "ens_end",
        ]
    ].to_string()
)
