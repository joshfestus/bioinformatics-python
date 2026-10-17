# What InterPro knows about the families and domains of two proteins.
import json
import urllib.request

INTERPRO = (
    "https://www.ebi.ac.uk/interpro/api/entry/interpro/protein/uniprot"
)

for accession, name in (
    ("P68871", "beta globin"),
    ("P08100", "rhodopsin"),
):
    with urllib.request.urlopen(f"{INTERPRO}/{accession}/") as response:
        entries = json.load(response)["results"]
    print(f"{name} ({accession}): {len(entries)} InterPro entries")
    for entry in entries:
        meta = entry["metadata"]
        pieces = [
            piece
            for location in entry["proteins"][0][
                "entry_protein_locations"
            ]
            for piece in location["fragments"]
        ]
        start = min(piece["start"] for piece in pieces)
        end = max(piece["end"] for piece in pieces)
        name = meta["name"][:26]
        print(f"  {meta['accession']}  {start:3} to {end:3}  {name}")
