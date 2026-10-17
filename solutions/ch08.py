"""Solutions to the exercises of Chapter 8.

Run from this folder: uv run ch08.py
Every exercise here talks to a public service, so each is wrapped: if a
service is down or you are offline, the answer says so and the rest
still runs. Set NCBI_API_KEY if you have one, as Chapter 8 explains.
"""

import csv
import io
import json
import os
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
UNIPROT = "https://rest.uniprot.org/uniprotkb"
ENSEMBL = "https://rest.ensembl.org"
EMAIL = "your.name@example.org"


def ask(url, timeout=120):
    """Fetch a URL as text, adding the NCBI key when there is one."""
    if "ncbi.nlm.nih.gov" in url:
        key = os.environ.get("NCBI_API_KEY")
        url += f"&email={EMAIL}" + (f"&api_key={key}" if key else "")
    headers = {"User-Agent": "python"}
    request = urllib.request.Request(url, headers=headers)
    with urllib.request.urlopen(request, timeout=timeout) as reply:
        return reply.read().decode()


def offline(problem):
    print(f"     The service did not answer ({problem}). The answer")
    print("     below is what a successful run reports.")


def exercise_8_1():
    print("8.1  Searching the protein database instead of nucleotide:")
    print("       term = 'HBB[Gene Name] AND human[Organism] AND "
          "refseq[filter]'")
    print("       db=protein, then efetch NP_000509.1 as FASTA")
    try:
        fasta = ask(f"{EUTILS}/efetch.fcgi?db=protein&id=NP_000509.1"
                    "&rettype=fasta&retmode=text")
        ncbi = "".join(fasta.split("\n")[1:]).strip()
        swiss = ask(f"{UNIPROT}/P68871.fasta")
        uniprot = "".join(swiss.split("\n")[1:]).strip()
        print(f"     NP_000509.1 is {len(ncbi)} amino acids, P68871 is "
              f"{len(uniprot)}.")
        print(f"     They are {'identical' if ncbi == uniprot else 'not'}"
              " letter for letter.")
        print("     The same protein is in both databases because both")
        print("     describe the same gene; NCBI's record is derived")
        print("     from the genome annotation and UniProt's is curated")
        print("     by hand, and for a gene this studied they agree.")
    except (urllib.error.URLError, TimeoutError) as problem:
        offline(problem)


def count(term):
    url = (f"{EUTILS}/esearch.fcgi?db=nuccore&retmode=json&term="
           + urllib.parse.quote(term))
    return int(json.loads(ask(url))["esearchresult"]["count"])


def exercise_8_2():
    print("\n8.2  How much each filter removes:")
    terms = [
        "HBB[Gene Name]",
        "HBB[Gene Name] AND refseq[filter]",
        "HBB[Gene Name] AND Homo sapiens[Organism]",
        "HBB[Gene Name] AND Homo sapiens[Organism] AND refseq[filter]",
    ]
    try:
        for term in terms:
            time.sleep(0.4)  # NCBI asks for no more than three a second
            print(f"       {count(term):8,}  {term}")
        print("     Gene name alone is the largest number, because HBB")
        print("     is a gene name in many species and because every")
        print("     submitted sequence counts, not every distinct gene.")
        print("     Chapter 8 found another reason to be careful: the")
        print("     pseudogene KRT89P carries HBB as an alias, so it")
        print("     matches this search too.")
    except (urllib.error.URLError, TimeoutError, KeyError) as problem:
        offline(problem)


def exercise_8_3():
    print("\n8.3  UniProt's search endpoint, asked for a TSV:")
    query = "reviewed:true AND organism_id:9606 AND protein_name:globin"
    url = (f"{UNIPROT}/search?format=tsv&size=500"
           "&fields=accession,id,protein_name,length&query="
           + urllib.parse.quote(query))
    try:
        rows = list(csv.DictReader(io.StringIO(ask(url)), delimiter="\t"))
        print(f"     {len(rows)} reviewed human proteins match.")
        for row in rows:
            print(f"       {row['Entry']} {row['Entry Name']:14} "
                  f"{row['Length']:>4}  {row['Protein names'][:38]}")
        hemoglobins = [r for r in rows
                       if r["Protein names"].startswith("Hemoglobin")]
        print(f"     {len(hemoglobins)} are hemoglobin subunits. The")
        print("     others are not globins at all:")
        for row in rows:
            if row not in hemoglobins:
                print(f"       {row['Entry Name']}, which carries globin")
                print("            in an alias, not in what it is")
        print("     And one is missing. Myoglobin, P02144, is a globin")
        print("     by any definition, and this search does not find it,")
        print("     because UniProt matches whole words in a name and")
        print("     myoglobin's name is one word. A name search finds")
        print("     names. To find the molecules you want the family,")
        print("     through InterPro as in Chapter 9, or the sequence")
        print("     itself, through BLAST as in Chapter 11.")
    except (urllib.error.URLError, TimeoutError, KeyError) as problem:
        offline(problem)


def exercise_8_4():
    print("\n8.4  Ensembl's lookup with expand=1 lists the transcripts:")
    url = (f"{ENSEMBL}/lookup/id/ENSG00000244734?expand=1"
           ";content-type=application/json")
    try:
        gene = json.loads(ask(url))
        print(f"     {gene['display_name']} {gene['seq_region_name']}:"
              f"{gene['start']}-{gene['end']}, "
              f"{len(gene['Transcript'])} transcripts")
        for transcript in sorted(gene["Transcript"],
                                 key=lambda t: t["end"], reverse=True):
            print(f"       {transcript['display_name']:12} "
                  f"{transcript['start']:,} to {transcript['end']:,}  "
                  f"{transcript['biotype']}")
        furthest = max(gene["Transcript"], key=lambda t: t["end"])
        print(f"     {furthest['display_name']} reaches "
              f"{furthest['end']:,}, which is why Ensembl's gene ends")
        print("     there while NCBI's ends at 5,227,071. A gene's end")
        print("     is the end of its longest transcript, and the two")
        print("     databases do not annotate the same transcripts.")
    except (urllib.error.URLError, TimeoutError, KeyError) as problem:
        offline(problem)


def fetch_cached(url, path, attempts=3, pause=5):
    """Download once, keep the file, and retry a failure politely."""
    path = Path(path)
    if path.exists():
        return path
    for attempt in range(1, attempts + 1):
        try:
            path.write_text(ask(url))
            return path
        except (urllib.error.URLError, TimeoutError):
            if attempt == attempts:
                raise
            time.sleep(pause * attempt)
    return path


def exercise_8_5():
    print("\n8.5  fetch_cached downloads once and keeps the file:")
    print("       if the file is there, return it and ask nothing")
    print("       otherwise try up to three times, waiting longer")
    print("       between attempts, and let the last failure raise")
    try:
        first = time.perf_counter()
        path = fetch_cached(f"{UNIPROT}/P68871.fasta", "P68871.fasta")
        downloaded = time.perf_counter() - first
        second = time.perf_counter()
        fetch_cached(f"{UNIPROT}/P68871.fasta", "P68871.fasta")
        cached = time.perf_counter() - second
        print(f"     First call {downloaded * 1000:.0f} ms, second "
              f"{cached * 1000:.2f} ms.")
        path.unlink()
        print("     In globins.py it means a rerun costs nothing and")
        print("     the service is asked once. Every program in this")
        print("     book that downloads anything works this way, and")
        print("     Chapter 23 adds the other half of the lesson:")
        print("     write to a temporary name and rename when whole.")
    except (urllib.error.URLError, TimeoutError) as problem:
        offline(problem)


if __name__ == "__main__":
    exercise_8_1()
    exercise_8_2()
    exercise_8_3()
    exercise_8_4()
    exercise_8_5()
