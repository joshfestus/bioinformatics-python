"""Solutions to the exercises of Project Two, Chapter 24.

Run from this folder: uv run ch24.py
Reads the project's files in ../projects/resistance. Exercises 24.1,
24.6 and 24.7 need BLAST+; 24.6 also downloads the catalog's DNA.
"""

import csv
import shutil
import subprocess
import sys
import urllib.error
import urllib.request
from collections import Counter
from pathlib import Path

from Bio import Align, SeqIO
from Bio.Align import substitution_matrices

WORK = Path("../projects/resistance").resolve()
AMR = ("https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance"
       "/AMRFinderPlus/database/latest")
PKPHS3 = "NC_016839.1"
MOBILE = ("transposase", "integrase", "recombinase", "resolvase")


def tool(name):
    beside = Path(sys.executable).parent / name
    return str(beside) if beside.exists() else shutil.which(name)


def read(name):
    with open(WORK / name) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def exercise_24_1():
    print("24.1  The same search at 80 percent identity over 60")
    print("      percent coverage.")
    rows = read("resistance.tsv")
    tight = {(row["gene"], row["locus_tag"]) for row in rows
             if row["subtype"] == "AMR"}
    if not tool("blastp"):
        print("      BLAST+ is not on the PATH; see Appendix A.")
        return
    done = subprocess.run(
        [sys.executable, "filter.py", "--min-identity", "80",
         "--min-coverage", "60"], cwd=WORK, capture_output=True,
        text=True)
    if done.returncode != 0:
        print(f"      filter.py failed: {done.stderr.strip()[:120]}")
        return
    loose = {(row["gene"], row["locus_tag"]) for row in read(
        "resistance.tsv") if row["subtype"] == "AMR"}
    # put the strict thresholds back before anything else reads the file
    subprocess.run([sys.executable, "filter.py"], cwd=WORK,
                   capture_output=True, text=True)
    added = sorted(name for name, _ in loose - tight)
    print(f"       {len(tight)} genes at 90/90, {len(loose)} at 80/60")
    print(f"       the ones that join: {', '.join(added)}")
    print("      None of the three should be reported as a finding.")
    print("      acrF is a full length protein 86 percent identical")
    print("      to the catalog's, which makes it a member of the")
    print("      same efflux family rather than that gene, and every")
    print("      Klebsiella has one. sul1 is a perfect match over 60")
    print("      percent of the reference because the copy here is")
    print("      truncated, and NCBI's own annotation says so. The")
    print("      aac(3)-IId match is 100 percent identical over 81")
    print("      percent, a protein 30 residues short of the")
    print("      reference, and whether a shortened enzyme still")
    print("      works cannot be settled from the sequence.")
    print("      Each would need a sentence of explanation in a")
    print("      report, which is the test: a finding you cannot")
    print("      state in one line is not yet a finding.")


def exercise_24_2():
    rows = read("resistance.tsv")
    kinds = Counter(row["subtype"] for row in rows)
    print("\n24.2  Everything the search found, by kind:")
    for subtype, count in kinds.most_common():
        print(f"       {subtype or 'unclassified':14} {count:>3}")
    names = {"NC_016845.1": "chromosome", "NC_016838.1": "pKPHS1",
             "NC_016846.1": "pKPHS2", "NC_016839.1": "pKPHS3"}
    for subtype in ("METAL", "VIRULENCE", "BIOCIDE"):
        found = [row for row in rows if row["subtype"] == subtype]
        if not found:
            continue
        print(f"\n       {subtype.lower()} hits:")
        for row in sorted(found, key=lambda r: int(r["start"])):
            where = names.get(row["replicon"], row["replicon"])
            print(f"         {row['gene']:12} {where:11} "
                  f"{int(row['start']):>9,}  {row['product'][:28]}")
    mercury = {names.get(row["replicon"], row["replicon"])
               for row in rows if row["subtype"] == "METAL"
               and row["gene"].lower().startswith("mer")}
    print(f"\n      The four mercury genes, the mer cluster, sit")
    print(f"      together on {', '.join(sorted(mercury))}. The fifth")
    print("      metal hit, fieF, is a chromosomal cation pump that")
    print("      every Klebsiella has and has nothing to do with")
    print("      mercury or with this plasmid.")
    print("      That is the same plasmid as blaKPC-2, and it is not")
    print("      a coincidence. Mercury resistance travels on")
    print("      transposons of the Tn21 family, which carry")
    print("      antibiotic resistance genes along with it, so a")
    print("      plasmid selected by mercury in the environment can")
    print("      carry a carbapenemase into a hospital. The")
    print("      virulence hits are iron scavenging genes, which are")
    print("      about causing disease rather than surviving drugs.")


def exercise_24_3():
    print("\n24.3  Running the project on a second isolate.")
    print("      Replace the seven accessions in download.py's")
    print("      REPLICONS with another complete genome's, from NCBI")
    print("      Datasets or the assembly page, and rerun:\n")
    print("        uv run download.py")
    print("        snakemake --cores 2 all\n")
    print("      Good comparisons are a carbapenemase producer from")
    print("      a different country, or a carriage isolate from a")
    print("      healthy person, which should carry the species'")
    print("      intrinsic blaSHV and fosA and little else.")
    print("      What to compare is not the gene lists but where the")
    print("      genes sit. Two isolates sharing blaKPC-2 on similar")
    print("      plasmids is an outbreak question; two sharing only")
    print("      chromosomal genes is two ordinary Klebsiella.")
    print("      Expect the counts to differ by more than the")
    print("      biology does: a draft genome in contigs cannot say")
    print("      whether a gene is on a plasmid, and the catalog")
    print("      version changes what is findable at all, which is")
    print("      why the report prints it.")


def exercise_24_4():
    record = next(r for r in SeqIO.parse(WORK / "genome.gb", "genbank")
                  if r.id == PKPHS3)
    resistance = {int(row["start"]): row["gene"]
                  for row in read("profile.tsv")
                  if row["accession"] == PKPHS3}
    print("\n24.4  Every gene of the island, 29,000 to 56,000 on "
          "pKPHS3:\n")
    print(f"       {'start':>7} {'kind':12} product")
    tally = Counter()
    for feature in record.features:
        if feature.type != "CDS":
            continue
        start = int(feature.location.start) + 1
        if not 29_000 <= start <= 56_000:
            continue
        product = feature.qualifiers.get("product", ["?"])[0]
        if start in resistance:
            kind = "resistance"
        elif any(word in product.lower() for word in MOBILE):
            kind = "mobility"
        else:
            kind = "neither"
        tally[kind] += 1
        label = resistance.get(start, product)
        print(f"       {start:>7,} {kind:12} {label[:46]}")
    print(f"\n      {tally['resistance']} resistance genes, "
          f"{tally['mobility']} mobility genes, {tally['neither']} "
          "neither.")
    print("      The neither pile is mostly hypothetical proteins,")
    print("      which is what an island looks like: a stretch")
    print("      assembled from pieces of other genomes, carrying")
    print("      things nobody has characterized along with the")
    print("      genes that made it worth keeping.")


def exercise_24_5():
    print("\n24.5  Checking gyrA's residues by hand.")
    proteins = SeqIO.to_dict(SeqIO.parse(WORK / "proteins.fasta",
                                         "fasta"))
    mutations = read("mutations.tsv")
    gyra = next((row for row in mutations
                 if row["mutation"].startswith("gyrA")), None)
    if gyra is None:
        print("      No gyrA mutation in mutations.tsv.")
        return
    isolate = proteins[gyra["locus_tag"]]
    # the catalog holds 26 gyrA proteins with different numbering, so
    # take the one the Klebsiella rules are written against, as
    # point.py does, rather than the first that matches the name
    wanted = None
    with open(WORK / "amr_mutations.tsv") as handle:
        header = handle.readline().lstrip("#").rstrip("\n").split("\t")
        for line in handle:
            row = dict(zip(header, line.rstrip("\n").split("\t")))
            if (row["taxgroup"] == "Klebsiella_pneumoniae"
                    and row["standard_mutation_symbol"].startswith(
                        "gyrA")):
                wanted = row["accession_version"]
                break
    reference = None
    for record in SeqIO.parse(WORK / "amr_proteins.fasta", "fasta"):
        if record.id.split("|")[0] == wanted:
            reference = record
            break
    if reference is None:
        print("      No gyrA reference in the catalog.")
        return
    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score, aligner.extend_gap_score = -11, -1
    alignment = aligner.align(isolate.seq, reference.seq)[0]
    columns = {}
    for (mine, mine_end), (theirs, theirs_end) in zip(*alignment.aligned):
        for offset in range(mine_end - mine):
            columns[theirs + offset + 1] = mine + offset
    print(f"       isolate {isolate.id}, {len(isolate.seq)} residues")
    print(f"       catalog reference {reference.id.split('|')[0]}, "
          f"{len(reference.seq)} residues")
    for position in (83, 87):
        index = columns.get(position)
        mine = isolate.seq[index] if index is not None else "?"
        theirs = reference.seq[position - 1]
        print(f"       position {position}: reference {theirs}, "
              f"isolate {mine}")
    print(f"      point.py reported {gyra['mutation']}, and reading")
    print("      the alignment by hand gives the same answer: the")
    print("      reference has serine at 83 and this isolate has")
    print("      isoleucine. Position 87 is aspartic acid in both,")
    print("      unchanged, which matters as much: a program that")
    print("      found a mutation at every position it checked")
    print("      would be broken.")
    print("      One trap on the way. The catalog holds 26 proteins")
    print("      called gyrA, and they do not share a numbering.")
    print("      Align against the wrong one and position 83 is a")
    print("      valine that was never mutated. The reference to use")
    print("      is the one the Klebsiella rules name, which is why")
    print("      point.py keys on the accession rather than the gene")
    print("      name, and why this answer does the same.")


def exercise_24_6():
    print("\n24.6  Searching the catalog's DNA with blastn.")
    cds = WORK / "amr_cds.fasta"
    if not cds.exists():
        try:
            request = urllib.request.Request(
                f"{AMR}/AMR_CDS.fa", headers={"User-Agent": "python"})
            with urllib.request.urlopen(request, timeout=600) as reply:
                cds.write_bytes(reply.read())
        except (urllib.error.URLError, TimeoutError) as problem:
            print(f"      Could not download AMR_CDS.fa ({problem}).")
            return
    if not tool("makeblastdb") or not tool("blastn"):
        print("      BLAST+ is not on the PATH; see Appendix A.")
        return
    subprocess.run([tool("makeblastdb"), "-in", "amr_cds.fasta",
                    "-dbtype", "nucl", "-out", "amr_cds"], cwd=WORK,
                   capture_output=True, text=True, check=True)
    nucleotide = WORK / "genome.fasta"
    if not nucleotide.exists():
        records = list(SeqIO.parse(WORK / "genome.gb", "genbank"))
        SeqIO.write(records, nucleotide, "fasta")
    done = subprocess.run(
        [tool("blastn"), "-query", "genome.fasta", "-db", "amr_cds",
         "-outfmt", "6 qseqid qstart qend sseqid pident length slen",
         "-evalue", "1e-20", "-perc_identity", "90"], cwd=WORK,
        capture_output=True, text=True, check=True)
    names, places = set(), set()
    for line in done.stdout.splitlines():
        query, start, end, subject, pident, length, slen = \
            line.split("\t")
        if int(length) < 0.9 * int(slen):
            continue
        fields = subject.split("|")
        names.add(fields[5] if len(fields) > 5 else subject)
        # round the position so the same locus hit by several alleles
        # of one gene counts once
        places.add((query, round(int(start) / 500)))
    protein_genes = {row["gene"] for row in read("resistance.tsv")
                     if row["subtype"] == "AMR"}
    print(f"       blastn reports {len(names)} different gene names")
    print(f"       at 90 percent identity over 90 percent of their")
    print(f"       length, but they sit at only {len(places)} places "
          "in the genome.")
    print(f"       the protein search named {len(protein_genes)} genes.")
    families = Counter(name.rstrip("0123456789-") for name in names)
    print("\n       the names that pile up on one locus:")
    for family, count in families.most_common(3):
        if count > 1:
            members = sorted(n for n in names
                             if n.rstrip("0123456789-") == family)
            print(f"         {family}: {count} alleles, "
                  f"{', '.join(members[:6])}")
    print("      That is the first lesson, and it is not the one the")
    print("      exercise asks for. The catalog carries every known")
    print("      allele of a gene, and at the DNA level they are")
    print("      almost identical, so one aadA cassette in this")
    print("      genome matches a dozen aadA entries at once.")
    print("      Count the places, not the names, or the profile")
    print("      doubles for no biological reason.")
    print("      The answer to what the exercise asks is still yes:")
    print("      a DNA search is the one that finds a gene whose")
    print("      reading frame has been broken. A frameshift or an")
    print("      early stop destroys the protein, so a protein")
    print("      search sees a fragment or nothing, while the DNA is")
    print("      still almost identical to the reference. Chapter")
    print("      16 met the same situation from the other side,")
    print("      where four variants sat in genes the annotation")
    print("      itself calls frameshifted.")
    print("      Which search you want depends on the question. For")
    print("      what the bacterium can do, use the protein: a")
    print("      broken gene confers nothing. For where a gene has")
    print("      been, use the DNA, because a broken copy is still a")
    print("      record of its arrival.")
    for name in WORK.glob("amr_cds.n*"):
        name.unlink()


def exercise_24_7():
    print("\n24.7  Confirming that blaNDM-1 is absent.")
    rows = [row for row in read("resistance.tsv")
            if "NDM" in row["gene"].upper()]
    print(f"       hits to any NDM gene in the profile: {len(rows)}")
    if not tool("blastp"):
        print("      BLAST+ is not on the PATH to search further.")
    else:
        reference = None
        for record in SeqIO.parse(WORK / "amr_proteins.fasta", "fasta"):
            if "NDM-1" in record.id:
                reference = record
                break
        if reference is not None:
            query = WORK / "ndm1.fasta"
            SeqIO.write([reference], query, "fasta")
            subprocess.run([tool("makeblastdb"), "-in",
                            "proteins.fasta", "-dbtype", "prot",
                            "-out", "isolate"], cwd=WORK,
                           capture_output=True, text=True, check=True)
            done = subprocess.run(
                [tool("blastp"), "-query", "ndm1.fasta", "-db",
                 "isolate", "-outfmt", "6 sseqid pident length evalue",
                 "-evalue", "10"], cwd=WORK, capture_output=True,
                text=True, check=True)
            lines = done.stdout.splitlines()
            total = sum(1 for _ in SeqIO.parse(
                WORK / "proteins.fasta", "fasta"))
            print(f"       blastp with blaNDM-1 against all "
                  f"{total:,} proteins")
            print(f"       of this isolate, at E up to 10: "
                  f"{len(lines)} hits")
            for line in lines[:3]:
                subject, pident, length, evalue = line.split("\t")
                print(f"         {subject} {pident}% over {length}, "
                      f"E {evalue}")
            query.unlink()
            for name in WORK.glob("isolate.p*"):
                name.unlink()
    print("      An absence is harder to prove than a presence, and")
    print("      three things have to hold before this one counts.")
    print("      The genome must be complete, which this one is:")
    print("      seven circular replicons, nothing in contigs, so")
    print("      there is no unassembled corner for a gene to hide")
    print("      in. The search must be able to find the gene if it")
    print("      were there, which a blastp against every protein")
    print("      at a generous threshold establishes. And the")
    print("      catalog must contain the gene, or you are testing")
    print("      your database rather than the isolate.")
    print("      With a draft genome none of that holds, and the")
    print("      honest sentence is that the gene was not found")
    print("      rather than that it is not there.")


if __name__ == "__main__":
    exercise_24_1()
    exercise_24_2()
    exercise_24_3()
    exercise_24_4()
    exercise_24_5()
    exercise_24_6()
    exercise_24_7()
