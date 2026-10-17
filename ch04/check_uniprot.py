# Do our translations match the proteins UniProt lists (Chapter 10)?
from fasta import gene_name, read_fasta
from seqtools import translate

genes = read_fasta("globin_cds.fasta")
proteins = read_fasta("../ch10/globins.fasta")

expected = {"HBB": "P68871", "HBA1": "P69905", "Hbb-bt": "P02088"}
uniprot = {header.split("|")[1]: seq for header, seq in proteins.items()}

for header, dna in genes.items():
    name = gene_name(header)
    ours, theirs = translate(dna), uniprot[expected[name]]
    print(f"{name:7} against {expected[name]}: ", end="")
    if ours == theirs:
        print("identical")
        continue
    print("differs at")
    for position, (a, b) in enumerate(zip(ours, theirs), start=1):
        if a != b:
            print(f"  position {position:3}: ours {a}, UniProt {b}")
