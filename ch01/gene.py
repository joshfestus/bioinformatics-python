# The coding sequence of human beta globin (HBB), NCBI NM_000518.5
gene = open("hbb.txt").read().strip()

print("Length:", len(gene), "bases")
for base in "ACGT":
    print(base, gene.count(base))

gc = (gene.count("G") + gene.count("C")) / len(gene)
print(f"GC content: {gc:.1%}")

codons = [gene[i : i + 3] for i in range(0, len(gene), 3)]
print("Codons:", len(codons))
print("First eight:", " ".join(codons[:8]))
print("Codon 7:", codons[6])

sickle = gene[:19] + "T" + gene[20:]
print("Sickle codon 7:", sickle[18:21])
