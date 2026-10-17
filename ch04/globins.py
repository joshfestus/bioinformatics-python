# Three globin coding sequences from NCBI: human HBB (NM_000518.5),
# human HBA1 (NM_000558.5) and mouse Hbb-bt (NM_008220.5).
from fasta import gene_name, read_fasta
from seqtools import gc_content, reverse_complement, translate

genes = read_fasta("globin_cds.fasta")

for header, dna in genes.items():
    protein = translate(dna)
    print(gene_name(header))
    print(f"  {len(dna)} bases, GC {gc_content(dna):.1%}")
    print(f"  protein: {len(protein)} amino acids, {protein[:12]}...")
    print(f"  other strand starts {reverse_complement(dna)[:12]}...")
