# Open reading frames in 8 kb of the E. coli K-12 genome around the lac
# genes (NC_000913.3, positions 360,000 to 368,000), on all six frames.
from Bio import SeqIO
from Bio.Seq import Seq

MIN_LENGTH = 100  # amino acids


def find_orfs(dna: Seq, min_length: int = MIN_LENGTH) -> list[tuple]:
    """(start, end, strand, length) of each ORF, counted from 1."""
    orfs = []
    for strand, sequence in ((+1, dna), (-1, dna.reverse_complement())):
        for frame in range(3):
            usable = (len(sequence) - frame) // 3 * 3
            protein = str(sequence[frame : frame + usable].translate(11))
            position = 0
            while (start := protein.find("M", position)) != -1:
                stop = protein.find("*", start)
                if stop == -1:
                    break
                if stop - start >= min_length:
                    first = frame + 3 * start  # zero based, on `sequence`
                    last = (
                        frame + 3 * stop + 3
                    )  # just past the stop codon
                    if strand == -1:
                        first, last = len(dna) - last, len(dna) - first
                    orfs.append((first + 1, last, strand, stop - start))
                position = stop + 1
    return sorted(orfs)


if __name__ == "__main__":
    region = SeqIO.read("lac_region.gb", "genbank")
    print(f"{region.id}, {len(region)} bp from E. coli K-12 MG1655")

    print("\nOpen reading frames of at least 100 amino acids")
    for start, end, strand, length in find_orfs(region.seq):
        print(f"  {start:5} to {end:5}  strand {strand:+}  {length:5} aa")

    print("\nGenes NCBI annotates in the same region")
    for feature in region.features:
        if feature.type == "CDS":
            start, end = feature.location.start + 1, feature.location.end
            strand = feature.location.strand
            name = feature.qualifiers["gene"][0]
            print(f"  {start:5} to {end:5}  strand {strand:+}  {name}")
