"""Tests for variants.py. Run with: uv run pytest"""

from pathlib import Path

import pytest
from Bio.Seq import Seq

from variants import amino_acid_change, annotate

# A made-up gene of three codons, ATG GAT TAA (Met, Asp, stop), placed at
# positions 101 to 109 of a genome, counting from 1.
GENE = Seq("ATGGATTAA")


def test_missense_on_the_plus_strand():
    # position 104 is the G of GAT; GAT to AAT turns Asp into Asn
    assert (
        amino_acid_change(GENE, "+", 101, 109, 104, "A") == "D2N missense"
    )


def test_synonymous_change():
    # GAT to GAC still codes for Asp
    assert amino_acid_change(GENE, "+", 101, 109, 106, "C") == (
        "D2D synonymous"
    )


def test_nonsense_change():
    # ATG TGG TAA (Met, Trp, stop) at positions 1 to 9: changing the last
    # G of TGG, position 6, to A makes TGA, a stop codon
    gene = Seq("ATGTGGTAA")
    assert amino_acid_change(gene, "+", 1, 9, 6, "A") == "W2* nonsense"


def test_minus_strand_gene():
    # The same made-up gene written on the minus strand: the genome holds
    # its reverse complement, TTAATCCAT, at positions 101 to 109.
    genome_side = GENE.reverse_complement()
    # The gene's G of GAT (gene position 4) sits at genome position 106,
    # where the plus strand shows its complement, C. G>A on the gene is
    # C>T on the genome.
    assert amino_acid_change(genome_side, "-", 101, 109, 106, "T") == (
        "D2N missense"
    )


DATA = Path("../data/rel606.fasta")


@pytest.mark.skipif(not DATA.exists(), reason="needs data/fetch.py ch15")
def test_pykf_in_the_real_genome():
    rows = annotate(
        "../ch16/haploid.vcf.gz",
        "../ch15/rel606.gff3",
        str(DATA),
        min_qual=30,
        min_depth=10,
    )
    effects = {row["position"]: row["effect"] for row in rows}
    assert len(effects) == 25
    assert effects[1_733_343] == "pykF: D127N missense"
