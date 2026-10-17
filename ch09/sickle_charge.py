# Normal and sickle beta globin compared by what their sequences predict.
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis

beta = str(SeqIO.read("hbb_human.fasta", "fasta").seq)
mature = beta[1:]  # the cell removes the first methionine
sickle = mature[:5] + "V" + mature[6:]  # glutamic acid 6 becomes valine

for name, protein in (("normal", mature), ("sickle", sickle)):
    analysis = ProteinAnalysis(protein)
    print(name)
    print(f"  length              {len(protein):8}")
    print(f"  molecular weight    {analysis.molecular_weight():8.1f} Da")
    print(f"  isoelectric point   {analysis.isoelectric_point():8.2f}")
    print(f"  charge at pH 7.4    {analysis.charge_at_pH(7.4):8.2f}")
    print(f"  GRAVY               {analysis.gravy():8.3f}")
