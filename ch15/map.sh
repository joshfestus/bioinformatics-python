# Map the E. coli reads to the ancestral REL606 genome, then sort and
# index the result. Run from the ch15 folder with: bash map.sh
REFERENCE=../data/rel606.fasta
READS_1=../data/SRR2584863_1.fastq.gz
READS_2=../data/SRR2584863_2.fastq.gz

minimap2 -ax sr -t 2 $REFERENCE $READS_1 $READS_2 2> map.log \
    | samtools sort -o reads.bam -
samtools index reads.bam
echo "mapped: $(samtools view -c -F 2308 reads.bam) reads"
