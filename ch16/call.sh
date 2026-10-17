# Call variants in the mapped reads of Chapter 15 with bcftools: once with
# its default of two copies of each chromosome, once with one.
REFERENCE=../data/rel606.fasta
READS=../ch15/reads.bam

bcftools mpileup -Ou -f $REFERENCE $READS 2> /dev/null \
    | bcftools call -mv -Oz -o calls.vcf.gz
bcftools mpileup -Ou -f $REFERENCE $READS 2> /dev/null \
    | bcftools call --ploidy 1 -mv -Oz -o haploid.vcf.gz
bcftools index -f calls.vcf.gz
bcftools index -f haploid.vcf.gz
echo "two copies: $(bcftools view -H calls.vcf.gz | wc -l) variants"
echo "one copy:   $(bcftools view -H haploid.vcf.gz | wc -l) variants"
