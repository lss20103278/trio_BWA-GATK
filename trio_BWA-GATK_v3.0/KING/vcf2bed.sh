python vcf_filter.py sample.vcf test6.vcf
java -Xmx15g  -jar /home/wei/bin/GenomeAnalysisTK.jar -R /home/wei/Documents/gatk/library/ucsc.hg19.fasta -T SelectVariants -V test6.vcf -env -ef -o test6-clean.vcf
/DATA/sslyu/soft/gatk-4.0.8.1/gatk SelectVariants -R /SSD750/PB1/db1/Homo/refseq/hg19.fa -V test.vcf --exclude-non-variants --exclude-filtered -O test6-clean.vcf
bcftools view --max-alleles 2 test6-clean.vcf > test7.vcf
./plink2 --vcf test7.vcf --make-bed --out ex7 --allow-extra-chr
./king -b ex7.bed --kinship --prefix Diabetes
