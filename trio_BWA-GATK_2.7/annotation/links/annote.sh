#!/bin/sh

soft_path="/SSD750/PB1/soft1"
DB="/SSD750/PB1/db1/Homo"
ref_genome="$DB/refseq/hg19.fa"
bwa_index="$DB/bwa_index/hg19.fa"
ref_1000G_indel="$DB/GATK/1000G_phase1.indels.hg19.sites.vcf"
ref_1000M_indel="$DB/GATK/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf"
ref_annotation="$DB/annovar"
ref_snp="$DB/GATK/dbsnp_138.hg19.vcf"
panel="$DB/design/SCH.capture.hg19.bed"


for i in `cat list`
do
table_annovar.pl $i.avinput $ref_annotation -buildver hg19 -out $i.ann -protocol refGene,snp138,1000g2014oct_all,esp6500_all,exac03,ljb26_all,clinvar,ensembl,hgmd -operation g,f,f,f,f,f,f,f,f -remove -nastring "NA" 
done
