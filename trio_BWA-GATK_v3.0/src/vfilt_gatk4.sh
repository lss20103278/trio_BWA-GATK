#!/bin/sh


dir=$1
sn=$2
mode=$3


##env:
source `pwd`/../../dbevn.sh


cd $dir

#mode:vqsr
if [ "$mode" = "vqsr" ];then
#snp
/DATA/sslyu/soft/gatk-4.0.8.1/gatk VariantRecalibrator -R $ref_genome -L $panel \
 -V $sn.raw.vcf \
 --resource hapmap,known=false,training=true,truth=true,prior=15.0:$ref_vqsr_hapmap \
 --resource omni,known=false,training=true,truth=false,prior=12.0:$ref_vqsr_omni \
 --resource 1000G,known=false,training=true,truth=false,prior=10.0:$ref_vqsr_1000G \
 --resource dbsnp,known=true,training=false,truth=false,prior=2.0:$ref_snp \
 -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR --max-gaussians 4 \
 -mode SNP \
 -O snp.recal --tranches-file snp.tranches --rscript-file snp.plots.R # gatk-4.1.0.0 doesn't work and don't know why

/DATA/sslyu/soft/gatk-4.1.0.0/gatk ApplyVQSR -R $ref_genome -L $panel \
 -V $sn.raw.vcf \
 -ts-filter-level 99.0 --recal-file snp.recal --tranches-file snp.tranches \
 -mode SNP \
 -O $sn.raw.vqsr.snp.vcf

#indel
/DATA/sslyu/soft/gatk-4.0.8.1/gatk VariantRecalibrator -R $ref_genome -L $panel \
 -V $sn.raw.vcf \
 --resource mills,known=false,training=true,truth=true,prior=12.0:$ref_1000M_indel \
 --resource dbsnp,know=true,training=false,truth=false,prior=2.0:$ref_snp \
 -an QD -an DP -an FS -an SOR -an MQRankSum -an ReadPosRankSum --max-gaussians 4 \
 -mode INDEL \
 -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
 -O indel.recal --rscript-file indels.plots.R  --tranches-file indel.tranches 

/DATA/sslyu/soft/gatk-4.1.0.0/gatk ApplyVQSR -R $ref_genome -L $panel \
 -V $sn.raw.vcf \
 -ts-filter-level 99.0 --recal-file indel.recal --tranches-file indel.tranches \
 -mode INDEL \
 -O $sn.raw.vqsr.indel.vcf

#pass:
awk -F "\t" '$1~/^#/ ||$7=="PASS" {print}' $sn.raw.vqsr.snp.vcf > $sn.raw.snp.pass.vcf
awk -F "\t" '$1~/^#/ ||$7=="PASS" {print}' $sn.raw.vqsr.indel.vcf > $sn.raw.indel.pass.vcf

fi


##mode:hard
if [ "$mode" = "hard" ];then
/DATA/sslyu/soft/gatk-4.1.0.0/gatk SelectVariants -R $ref_genome -select-type SNP --variant $sn.raw.vcf -O $sn.raw.snp.vcf

/DATA/sslyu/soft/gatk-4.1.0.0/gatk VariantFiltration -R $ref_genome --variant $sn.raw.snp.vcf -O $sn.raw.snp.fil.vcf --filter-name "REJECT" --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0"

/DATA/sslyu/soft/gatk-4.1.0.0/gatk SelectVariants -R $ref_genome -select-type INDEL --variant $sn.raw.vcf -O $sn.raw.indel.vcf

/DATA/sslyu/soft/gatk-4.1.0.0/gatk VariantFiltration -R $ref_genome --variant $sn.raw.indel.vcf -O $sn.raw.indel.fil.vcf --filter-name "REJECT" --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0"

awk '/^#/ || $7 == "PASS"' $sn.raw.snp.fil.vcf >$sn.raw.snp.pass.vcf
awk '/^#/ || $7 == "PASS"' $sn.raw.indel.fil.vcf >$sn.raw.indel.pass.vcf

fi


cd ../
