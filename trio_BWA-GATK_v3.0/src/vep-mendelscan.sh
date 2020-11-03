#!/bin/sh

sn=$1
ped=$2

##env:
source `pwd`/../../dbevn.sh


outdir="segtrio"
[ ! -d $outdir ] && mkdir $outdir


#seperate multiple alles
#bcftools norm -m - $sn.gt.vcf | bcftools norm -f $ref_genome >$outdir/$sn.vcf
bcftools norm -m - $sn.gt.vcf | bcftools norm -f $ref_genome |sed -n '/^chrM/!p' >$outdir/$sn.vcf # 去掉线粒体变异，CADD不再支持对线粒体变异打分

#annotation:vep
#variant_effect_predictor.pl -i $outdir/$sn.vcf  -o $outdir/$sn.vep --cache --offline --dir /SSD750/PB3/db3/Homo/vep/ --sift b --polyphen b --symbol --canonical --force # Can't locate CGI.pm in @INC...
perl /DATA/sslyu/soft/variant_effect_predictor/variant_effect_predictor.pl -i $outdir/$sn.vcf  -o $outdir/$sn.vep --cache --offline --dir /SSD750/PB3/db3/Homo/vep/ --sift b --polyphen b --symbol --canonical --force

#annotation:mendelscan
java -jar $soft_path/MendelScan/MendelScan.v1.2.1.fix.jar score $outdir/$sn.vcf --vep-file $outdir/$sn.vep --ped-file $ped --output-file $outdir/$sn.mendel.tsv --output-vcf $outdir/$sn.mendel.vcf --inheritance recessive

#perl /DATA/ypliu/opt/mskcc-vcf2maf-1b16b35/vcf2maf.pl --input-vcf $outdir/$sn.vcf --output-maf $outdir/$sn.maf --vep-path /DATA/ypliu/opt/mskcc-vcf2maf-1b16b35/vep --vep-data /DATA/ypliu/opt/mskcc-vcf2maf-1b16b35/.vep --ref-fasta /DATA/sslyu/refGene/hg19.fa --filter-vcf /DATA/ypliu/opt/mskcc-vcf2maf-1b16b35/.vep/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz # Can't locate CGI.pm in @INC...
perl /DATA/ypliu/opt/mskcc-vcf2maf-1b16b35/vcf2maf.pl --input-vcf $outdir/$sn.vcf --output-maf $outdir/$sn.maf --vep-path /DATA/sslyu/soft/mskcc-vcf2maf-1b16b35/vep/ --vep-data /DATA/ypliu/opt/mskcc-vcf2maf-1b16b35/.vep --ref-fasta /DATA/sslyu/refGene/hg19.fa --filter-vcf /DATA/ypliu/opt/mskcc-vcf2maf-1b16b35/.vep/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz 
sed -i '1d' $outdir/$sn.maf

#/DATA/sslyu/soft/annovar/convert2annovar.pl -format vcf4 $outdir/$sn.mendel.vcf -outfile $outdir/$sn.link --includeinfo # annovar 20180416
/DATA/sslyu/soft/annovar/convert2annovar.pl -format vcf4 -allsample -withfreq $outdir/$sn.mendel.vcf -outfile $outdir/$sn.link --includeinfo # annovar 20180416 no longer does line-to-line conversion for multi-sample VCF files, include all variants in output, use '-format vcf4old' or use '-format vcf4 -allsample -withfreq'
#sed 's/^chr//' $outdir/$sn.vcf|gzip -c ->$outdir/$sn.vcf.gz
gzip -c $outdir/$sn.mendel.vcf>$outdir/$sn.mendel.vcf.gz
#score.sh $outdir/$sn.mendel.vcf.gz $outdir/$sn.CADD.gz # Can't locate CGI.pm in @INC...
export VEPpath="/DATA/sslyu/soft/variant_effect_predictor/variant_effect_predictor.pl" 
/DATA/sslyu/soft/CADD_v1.3/bin/score.sh $outdir/$sn.mendel.vcf.gz $outdir/$sn.CADD.gz  
bgzip -df $outdir/$sn.CADD.gz
sed -i '1d' $outdir/$sn.CADD

#annotation:annovar
#./mendel_to_annovar.sh $sn
/DATA/sslyu/trio_BWA-GATK_v3.0/src/mendel_to_annovar.py $sn

/DATA/sslyu/soft/annovar/table_annovar.pl $outdir/$sn.avinput $ref_annotation -buildver hg19 -out $outdir/$sn.ann -protocol refGene,avsnp147,1000g2015aug_all,1000g2015aug_eas,gnomad_genome_eas,gnomad_exomes_hom,gnomad_exomes_hom_all,gnomad_exomes_hemi,exac03_eas,esp6500_all,dbscsnv11,dbnsfp33a,revel,intervar_20180118,clinvar_20190211,ensembl,hgmd_2018_spring,bed,bed,bed,bed -operation g,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,r,r,r,r -bedfile hg19_Dead_Zone.txt,hg19_Problem_High.txt,hg19_Problem_Low.txt,hg19_HI_Predictions_Version3.bed -arg ',,,,,,,,,,,,,,,,,-colsWanted 4,-colsWanted 4,-colsWanted 4,-colsWanted 4' -remove -nastring "NA" -otherinfo # annovar:20180416 intervar:20170202->20180118 clinvar:20180603->20181225->20190211 
#/DATA/sslyu/soft/annovar/table_annovar.pl $outdir/$sn.avinput $ref_annotation -buildver hg19 -out $outdir/$sn.ann -protocol refGene,avsnp147,1000g2015aug_all,1000g2015aug_eas,gnomad_genome_eas,gnomad_exomes_hom,gnomad_exomes_hom_all,gnomad_exomes_hemi,gnomad_exomes_r2.1,exac03_eas,esp6500_all,dbscsnv11,dbnsfp33a,revel,intervar_20180118,clinvar_20190211,ensembl,hgmd_2018_spring,bed,bed,bed,bed -operation g,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,r,r,r,r -bedfile hg19_Dead_Zone.txt,hg19_Problem_High.txt,hg19_Problem_Low.txt,hg19_HI_Predictions_Version3.bed -arg ',,,,,,,,,,,,,,,,,,-colsWanted 4,-colsWanted 4,-colsWanted 4,-colsWanted 4' -remove -nastring "NA" -otherinfo # annovar:20180416 intervar:20170202->20180118 clinvar:20180603->20181225->20190211 gnomad:r2.1
sed -i "1s/bed\tbed2\tbed3\tbed4/Dead_Zone\tProblem_High\tProblem_Low\tHI_Predictions/" $outdir/$sn.ann\.hg19_multianno.txt

otherinfo=`cat $outdir/$sn.avinput.info`
sed -i "s/Otherinfo/`echo $otherinfo|sed "s/ /\t/g"`/" $outdir/$sn.ann.hg19_multianno.txt


