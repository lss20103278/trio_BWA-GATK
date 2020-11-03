#!/bin/sh

:<<Note_script_version
ver:2.6
for:WES/panel variants detection
update:2017/9/19
by:ykl
Note_script_version

:<<HELP
how to run?:
HELP

sample=$1
proc=$2

basedir=`cd $(dirname $0); pwd -P`
basedir="/DATA/sslyu/trio_BWA-GATK_v3.0/src/"

##env:
#source `pwd`/dbevn.sh
source ../dbevn.sh

cmd="clean mapping sorting precalling gvcf index_single gtgvcf hardfilt left-normalize depth qcsum qcvcf annotation"
##-->Workflow Control
#usage:select cmds from following
#full functions : "qc clean mapping sorting precalling calling gvcf gtgvcf vqsr hardfilt annotation depth abase qcsum qcvcf"
# for panel single
#cmd="clean mapping sorting precalling gvcf index_single gtgvcf hardfilt left-normalize depth qcsum qcvcf annotation"
# for trio
#cmd="clean mapping sorting precalling gvcf index_trio depth qcsum"
# for wxs single
#cmd="clean mapping sorting precalling gvcf index_single gtgvcf vqsr left-normalize depth qcsum qcvcf annotation"
:<<Note_pepline_main
rawdata -> clean -> mapping -> sorting -> precalling (bqsr)-> calling -> vqsr(WXS only) -> annotation

Thus,use:
cmd="clean mapping sorting precalling calling vqsr annotation"
Note_pepline_main

:<<Note_pepline_sub1
rawdata -> clean -> mapping -> sorting -> precalling (bqsr)-> gvcf -> gtgvcf -> vqsr -> annotation

Thus,use:
cmd="clean mapping sorting precalling gvcf vqsr annotation"
Note_pepline_sub1

:<<Note_pepline_sub2
rawdata -> clean -> mapping -> sorting -> precalling (bqsr)-> calling -> hardfilt(PANEL only) -> annotation
Note_pepline_sub2



##-->Name of raw files
##tips:make sure the rules of sample names
raw_data1=../raw/$sample\_*R1*.fastq.gz
raw_data2=../raw/$sample\_*R2*.fastq.gz


##----------------------------------------------------------------check functions
function runlog(){
echo -e "\033[31m##$1 \n\033[0m"
}
#usage:runlog "loginfos"

function check_file(){
if [ ! -e $1 ];then
echo -e "\033[31mError:$1 not exit,exiting...\033[0m"
exit
fi
}
#usage:check_file filename
##----------------------------------------------------------------check functions



##-->Check Path parameters:
check_file $ref_genome
check_file $panel
check_file $ref_1000G_indel
check_file $ref_1000M_indel
check_file $ref_snp
check_file $adapter
##
echo "start `date`" >> finished
sed -i "/^$/d" finished



##-->starting:
##cmd:qc
if [ "`echo $cmd|sed "s/ /\n/g"|awk '$1=="qc"{print "on"}'`" = "on" ];then
outdir="1_data"
runlog "running:qc_raw"
[ ! -e $outdir/qc_raw ] && mkdir -p $outdir/qc_raw
fastqc $raw_data1 -o $outdir/qc_raw > $outdir/qc_raw/qc_1.log 
fastqc $raw_data2 -o $outdir/qc_raw > $outdir/qc_raw/qc_2.log 
echo "qc `date`" >> finished
fi


##cmd:clean
if [ "`echo $cmd|sed "s/ /\n/g"|awk '$1=="clean"{print "on"}'`" = "on" ];then
sed -i "/^clean/d" finished
outdir="1_data"
runlog "running:clean & qc_clean"
[ ! -e $outdir/qc_clean ] && mkdir -p $outdir/qc_clean
#trim_galore -q 20 --phred33 --fastqc --stringency 1 --fastqc_args "--outdir $outdir/qc_clean" -e 0.1 --length 35 -o $outdir --paired $raw_data1 $raw_data2 > $outdir/run_log_trim
#mv $outdir/*$sample*R1*.fq $outdir/$sample\_R1.fastq
#mv $outdir/*$sample*R2*.fq $outdir/$sample\_R2.fastq

fastq-mcf -l 50 -o $outdir/$sample\_R1.fastq -o $outdir/$sample\_R2.fastq $adapter $raw_data1 $raw_data2
echo "clean `date`" >> finished
fi


##cmd:mapping
if [ "`echo $cmd|sed "s/ /\n/g"|awk '$1=="mapping"{print "on"}'`" = "on" ];then
sed -i "/^mapping/d" finished
outdir="2_mapping"
runlog "mapping"
[ ! -e $outdir ] && mkdir -p $outdir
#inRead1=$raw_data1
#inRead2=$raw_data2
inRead1=1_data/$sample\_R1.fastq
inRead2=1_data/$sample\_R2.fastq

check_file $inRead1
check_file $inRead2
bwa mem -t $proc -M $bwa_index -R "@RG\tID:$sample\tSM:$sample\tPL:Illumina" $inRead1 $inRead2 > $outdir/$sample.sam
echo "mapping `date`" >> finished
fi


##cmd:sorting
if [ "`echo $cmd|sed "s/ /\n/g"|awk '$1=="sorting"{print "on"}'`" = "on" ];then
sed -i "/^sorting/d" finished
outdir="2_mapping"
runlog "running:sorting:sam->bam"
check_file $outdir/$sample.sam
samtools view -b $outdir/$sample.sam > $outdir/$sample.bam

runlog "running:sorting:sort"
samtools sort $outdir/$sample.bam -o $outdir/$sample.sort.bam

runlog "running:sorting:index"
samtools index $outdir/$sample.sort.bam

runlog "running:sorting:mark-duplicates"
java -jar $soft_path/picard.jar MarkDuplicates INPUT=$outdir/$sample.sort.bam OUTPUT=/BP12_share/sslyu/bam/$sample.sort.mkdup.bam METRICS_FILE=$outdir/$sample.metrics

runlog "running:sorting:index-mkduplicated"
samtools index /BP12_share/sslyu/bam/$sample.sort.mkdup.bam

rm $outdir/$sample.sam
rm $outdir/$sample.bam
rm $outdir/$sample.sort.bam
rm $outdir/$sample.sort.bam.bai
echo "sorting `date`" >> finished
fi


##cmd:depth
if [ "`echo $cmd|sed "s/ /\n/g"|awk '$1=="depth"{print "on"}'`" = "on" ];then
sed -i "/^depth/d" finished
runlog "running:depth"
outdir="2_mapping"
[ ! -e $outdir ] && mkdir -p $outdir
check_file /BP12_share/sslyu/bam/$sample.sort.mkdup.bam
java -jar $soft_path/GATK/3.7/GenomeAnalysisTK.jar \
    -T DepthOfCoverage \
    -R $ref_genome \
    -L $panel \
    -o $outdir/$sample.depth \
    -I /BP12_share/sslyu/bam/$sample.sort.mkdup.bam \
    -geneList $ref_annotation/refseq.sorted.txt -ct 1 -ct 10 -ct 20 # It seems that GATK4 hasn't wrapped this function in
#samtools depth -b $panel /BP12_share/sslyu/bam/$sample.sort.mkdup.bam > $outdir/$sample.depth    
echo "depth `date`" >> finished
fi


##cmd:abase
if [ "`echo $cmd|sed "s/ /\n/g"|awk '$1=="abase"{print "on"}'`" = "on" ];then
sed -i "/^abase/d" finished
runlog "running:abase"
outdir="2_mapping"
[ ! -e $outdir ] && mkdir -p $outdir
check_file /BP12_share/sslyu/bam/$sample.sort.mkdup.bam
java -jar $soft_path/picard.jar CollectAlignmentSummaryMetrics R=$ref_genome \
 INPUT=/BP12_share/sslyu/bam/$sample.sort.mkdup.bam \
 O=$outdir/$sample.aln

echo "abase `date`" >> finished
fi


##cmd:qcsum
if [ "`echo $cmd|sed "s/ /\n/g"|awk '$1=="qcsum"{print "on"}'`" = "on" ];then
sed -i "/^qcsum/d" finished
runlog "running:qcsum"
outdir="2_mapping"
check_file $outdir/$sample.depth
echo -e "sample\t1Xcov\t20Xcov\tAvgDepth\tduplication" > $outdir/$sample.qcsum
printf $sample"\t" >> $outdir/$sample.qcsum
awk -F "\t" '{depth+=$4} $4>=1{sum1++} $4>=20{sum20++} END{printf (sum1+0.0)/(NR-1)*100"\t"(sum20+0.0)/(NR-1)*100"\t"(depth+0.0)/(NR-1)"\t"}' $outdir/$sample.depth >> $outdir/$sample.qcsum # gatk
#awk -F "\t" '{depth+=$3} $3>=1{sum1++} $3>=20{sum20++} END{printf (sum1+0.0)/(NR-1)*100"\t"(sum20+0.0)/(NR-1)*100"\t"(depth+0.0)/(NR-1)"\t"}' $outdir/$sample.depth >> $outdir/$sample.qcsum # samtools
[ -f $outdir/$sample.metrics ] && grep PERCENT_DUPLICATION -A 1 $outdir/$sample.metrics | tail -n 1 | awk '{print $(NF-1)}' >> $outdir/$sample.qcsum
echo "qcsum `date`" >> finished
fi


##cmd:precalling
if [ "`echo $cmd|sed "s/ /\n/g"|awk '$1=="precalling"{print "on"}'`" = "on" ];then
sed -i "/^precalling/d" finished
outdir="3_variants"
[ ! -e $outdir ] && mkdir -p $outdir
#runlog "running:precalling:realignment"
#java -jar $soft_path/GATK/3.7/GenomeAnalysisTK.jar  -T RealignerTargetCreator  -U ALLOW_N_CIGAR_READS  -nt $proc  -R $ref_genome -L $panel  -known $ref_1000G_indel  -known $ref_1000M_indel  -I /BP12_share/sslyu/bam/$sample.sort.mkdup.bam  -o $outdir/$sample.intervals # RealignerTargetCreator has been retired from GATK4.0
#
#check_file $outdir/$sample.intervals
#
#java -jar $soft_path/GATK/3.7/GenomeAnalysisTK.jar  -T IndelRealigner  -R $ref_genome  -known $ref_1000G_indel  -known $ref_1000M_indel  -I /BP12_share/sslyu/bam/$sample.sort.mkdup.bam  -targetIntervals $outdir/$sample.intervals  -o $outdir/$sample.realn.bam # IndelRealigner has been retired from GATK4.0

runlog "running:precalling:BQSR"
#check_file $outdir/$sample.realn.bam
#
#java -jar $soft_path/GATK/3.7/GenomeAnalysisTK.jar -T BaseRecalibrator -nct $proc -knownSites $ref_snp -knownSites $ref_1000G_indel -knownSites $ref_1000M_indel -R $ref_genome -L $panel -I $outdir/$sample.realn.bam -o $outdir/$sample.realn.recal # In GATK4.0, BaseRecalibrator is no longer able to do on the fly recalibration, so you have to run it on the recalibrated bam to do the QC cycle
#/DATA/sslyu/soft/gatk-4.1.0.0/gatk BaseRecalibratorSpark -I /BP12_share/sslyu/bam/$sample.sort.mkdup.bam -R $ref_genome_spark --known-sites $ref_snp --known-sites $ref_1000G_indel --known-sites $ref_1000M_indel -O $outdir/$sample.realn.recal -L $panel --spark-master local[4]
/DATA/sslyu/soft/gatk-4.1.0.0/gatk BaseRecalibrator -I /BP12_share/sslyu/bam/$sample.sort.mkdup.bam -R $ref_genome --known-sites $ref_snp --known-sites $ref_1000G_indel --known-sites $ref_1000M_indel -O $outdir/$sample.realn.recal -L $panel

check_file $outdir/$sample.realn.recal

runlog "running:precalling:BQSR:PrintReads"
#java -jar $soft_path/GATK/3.7/GenomeAnalysisTK.jar -T PrintReads -nct $proc -R $ref_genome -I $outdir/$sample.realn.bam -o $outdir/$sample.realn.recal.bam --BQSR $outdir/$sample.realn.recal # In GATK4.0 PrintReads --BQSR is replaced by ApplyBQSR
#/DATA/sslyu/soft/gatk-4.1.0.0/gatk ApplyBQSRSpark -I /BP12_share/sslyu/bam/$sample.sort.mkdup.bam -bqsr $outdir/$sample.realn.recal -O $outdir/$sample.realn.recal.bam -R $ref_genome_spark --spark-master local[4] # It seems that spark runs more slowly, don't know why
/DATA/sslyu/soft/gatk-4.1.0.0/gatk ApplyBQSR -I /BP12_share/sslyu/bam/$sample.sort.mkdup.bam -bqsr $outdir/$sample.realn.recal -O $outdir/$sample.realn.recal.bam -R $ref_genome

rm $outdir/$sample.realn.ba*

echo "precalling `date`" >> finished
fi


##cmd:calling
if [ "`echo $cmd|sed "s/ /\n/g"|awk '$1=="calling"{print "on"}'`" = "on" ];then
sed -i "/^calling/d" finished
outdir="3_variants"
runlog "running:calling"
check_file $outdir/$sample.realn.recal.bam
java -jar $soft_path/GATK/3.7/GenomeAnalysisTK.jar -T HaplotypeCaller -R $ref_genome -nct $proc -L $panel --dbsnp $ref_snp -stand_call_conf 10 -newQual -I $outdir/$sample.realn.recal.bam -o $outdir/$sample.raw.vcf
#java -Xmx32g -jar $soft_path/GenomeAnalysisTK.jar -T UnifiedGenotyper -nct $proc -R $ref_path/$ref_genome -glm BOTH -L $panel --dbsnp $ref_snp -stand_call_conf 30 -stand_emit_conf 10 -I $outdir/$sample.realn.recal.bam -o $outdir/$sample.raw.vcf

echo "calling `date`" >> finished
fi


##cmd:gvcf
if [ "`echo $cmd|sed "s/ /\n/g"|awk '$1=="gvcf"{print "on"}'`" = "on" ];then
sed -i "/^gvcf/d" finished
outdir="3_variants"
[ ! -e $outdir/gvcf ] && mkdir -p $outdir/gvcf
runlog "running:calling"
check_file $outdir/$sample.realn.recal.bam
/DATA/sslyu/soft/gatk-4.1.0.0/gatk HaplotypeCaller -R $ref_genome -I $outdir/$sample.realn.recal.bam -O $outdir/gvcf/$sample.g.vcf -ERC GVCF --dbsnp $ref_snp -L $panel -new-qual # It seems that in GATK4.0 there is no --variant_index_type and --variant_index_parameter options
rm $outdir/$sample.realn.recal.bam*
echo "gvcf `date`" >> finished
fi

##cmd:index_single
if [ "`echo $cmd|sed "s/ /\n/g"|awk '$1=="index_single"{print "on"}'`" = "on" ];then
sed -i "/^index_single/d" finished
outdir="3_variants"
runlog "running:index_single"
[ ! -d ../gvcf ] && mkdir -p ../gvcf
if [ ! -e ../gvcf/$sample.g.vcf ];then
cp $outdir/gvcf/$sample.g.vcf* ../gvcf/
fi
if [ ! -e ../gvcf/$sample.g.vcf.idx ];then
/DATA/sslyu/soft/gatk-4.1.0.0/gatk IndexFeatureFile -F ../gvcf/$sample.g.vcf
fi
echo "index_single `date`" >> finished
fi

##cmd:index_trio
if [ "`echo $cmd|sed "s/ /\n/g"|awk '$1=="index_trio"{print "on"}'`" = "on" ];then
sed -i "/^index_trio/d" finished
outdir="3_variants"
runlog "running:index_trio"
[ ! -d ../gvcf ] && mkdir -p ../gvcf
if [ ! -e ../gvcf/$sample.raw.g.vcf ];then
cp $outdir/gvcf/$sample.g.vcf ../gvcf/$sample.raw.g.vcf
fi
$basedir/gvcf.py ../gvcf/$sample
if [ ! -e ../gvcf/$sample.g.vcf.idx ];then
/DATA/sslyu/soft/gatk-4.1.0.0/gatk IndexFeatureFile -F ../gvcf/$sample.g.vcf #when using GATK4.0 and sentieon
fi
echo "index_trio `date`" >> finished
fi


##cmd:gtgvcf
if [ "`echo $cmd|sed "s/ /\n/g"|awk '$1=="gtgvcf"{print "on"}'`" = "on" ];then
sed -i "/^gtgvcf/d" finished
outdir="3_variants"
/DATA/sslyu/soft/gatk-4.1.0.0/gatk GenotypeGVCFs -R $ref_genome -L $panel --dbsnp $ref_snp -stand-call-conf 10 -O $outdir/$sample.raw.vcf -V $outdir/gvcf/$sample.g.vcf

echo "gtgvcf `date`" >> finished
fi

##cmd:hardfilt
if [ "`echo $cmd|sed "s/ /\n/g"|awk '$1=="hardfilt"{print "on"}'`" = "on" ];then
sed -i "/^hardfilt/d" finished

runlog "running:calling:select-snp"
outdir="3_variants"
/DATA/sslyu/soft/gatk-4.1.0.0/gatk IndexFeatureFile -F $outdir/$sample.raw.vcf
/DATA/sslyu/soft/gatk-4.1.0.0/gatk SelectVariants -R $ref_genome -select-type SNP --variant $outdir/$sample.raw.vcf -O $outdir/$sample.raw.snp.vcf

runlog "running:calling:filt-snp"
/DATA/sslyu/soft/gatk-4.1.0.0/gatk VariantFiltration -R $ref_genome --variant $outdir/$sample.raw.snp.vcf -O $outdir/$sample.raw.snp.fil.vcf --filter-name "REJECT" --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0"

runlog "running:calling:select-indel"
/DATA/sslyu/soft/gatk-4.1.0.0/gatk SelectVariants -R $ref_genome -select-type INDEL --variant $outdir/$sample.raw.vcf -O $outdir/$sample.raw.indel.vcf

runlog "running:calling:filt-indel"
/DATA/sslyu/soft/gatk-4.1.0.0/gatk VariantFiltration -R $ref_genome --variant $outdir/$sample.raw.indel.vcf -O $outdir/$sample.raw.indel.fil.vcf --filter-name "REJECT" --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0"

runlog "running:calling:select-passed-snp/indel"
awk '/^#/ || $7 == "PASS"' $outdir/$sample.raw.snp.fil.vcf >$outdir/$sample.pass.snp.vcf
awk '/^#/ || $7 == "PASS"' $outdir/$sample.raw.indel.fil.vcf >$outdir/$sample.pass.indel.vcf

runlog "running:calling:catvariants"
/DATA/sslyu/soft/gatk-4.1.0.0/gatk SortVcf -R $ref_genome -I $outdir/$sample.pass.snp.vcf -I $outdir/$sample.pass.indel.vcf -O $outdir/$sample.mul.vcf

echo "hardfilt `date`" >> finished
fi

##cmd:vqsr
if [ "`echo $cmd|sed "s/ /\n/g"|awk '$1=="vqsr"{print "on"}'`" = "on" ];then
sed -i "/^vqsr/d" finished
runlog "running:vqsr"
outdir="3_variants/vqsr"
[ ! -e $outdir ] && mkdir -p $outdir
/DATA/sslyu/soft/gatk-4.0.8.1/gatk VariantRecalibrator -R $ref_genome -L $panel \
 -V 3_variants/$sample.raw.vcf \
 --resource hapmap,known=false,training=true,truth=true,prior=15.0:$ref_vqsr_hapmap \
 --resource omni,known=false,training=true,truth=false,prior=12.0:$ref_vqsr_omni \
 --resource 1000G,known=false,training=true,truth=false,prior=10.0:$ref_vqsr_1000G \
 --resource dbsnp,known=true,training=false,truth=false,prior=2.0:$ref_snp \
 -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR --max-gaussians 4 \
 -mode SNP \
 -O $outdir/snp.recal --tranches-file $outdir/snp.tranches --rscript-file $outdir/snp.plots.R # In GATK4, --recal-file in VariantRcalibrator is replace by -O

/DATA/sslyu/soft/gatk-4.1.0.0/gatk ApplyVQSR -R $ref_genome -L $panel \
 -V 3_variants/$sample.raw.vcf \
 -ts-filter-level 99.0 --recal-file $outdir/snp.recal --tranches-file $outdir/snp.tranches \
 -mode SNP \
 -O $outdir/$sample.raw.vqsr.snp.vcf 

/DATA/sslyu/soft/gatk-4.0.8.1/gatk VariantRecalibrator -R $ref_genome -L $panel \
 -V 3_variants/$sample.raw.vcf \
 --resource mills,known=false,training=true,truth=true,prior=12.0:$ref_1000M_indel \
 --resource dbsnp,know=true,training=false,truth=false,prior=2.0:$ref_snp \
 -an QD -an DP -an FS -an SOR -an MQRankSum -an ReadPosRankSum --max-gaussians 4 \
 -mode INDEL \
 -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
 -O $outdir/indel.recal --rscript-file $outdir/indels.plots.R  --tranches-file $outdir/indel.tranches

/DATA/sslyu/soft/gatk-4.1.0.0/gatk ApplyVQSR -R $ref_genome -L $panel \
 -V 3_variants/$sample.raw.vcf \
 -ts-filter-level 99.0 --recal-file $outdir/indel.recal --tranches-file $outdir/indel.tranches \
 -mode INDEL \
 -O $outdir/$sample.raw.vqsr.indel.vcf
 
#filt pass:
awk -F "\t" '$1~/^#/ ||$7=="PASS" {print}' $outdir/$sample.raw.vqsr.snp.vcf > $outdir/$sample.raw.vqsr.snp.pass.vcf
awk -F "\t" '$1~/^#/ ||$7=="PASS" {print}' $outdir/$sample.raw.vqsr.indel.vcf > $outdir/$sample.raw.vqsr.indel.pass.vcf

#combine:
/DATA/sslyu/soft/gatk-4.1.0.0/gatk SortVcf -R $ref_genome -I $outdir/$sample.raw.vqsr.snp.pass.vcf -I $outdir/$sample.raw.vqsr.indel.pass.vcf -O 3_variants/$sample.mul.vcf

echo "vqsr `date`" >> finished
fi

##cmd:left-normalize
if [ "`echo $cmd|sed "s/ /\n/g"|awk '$1=="left-normalize"{print "on"}'`" = "on" ];then
runlog "running:left-normalize"
bcftools norm -m - 3_variants/$sample.mul.vcf | bcftools norm -f $ref_genome > 3_variants/$sample.vcf
echo "left-normalize `date`" >> finished
fi

##cmd:annotation
if [ "`echo $cmd|sed "s/ /\n/g"|awk '$1=="annotation"{print "on"}'`" = "on" ];then
sed -i "/^annotation/d" finished
runlog "running:annotation"
outdir="3_variants"
check_file $outdir/$sample.vcf
#perl /DATA/ypliu/opt/mskcc-vcf2maf-1b16b35/vcf2maf.pl --input-vcf $outdir/$sample.vcf --output-maf $outdir/$sample.maf --vep-path /DATA/ypliu/opt/mskcc-vcf2maf-1b16b35/vep --vep-data /DATA/ypliu/opt/mskcc-vcf2maf-1b16b35/.vep --ref-fasta /DATA/sslyu/refGene/hg19.fa --filter-vcf /DATA/ypliu/opt/mskcc-vcf2maf-1b16b35/.vep/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz
perl /DATA/ypliu/opt/mskcc-vcf2maf-1b16b35/vcf2maf.pl --input-vcf $outdir/$sample.vcf --output-maf $outdir/$sample.maf --vep-path /DATA/sslyu/soft/mskcc-vcf2maf-1b16b35/vep/ --vep-data /DATA/ypliu/opt/mskcc-vcf2maf-1b16b35/.vep --ref-fasta /DATA/sslyu/refGene/hg19.fa --filter-vcf /DATA/ypliu/opt/mskcc-vcf2maf-1b16b35/.vep/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz
sed -i '1d' $outdir/$sample.maf
/DATA/sslyu/soft/annovar/convert2annovar.pl -format vcf4 $outdir/$sample.vcf -outfile $outdir/$sample.avinput
/DATA/sslyu/soft/annovar/convert2annovar.pl -format vcf4 $outdir/$sample.vcf -outfile $outdir/$sample.link --includeinfo
sed 's/^chr//' $outdir/$sample.vcf|gzip -c ->$outdir/$sample.vcf.gz
#score.sh $outdir/$sample.vcf.gz $outdir/$sample.CADD.gz # Can't locate CGI.pm in @INC...
export VEPpath="/DATA/sslyu/soft/variant_effect_predictor/variant_effect_predictor.pl"
/DATA/sslyu/soft/CADD_v1.3/bin/score.sh $outdir/$sample.vcf.gz $outdir/$sample.CADD.gz
bgzip -df $outdir/$sample.CADD.gz
sed -i '1d' $outdir/$sample.CADD
/DATA/sslyu/soft/annovar/table_annovar.pl $outdir/$sample.avinput $ref_annotation -buildver hg19 -out $outdir/$sample.ann -protocol refGene,avsnp147,1000g2015aug_all,1000g2015aug_eas,gnomad_genome_eas,gnomad_exomes_hom,gnomad_exomes_hom_all,gnomad_exomes_hemi,gnomad_exomes_r2.1,exac03_eas,esp6500_all,dbscsnv11,dbnsfp33a,revel,intervar_20180118,clinvar_20190211,ensembl,hgmd_2018_spring,bed,bed,bed,bed -operation g,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,r,r,r,r -bedfile hg19_Dead_Zone.txt,hg19_Problem_High.txt,hg19_Problem_Low.txt,hg19_HI_Predictions_Version3.bed -arg ',,,,,,,,,,,,,,,,,,-colsWanted 4,-colsWanted 4,-colsWanted 4,-colsWanted 4' -remove -nastring "NA" -otherinfo # annovar 20180416 intervar_201702->20180118 clinvar_20180603->20181225->20190211 gnomad:r2.1
sed -i "1s/bed\tbed2\tbed3\tbed4/Dead_Zone\tProblem_High\tProblem_Low\tHI_Predictions/" $outdir/$sample.ann\.hg19_multianno.txt
sed -i "s/Otherinfo/Heterozygosity\tQual\tdepth/" $outdir/$sample.ann\.hg19_multianno.txt #GenotypingQuality-->Qual
echo "annotation `date`" >> finished
fi

##cmd:qcvcf
if [ "`echo $cmd|sed "s/ /\n/g"|awk '$1=="qcvcf"{print "on"}'`" = "on" ];then
sed -i "/^qcvcf/d" finished
runlog "running:qcvcf"
check_file 3_variants/$sample.vcf
echo -e "sample\tindel_genomic\tsnp_genomic\tts/tv_genomic\tsnp_exonic\tts/tv_exonic" > 3_variants/$sample.qcvcf
printf "$sample\t" >> 3_variants/$sample.qcvcf
[ ! -d tmprun ] && mkdir tmprun
/DATA/sslyu/soft/annovar/convert2annovar.pl -format vcf4 3_variants/$sample.vcf -outfile tmprun/$sample.avinput
/DATA/sslyu/soft/annovar/table_annovar.pl tmprun/$sample.avinput $ref_annotation -buildver hg19 -out tmprun/$sample -protocol refGene -operation g -remove

awk -F "\t" 'BEGIN{indel=0}  $4=="-"||$5=="-"{indel++} $4!="-"&&$5!="-"{print $4$5} END{print indel"\t"NR-1-indel"\t"}' tmprun/$sample.hg19_multianno.txt | awk 'BEGIN{ts=0} $1=="AG"||$1=="GA"||$1=="CT"||$1=="TC"{ts++} END{printf $0;tv=NR-2-ts;printf (ts+0.0)/tv"\t"}' >> 3_variants/$sample.qcvcf

awk -F "\t" '$6=="exonic"||$6=="exonic;splicing"{print $4$5}' tmprun/$sample.hg19_multianno.txt | awk '$1!~/-/{print}'| awk 'BEGIN{ts=0} $1=="AG"||$1=="GA"||$1=="CT"||$1=="TC"{ts++} END{tv=NR-1-ts;printf ts+tv"\t"(ts+0.0)/tv}' >> 3_variants/$sample.qcvcf

rm -rf tmprun
echo "qcvcf `date`" >> finished
fi
