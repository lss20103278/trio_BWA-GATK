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

##env:
#source `pwd`/dbevn.sh
source ../dbevn.sh


##-->Workflow Control
#usage:select cmds from following
#full functions : "qc clean mapping sorting precalling calling gvcf gtgvcf vqsr hardfilt left-normalize annotation depth abase qcsum qcvcf"
# for panel single
#cmd="clean mapping precalling gvcf gtgvcf hardfilt left-normalize depth qcsum qcvcf annotation"
# for trio
#cmd="clean mapping precalling gvcf depth qcsum"
# for wxs single
#cmd="clean mapping precalling gvcf gtgvcf vqsr left-normalize depth qcsum qcvcf annotation"
:<<Note_pepline_main
rawdata -> clean -> mapping -> sorting -> precalling (bqsr)-> calling -> vqsr(WXS only) -> left-normalize -> annotation

Thus,use:
cmd="clean mapping sorting precalling calling vqsr left-normalize annotation"
Note_pepline_main

:<<Note_pepline_sub1
rawdata -> clean -> mapping -> sorting -> precalling (bqsr)-> gvcf -> gtgvcf -> vqsr -> left-normalize -> annotation

Thus,use:
cmd="clean mapping sorting precalling gvcf vqsr left-normalize annotation"
Note_pepline_sub1

:<<Note_pepline_sub2
rawdata -> clean -> mapping -> sorting -> precalling (bqsr)-> calling -> hardfilt(PANEL only) -> left-normalize -> annotation
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
echo -e "\033[31mError:$1 not exist,exiting...\033[0m"
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
echo "" >> finished
sed -i "/^$/d" finished
echo "start `date`" >> finished

##-->starting:
##cmd:qc
if [ "`echo $cmd|sed "s/ /\n/g"|awk '$1=="qc"{print "on"}'`" = "on" ];then
outdir="1_data"
runlog "running:qc_raw"
[ ! -e $outdir/qc_raw ] && mkdir -p $outdir/qc_raw
fastqc $raw_data1 -o $outdir/qc_raw > $outdir/qc_raw/qc_1.log 
fastqc $raw_data2 -o $outdir/qc_raw > $outdir/qc_raw/qc_2.log 
fi

##cmd:clean
if [ "`echo $cmd|sed "s/ /\n/g"|awk '$1=="clean"{print "on"}'`" = "on" ];then
sed -i "/^clean/d" finished
outdir="1_data"
runlog "running:clean & qc_clean"
[ ! -e $outdir/qc_clean ] && mkdir -p $outdir/qc_clean

fastq-mcf -l 50 -o $outdir/$sample\_R1.fastq.gz -o $outdir/$sample\_R2.fastq.gz $adapter $raw_data1 $raw_data2

echo "clean `date`" >> finished
fi

##cmd:mapping and sorting
if [ "`echo $cmd|sed "s/ /\n/g"|awk '$1=="mapping"{print "on"}'`" = "on" ];then
sed -i "/^mapping/d" finished
outdir="2_mapping"
runlog "running:sorting:mapping,sam->bam"
[ ! -e $outdir ] && mkdir -p $outdir

inRead1=1_data/$sample\_R1.fastq.gz
inRead2=1_data/$sample\_R2.fastq.gz

check_file $inRead1
check_file $inRead2

(bwa mem -t $proc -M $bwa_index -R "@RG\tID:$sample\tSM:$sample\tPL:Illumina" $inRead1 $inRead2 || echo -n 'error') | sentieon util sort -o $outdir/$sample.sort.bam -t $proc --sam2bam -i - # 由bwa生成的sam文件是按字典式排序法进行的排序（chr10,chr11..chr19,chr1,chr20..chr22,chr2,chr3..chrM,chrX,chrY),但那会死gatk在进行call snp的时候是按照染色体组型进行的（chrM,chr1,chr2..chr22,chrX,chrY),因此要对原始sam文件进行reorder; -R参数对文件进行加头处理，因为gatk2.0以上版本不再支持无头文件的变异检测

sentieon driver -t $proc -r $ref_genome -i $outdir/$sample.sort.bam \
	--algo GCBias --summary $outdir/gc_summary.txt $outdir/gc_metric.txt \
	--algo MeanQualityByCycle $outdir/mq_metric.txt \
	--algo QualDistribution $outdir/qd_metric.txt \
	--algo InsertSizeMetricAlgo $outdir/is_metric.txt \
	--algo AlignmentStat $outdir/aln_metric.txt
	
sentieon plot metrics -o $outdir/metrics.pdf gc=$outdir/gc_metric.txt mq=$outdir/mq_metric.txt qd=$outdir/qd_metric.txt isize=$outdir/is_metric.txt

echo "mapping `date`" >> finished

runlog "running:sorting:rm-duplicates"  #PCR重复：如果两条reads具有相同的长度而且比对到基因组的同一位置，那么就认为这样的reads是由PCR扩增而来，就会被标记；在制备文库的过程中，由于PCR扩增过程中会存在一些偏差，也就是说有的序列会被过量扩增。这样，在比对的时候，这些过量扩增出来的完全相同的序列就会比对到基因组的相同位置。而这些过量扩增的reads并不是基因组自身固有序列，不能作为变异检测的证据，因此，要尽量去除这些由PCR扩增所形成的duplicates，

sentieon driver -t $proc -i $outdir/$sample.sort.bam --algo LocusCollector --fun score_info $outdir/score.txt

sentieon driver -t $proc -i $outdir/$sample.sort.bam --algo Dedup --rmdup --score_info $outdir/score.txt --metrics $outdir/$sample.metrics $outdir/$sample.sort.mkdup.bam
echo "sorting `date`" >> finished
fi

##cmd:depth
if [ "`echo $cmd|sed "s/ /\n/g"|awk '$1=="depth"{print "on"}'`" = "on" ];then
sed -i "/^depth/d" finished
runlog "running:depth"
outdir="2_mapping"
check_file $outdir/$sample.sort.mkdup.bam
sentieon driver -r $ref_genome -i $outdir/$sample.sort.mkdup.bam --interval $panel \
	 --algo CoverageMetrics --gene_list $ref_annotation/refseq.sorted.txt \
	 --cov_thresh 1 --cov_thresh 10 --cov_thresh 20 $outdir/$sample.depth

echo "depth `date`" >> finished
fi

##cmd:qcsum
if [ "`echo $cmd|sed "s/ /\n/g"|awk '$1=="qcsum"{print "on"}'`" = "on" ];then
sed -i "/^qcsum/d" finished
runlog "running:qcsum"
outdir="2_mapping"
check_file $outdir/$sample.depth
echo -e "sample\t1Xcov\t20Xcov\tAvgDepth\tduplication" > $outdir/$sample.qcsum
printf $sample"\t" >> $outdir/$sample.qcsum
awk -F "\t" '{depth+=$4} $4>=1{sum1++} $4>=20{sum20++} END{printf (sum1+0.0)/(NR-1)*100"\t"(sum20+0.0)/(NR-1)*100"\t"(depth+0.0)/(NR-1)"\t"}' $outdir/$sample.depth >> $outdir/$sample.qcsum
[ -f $outdir/$sample.metrics ] && grep PERCENT_DUPLICATION -A 1 $outdir/$sample.metrics | tail -n 1 | awk '{print $(NF-1)}' >> $outdir/$sample.qcsum  # metrics文件中的PERCENT_DUPLICATION是指PCR重复吗？
#rm $outdir/$sample.depth
echo "qcsum `date`" >> finished
fi

##cmd:precalling
if [ "`echo $cmd|sed "s/ /\n/g"|awk '$1=="precalling"{print "on"}'`" = "on" ];then
sed -i "/^precalling/d" finished
outdir="3_variants"
[ ! -e $outdir ] && mkdir -p $outdir
runlog "running:precalling:realignment"  #GATK4.0为了追求效率，不如3.8的准确性高；GATK速度慢并且准的原因是因为有一步local assembly，因此只对SNP感兴趣的话，我个人认为可以不用GATK

sentieon driver -t $proc -r $ref_genome -i 2_mapping/$sample.sort.mkdup.bam --algo Realigner -k $ref_1000G_indel -k $ref_1000M_indel --interval_list $panel $outdir/$sample.realn.bam  # realignment将比对到indel附近的reads进行局部重新比对，将比对的错误率降到最低。一般来说，绝大部分需要进行重新比对的基因组区域，都是因为插入/缺失的存在，因为在indel附近的比对会出现大量的碱基错配，这些碱基的错配很容易被误认为SNP。还有，在比对过程中，比对算法对于每一条read的处理都是独立的，不可能同时把多条reads与参考基因组比对来拍错。因此，即使有一些 reads能够正确的比对到indel，但那些恰恰比对到indel开始或者结束为止的read也会有很高的比对错误率，这都是需要重新比对的。

runlog "running:precalling:BQSR" #对bam文件里reads的碱基质量值进行重新校正，使最后输出的bam文件中reads中碱基的质量值能够更加接近真实值？在reads碱基质量值被校正之前，我们要保留质量值在Q25以上的碱基，但是实际上质量值在Q25的这些碱基的错误率在1%，也就是说质量值只有Q20，这样就会对后续的变异检测的可信度造成影响。还有，在边合成边测序的测序过程中，在reads末端碱基的错误率往往要比起始部位更高。另外，AC的质量值往往要低于TG。BQSR就是要对这些质量值进行校正。
#根据一些known sites，生成一个校正质量值所需要的数据文件，利用这个文件来生成校正后的数据文件，最后生成碱基质量值校正前后的比较图，利用工具将经过质量值校正的数据输出到新的bam文件中，用于后续的变异检测
check_file $outdir/$sample.realn.bam

sentieon driver -t $proc -r $ref_genome -i $outdir/$sample.realn.bam --interval $panel --algo QualCal -k $ref_snp -k $ref_1000G_indel -k $ref_1000M_indel $outdir/recal_data.table  # QualCal是干什么用的?

sentieon driver -t $proc -r $ref_genome -i $outdir/$sample.realn.bam --interval $panel -q $outdir/recal_data.table --algo QualCal -k $ref_snp -k $ref_1000G_indel -k $ref_1000M_indel $outdir/recal_data.table.post

sentieon driver -t $proc --algo QualCal --plot --before $outdir/recal_data.table --after $outdir/recal_data.table.post recal.csv

sentieon plot bqsr -o recal_plots.pdf recal.csv

echo "precalling `date`" >> finished
fi

##cmd:calling
if [ "`echo $cmd|sed "s/ /\n/g"|awk '$1=="calling"{print "on"}'`" = "on" ];then
sed -i "/^calling/d" finished
outdir="3_variants"
runlog "running:calling"
check_file $outdir/$sample.realn.bam

sentieon driver -t $proc -r $ref_genome --interval $panel -i $outdir/$sample.realn.bam -q $outdir/recal_data.table --algo Haplotyper -d $ref_snp --emit_conf=10 --call_conf=30 $outdir/$sample.raw.vcf  # germline call snp 用Haplotyper 

echo "calling `date`" >> finished
fi

##cmd:gvcf
if [ "`echo $cmd|sed "s/ /\n/g"|awk '$1=="gvcf"{print "on"}'`" = "on" ];then
sed -i "/^gvcf/d" finished
outdir="3_variants"
[ ! -e $outdir/gvcf ] && mkdir -p $outdir/gvcf
runlog "running:calling"
check_file $outdir/$sample.realn.bam

sentieon driver -t $proc -r $ref_genome --interval $panel -i $outdir/$sample.realn.bam -q $outdir/recal_data.table --algo Haplotyper \
	-d $ref_snp --emit_mode gvcf $outdir/gvcf/$sample.g.vcf  # 用Haplotyper call snp，输出gvcf；使用贝叶斯最大似然模型，估计基因型和基因频率，给出后验概率
# --interval 指定对基因组的某一区域进行变异检测； gvcf文件与vcf文件都是vcf文件，不同的是gvcf文件同时会记录为蜕变的微店的覆盖情况，在多个样本的vcf文件进行合并的时候，需要区分./.和0/0的情况，./.是未检出的基因型，而0/0时未突变的基因型，如果仅使用普通的v 吃饭文件进行合并，那么就无法区分出这两种情况，进而对合并结果产生偏差。

#sentieon driver -t $proc -r $ref_genome --interval $panel -i $outdir/$sample.realn.bam -q $outdir/recal_data.table --algo Haplotyper \
#	-d $ref_snp --emit_conf=10 --call_conf=10 --emit_mode gvcf $outdir/gvcf/$sample.g.vcf

echo "gvcf `date`" >> finished
fi

##cmd:gtgvcf
if [ "`echo $cmd|sed "s/ /\n/g"|awk '$1=="gtgvcf"{print "on"}'`" = "on" ];then
sed -i "/^gtgvcf/d" finished
outdir="3_variants"

sentieon driver -t $proc -r $ref_genome --interval $panel --algo GVCFtyper --call_conf=10 -d $ref_snp -v $outdir/gvcf/$sample.g.vcf $outdir/$sample.raw.vcf # GVCFtyper将多个gvcf文件合并起来，但是这里只有一个gvcf文件，所以这一步？

echo "gtgvcf `date`" >> finished
fi

##cmd:hardfilt
if [ "`echo $cmd|sed "s/ /\n/g"|awk '$1=="hardfilt"{print "on"}'`" = "on" ];then
sed -i "/^hardfilt/d" finished
runlog "running:calling:select-snp"
outdir="3_variants"
java -jar $soft_path/GATK/3.7/GenomeAnalysisTK.jar -T SelectVariants -R $ref_genome -selectType SNP --variant $outdir/$sample.raw.vcf -o $outdir/$sample.raw.snp.vcf

runlog "running:calling:filt-snp"
java -jar $soft_path/GATK/3.7/GenomeAnalysisTK.jar -T VariantFiltration -R $ref_genome --variant $outdir/$sample.raw.snp.vcf -o $outdir/$sample.raw.snp.fil.vcf --filterName "REJECT" --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" #将SNP选出来，然后根据指标过滤SNP

runlog "running:calling:select-indel"
java -jar $soft_path/GATK/3.7/GenomeAnalysisTK.jar -T SelectVariants -R $ref_genome -selectType INDEL --variant $outdir/$sample.raw.vcf -o $outdir/$sample.raw.indel.vcf

runlog "running:calling:filt-indel"
java -jar $soft_path/GATK/3.7/GenomeAnalysisTK.jar -T VariantFiltration -R $ref_genome --variant $outdir/$sample.raw.indel.vcf -o $outdir/$sample.raw.indel.fil.vcf --filterName "REJECT" --filterExpression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" #将indel选出来，根据指标过滤indel

runlog "running:calling:select-passed-snp/indel"
awk '/^#/ || $7 == "PASS"' $outdir/$sample.raw.snp.fil.vcf >$outdir/$sample.pass.snp.vcf
awk '/^#/ || $7 == "PASS"' $outdir/$sample.raw.indel.fil.vcf >$outdir/$sample.pass.indel.vcf

runlog "running:calling:catvariants"
java -cp $soft_path/GATK/3.7/GenomeAnalysisTK.jar org.broadinstitute.gatk.tools.CatVariants -R $ref_genome -V $outdir/$sample.pass.snp.vcf -V $outdir/$sample.pass.indel.vcf -out $outdir/$sample.mul.vcf -assumeSorted

#runlog "running:left-normalize"
##bcftools norm -m - $outdir/$sample.vcf | bcftools norm -f $ref_genome >$outdir/$sample.vcf
#bcftools norm -m - $outdir/$sample.mul.vcf | bcftools norm -f $ref_genome > $outdir/$sample.vcf

echo "hardfilt `date`" >> finished
fi

##cmd:vqsr
if [ "`echo $cmd|sed "s/ /\n/g"|awk '$1=="vqsr"{print "on"}'`" = "on" ];then
sed -i "/^vqsr/d" finished
runlog "running:vqsr" # 使用高斯模型，根据已有的真实变异位点（HapMap，omni...）来训练，最后的到一个训练好的能够评估真伪的错误评估模型，评估每一个变异位点发生错误的概率，给出一个得分，就是在训练好的混合高斯模型下，一个位点是真实的概率比上这个位点是假阳性的概率的log值，这个值越大就越好；这个模型首先要拿到真实变异数据集和待过滤的原始变异数据集的交集，然后对这些变异相对于具体注释信息的分布情况进行模拟，将这些变异位点进行聚类，最后根据聚类结果赋予所有变异位点相应的log值，越接近聚类核心的变异位点得到的值越高。设置阈值，log值高于阈值被认为是可信的，低于阈值的就会被过滤掉，设置阈值要兼顾敏感度和质量值
outdir="3_variants/vqsr"
[ ! -e $outdir ] && mkdir -p $outdir

#SNP
resource_text="--resource $ref_vqsr_hapmap --resource_param hapmap,known=false,training=true,truth=true,prior=15.0 "
resource_text="$resource_text --resource $ref_vqsr_omni --resource_param omni,known=false,training=true,truth=false,prior=12.0 "
resource_text="$resource_text --resource $ref_vqsr_1000G --resource_param 1000G,known=false,training=true,truth=false,prior=10.0 "
resource_text="$resource_text --resource $ref_snp --resource_param dbsnp,known=true,training=false,truth=false,prior=2.0"

annotation_array="QD MQ MQRankSum ReadPosRankSum FS SOR"
for annotation in $annotation_array; do 
	annotate_text="$annotate_text --annotation $annotation"
done

sentieon driver -r $ref_genome --algo VarCal -v 3_variants/$sample.raw.vcf $resource_text $annotate_text --max_gaussians 4 --var_type SNP --plot_file $outdir/snp.plots.R \
	--tranches_file $outdir/snp.tranches $outdir/snp.recal
	
sentieon driver -r $ref_genome --algo ApplyVarCal -v 3_variants/$sample.raw.vcf --sensitivity 99.0 --var_type SNP --tranches_file $outdir/snp.tranches --recal $outdir/snp.recal $outdir/$sample.raw.vqsr.snp.vcf

sentieon plot vqsr -o $outdir/snp.VQSR.pdf $outdir/snp.plots.R 

#INDEL
resource_text="--resource $ref_1000M_indel --resource_param mills,known=false,training=true,truth=true,prior=12.0 "
resource_text="$resource_text --resource $ref_snp --resource_param dbsnp,known=true,training=false,truth=false,prior=2.0"

annotate_text=""
annotation_array="QD DP MQRankSum ReadPosRankSum FS SOR"
for annotation in $annotation_array; do 
	annotate_text="$annotate_text --annotation $annotation"
done

sentieon driver -r $ref_genome --algo VarCal -v 3_variants/$sample.raw.vcf $resource_text $annotate_text --max_gaussians 4 --var_type INDEL --plot_file $outdir/indel.plots.R \
	--tranche 100.0 --tranche 99.9 --tranche 99.0 --tranche 90.0 \
	--tranches_file $outdir/indel.tranches $outdir/indel.recal
	
sentieon driver -r $ref_genome --algo ApplyVarCal -v 3_variants/$sample.raw.vcf --sensitivity 99.0 --var_type INDEL --tranches_file $outdir/indel.tranches --recal $outdir/indel.recal $outdir/$sample.raw.vqsr.indel.vcf

sentieon plot vqsr -o $outdir/indel.VQSR.pdf $outdir/indel.plots.R 
 
#filt pass:
awk -F "\t" '$1~/^#/ ||$7=="PASS" {print}' $outdir/$sample.raw.vqsr.snp.vcf > $outdir/$sample.raw.vqsr.snp.pass.vcf
awk -F "\t" '$1~/^#/ ||$7=="PASS" {print}' $outdir/$sample.raw.vqsr.indel.vcf > $outdir/$sample.raw.vqsr.indel.pass.vcf

#combine:
java -cp $soft_path/GATK/3.7/GenomeAnalysisTK.jar org.broadinstitute.gatk.tools.CatVariants -R $ref_genome -V $outdir/$sample.raw.vqsr.snp.pass.vcf -V $outdir/$sample.raw.vqsr.indel.pass.vcf -out 3_variants/$sample.mul.vcf -assumeSorted

#runlog "running:left-normalize"
#bcftools norm -m - 3_variants/$sample.mul.vcf | bcftools norm -f $ref_genome > 3_variants/$sample.vcf
echo "vqsr `date`" >> finished
fi

##cmd:left-normalize
runlog "running:left-normalize"
#如果不进行归一化，相同的变异可能会因为形式的不同而导致无法比较，例如：
#    #CHROM    POS    ID    REF    ALT
#    1         900010 .     GC     GCC
#    1         900010 .     G      GC
#这两个位点其实是相同的位点，如果不进行处理，对后续的注释和分析会带来很大的麻烦
#但是怎么确定是G后面插入了一个C，不是C后面插入了一个C？
# indel的左对齐
#    #CHROM    POS        ID    REF    ALT
#    1         900010     .     GC     GCC
#    1         12255506   .     TCCCCC TCCC
#通过左对齐后可以变成
#    #CHROM    POS        ID    REF    ALT
#    1         900010     .     G      GC 
#    1         12255506   .     TCC    T

bcftools norm -m - 3_variants/$sample.mul.vcf | bcftools norm -f $ref_genome > 3_variants/$sample.vcf
echo "left-normalize `date`" >> finished

##cmd:annotation
if [ "`echo $cmd|sed "s/ /\n/g"|awk '$1=="annotation"{print "on"}'`" = "on" ];then
sed -i "/^annotation/d" finished
runlog "running:annotation" # 使用vep、annovar和CADD进行注释
outdir="3_variants"
check_file $outdir/$sample.vcf
perl /DATA/ypliu/opt/mskcc-vcf2maf-1b16b35/vcf2maf.pl --input-vcf $outdir/$sample.vcf --output-maf $outdir/$sample.maf --vep-path /DATA/ypliu/opt/mskcc-vcf2maf-1b16b35/vep --vep-data /DATA/ypliu/opt/mskcc-vcf2maf-1b16b35/.vep --ref-fasta /DATA/sslyu/refGene/hg19.fa --filter-vcf /DATA/ypliu/opt/mskcc-vcf2maf-1b16b35/.vep/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz
sed -i '1d' $outdir/$sample.maf
/SSD750/PB1/soft1/annovar/convert2annovar.pl -format vcf4 $outdir/$sample.vcf -outfile $outdir/$sample.avinput
/SSD750/PB1/soft1/annovar/convert2annovar.pl -format vcf4 $outdir/$sample.vcf -outfile $outdir/$sample.link --includeinfo
sed 's/^chr//' $outdir/$sample.vcf|gzip -c ->$outdir/$sample.vcf.gz
score.sh $outdir/$sample.vcf.gz $outdir/$sample.CADD.gz
bgzip -df $outdir/$sample.CADD.gz
sed -i '1d' $outdir/$sample.CADD
#/SSD750/PB1/soft1/annovar/table_annovar.pl $outdir/$sample.avinput $ref_annotation -buildver hg19 -out $outdir/$sample.ann -protocol refGene,avsnp147,1000g2015aug_all,esp6500_all,exac03,dbnsfp33a,clinvar_20170130,ensembl,hgmd -operation g,f,f,f,f,f,f,f,f -remove -nastring "NA" -otherinfo
/SSD750/PB1/soft1/annovar/table_annovar.pl $outdir/$sample.avinput $ref_annotation -buildver hg19 -out $outdir/$sample.ann -protocol refGene,avsnp147,1000g2015aug_all,1000g2015aug_eas,gnomad_genome_eas,gnomad_exomes_hom,gnomad_exomes_hom_all,gnomad_exomes_hemi,exac03_eas,esp6500_all,dbscsnv11,dbnsfp33a,revel,intervar_20170202,clinvar_20170905,ensembl,hgmd_2018_spring,bed,bed,bed,bed -operation g,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,r,r,r,r -bedfile hg19_Dead_Zone.txt,hg19_Problem_High.txt,hg19_Problem_Low.txt,hg19_HI_Predictions_Version3.bed -arg ',,,,,,,,,,,,,,,,,-colsWanted 4,-colsWanted 4,-colsWanted 4,-colsWanted 4' -remove -nastring "NA" -otherinfo
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
convert2annovar.pl -format vcf4 3_variants/$sample.vcf -outfile tmprun/$sample.avinput
table_annovar.pl tmprun/$sample.avinput $ref_annotation -buildver hg19 -out tmprun/$sample -protocol refGene -operation g -remove

awk -F "\t" 'BEGIN{indel=0}  $4=="-"||$5=="-"{indel++} $4!="-"&&$5!="-"{print $4$5} END{print indel"\t"NR-1-indel"\t"}' tmprun/$sample.hg19_multianno.txt | awk 'BEGIN{ts=0} $1=="AG"||$1=="GA"||$1=="CT"||$1=="TC"{ts++} END{printf $0;tv=NR-2-ts;printf (ts+0.0)/tv"\t"}' >> 3_variants/$sample.qcvcf

awk -F "\t" '$6=="exonic"||$6=="exonic;splicing"{print $4$5}' tmprun/$sample.hg19_multianno.txt | awk '$1!~/-/{print}'| awk 'BEGIN{ts=0} $1=="AG"||$1=="GA"||$1=="CT"||$1=="TC"{ts++} END{tv=NR-1-ts;printf ts+tv"\t"(ts+0.0)/tv}' >> 3_variants/$sample.qcvcf

rm -rf tmprun
echo "qcvcf `date`" >> finished
fi
