#!/bin/sh

:<<Note_script_version
ver:2.1
for:
update:2017/9/19
by:ykl
Note_script_version

:<<HELP
triomode="gvcf | combine"
filtmode="vqsr | hardfilt"

pre-work:
./ped: prepare *.ped for each family
./gvcf: prepare *.g.vcf for each sample
HELP


sn=$1
proc=$2
filtmode=$3

if [ "$filtmode" != "hard" -a "$filtmode" != "vqsr" ];then
    echo "ERROR: filter method not specified or wrong type('hard' / 'vqsr' only). Exiting ..."
    exit
fi

##env:
source `pwd`/../../dbevn.sh




cp ../ped/$sn.ped ./
samples=`awk '{print $2}' $sn.ped`

#for i in $samples
#do
#[ ! -e ../gvcf/$i.raw.g.vcf ] && mv ../gvcf/$i.g.vcf ../gvcf/$i.raw.g.vcf
#./gvcf.sh $i
#done

echo "" > $sn.gvcf.list
for i in $samples
do
echo "../gvcf/$i.g.vcf" >> $sn.gvcf.list
done 
sed -i '/^$/d' $sn.gvcf.list



outdir="triotmp"
[ ! -d $outdir ] && mkdir $outdir
#combine gvcfs:
java -jar $soft_path/GATK/3.7/GenomeAnalysisTK.jar -T GenotypeGVCFs -R $ref_genome -L $panel --dbsnp $ref_snp \
 -stand_call_conf 10 \
 -o $outdir/$sn.raw.vcf \
 --variant $sn.gvcf.list
 

#vqsr:
sh vfilt.sh $outdir $sn $filtmode


#CalculateGenotypePosteriors:
java -jar $soft_path/GATK/3.7/GenomeAnalysisTK.jar -T CalculateGenotypePosteriors \
 -R $ref_genome \
 -V $outdir/$sn.raw.snp.pass.vcf \
 -ped ../ped/$sn.ped \
 -o $outdir/$sn.snp.gt.vcf 

java -jar $soft_path/GATK/3.7/GenomeAnalysisTK.jar -T CalculateGenotypePosteriors \
 -R $ref_genome \
 -V $outdir/$sn.raw.indel.pass.vcf \
 -ped ../ped/$sn.ped \
 -o $outdir/$sn.indel.gt.vcf

#lowGQ:
java -jar $soft_path/GATK/3.7/GenomeAnalysisTK.jar -T VariantFiltration \
 -R $ref_genome \
 --genotypeFilterExpression "DP<20 && GQ<20" \
 --genotypeFilterName "lowGQ" \
 -V $outdir/$sn.snp.gt.vcf \
 -o $outdir/$sn.snp.gt.flt.vcf

java -jar $soft_path/GATK/3.7/GenomeAnalysisTK.jar -T VariantFiltration \
 -R $ref_genome \
 --genotypeFilterExpression "DP<20 && GQ<20" \
 --genotypeFilterName "lowGQ" \
 -V $outdir/$sn.indel.gt.vcf \
 -o $outdir/$sn.indel.gt.flt.vcf

#denovo_anotation:
java -jar $soft_path/GATK/3.7/GenomeAnalysisTK.jar -T VariantAnnotator -L $panel \
 -R $ref_genome \
 -A PossibleDeNovo \
 -ped ../ped/$sn.ped \
 -V $outdir/$sn.snp.gt.flt.vcf \
 -o $outdir/$sn.snp.gt.novo.vcf

java -jar $soft_path/GATK/3.7/GenomeAnalysisTK.jar -T VariantAnnotator -L $panel \
 -R $ref_genome \
 -A PossibleDeNovo \
 -ped ../ped/$sn.ped \
 -V $outdir/$sn.indel.gt.flt.vcf \
 -o $outdir/$sn.indel.gt.novo.vcf

#combine
java -cp $soft_path/GATK/3.7/GenomeAnalysisTK.jar org.broadinstitute.gatk.tools.CatVariants -R $ref_genome -V $outdir/$sn.snp.gt.novo.vcf -V $outdir/$sn.indel.gt.novo.vcf -out $sn.gt.vcf -assumeSorted


grep -v ^## $sn.gt.vcf>$sn.tmp.vcf
grep ^## $sn.gt.vcf>$sn.head.vcf

./sort_sample.sh $sn
if [ -e $sn.sort.tmp.vcf ]
then
cat $sn.head.vcf $sn.sort.tmp.vcf>$sn.gt.vcf
else
cat $sn.head.vcf $sn.tmp.vcf>$sn.gt.vcf
fi

#vep-mendelscan-annovar
sh vep-mendelscan.sh $sn ../ped/$sn.mendel.ped



#sep:
outdir='sep'
[ ! -e $outdir ] && mkdir $outdir
for sepsn in $samples
do
#java -jar $soft_path/GenomeAnalysisTK.jar -T SelectVariants -R $ref_genome -L $panel -V $sn.gt.vcf -sn $sepsn -o $outdir/$sepsn.vcf

vcf-subset -c $sepsn $sn.gt.vcf > $outdir/$sepsn.vcf
done


#qcvcf:
echo -e "sample\tindel_genomic\tsnp_genomic\tts/tv_genomic\tsnp_exonic\tts/tv_exonic" > triotmp/$sn.qcvcf
[ ! -d tmprun ] && mkdir tmprun
for sepsn in $samples
do
printf "$sepsn\t" >> triotmp/$sn.qcvcf
convert2annovar.pl -format vcf4 sep/$sepsn.vcf -outfile tmprun/$sepsn.avinput
table_annovar.pl tmprun/$sepsn.avinput $ref_annotation -buildver hg19 -out tmprun/$sepsn -protocol refGene -operation g -remove

awk -F "\t" 'BEGIN{indel=0}  $4=="-"||$5=="-"{indel++} $4!="-"&&$5!="-"{print $4$5} END{print indel"\t"NR-1-indel"\t"}' tmprun/$sepsn.hg19_multianno.txt | awk 'BEGIN{ts=0} $1=="AG"||$1=="GA"||$1=="CT"||$1=="TC"{ts++} END{printf $0;tv=NR-2-ts;printf (ts+0.0)/tv"\t"}' >> triotmp/$sn.qcvcf

awk -F "\t" '$6=="exonic"||$6=="exonic;splicing"{print $4$5}' tmprun/$sepsn.hg19_multianno.txt | awk '$1!~/-/{print}'| awk 'BEGIN{ts=0} $1=="AG"||$1=="GA"||$1=="CT"||$1=="TC"{ts++} END{tv=NR-1-ts;printf ts+tv"\t"(ts+0.0)/tv"\n"}' >> triotmp/$sn.qcvcf
done
rm -rf tmprun



