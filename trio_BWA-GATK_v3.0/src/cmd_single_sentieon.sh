#!/bin/sh

#SBATCH -J jobname
#SBATCH -p BP10
#SBATCH -N 1
#SBATCH -n 4

sample=list

sh run2_sentieon.sh $sample $SLURM_NPROCS
#samtools depth -b /DATA/sslyu/Project/Genetics_children_hospital/all_ch_freq_gene.bed /BP12_share/sslyu/bam/$sample.sort.mkdup.bam > ../depth/$sample.depth
scp /BP12_share/sslyu/bam/$sample.sort.mkdup.bam* klyang@122.112.248.194:/media/bpdata/data/exon_temp2/
[ ! -e ../annotation ] && mkdir -p ../annotation
cp 2_mapping/$sample.depth.sample_gene_summary ../annotation
cp 3_variants/$sample.{ann.hg19_multianno.txt,CADD,maf,link} ../annotation
cd ../annotation
echo $sample >> list
python /DATA/sslyu/trio_BWA-GATK_3.0_20190122/src/score_re.py -sn $sample
python /DATA/sslyu/trio_BWA-GATK_3.0_20190122/src/annotation_filt.py -sn $sample -c_f_m c
cp $sample\_children_hospital*xlsx ../98_20190105_IDT-PANEL
cp ../$sample/3_variants/$sample.vcf ../98_20190105_IDT-PANEL/vcf
echo "end `date`" >> ../$sample/finished
curl 'https://oapi.dingtalk.com/robot/send?access_token=dc5bf30a873950e8c73bb37abf1851ebaf5d83f2792f08684ace2be621a12563' -H 'Content-Type: application/json' -d '{"msgtype": "text","text": {"content": "'"$sample"' annotation excel done"}}'
