#!/bin/sh

#SBATCH -J peddy
#SBATCH -p BP10
#SBATCH -N 1
#SBATCH -n 4

sample=peddy


#select from 'hard'/'vqsr'
filtmode='vqsr'

basedir=/DATA/sslyu/trio_BWA-GATK_v3.0/src

path=`pwd`
var=`echo $path | awk -F '/' '{print $NF}'`
#echo $path
#if echo $path | grep peddy > /dev/null
if [ $var = 'peddy' ]
then
sh $basedir/trio_sentieon.sh $sample $SLURM_NPROCS
curl 'https://oapi.dingtalk.com/robot/send?access_token=dc5bf30a873950e8c73bb37abf1851ebaf5d83f2792f08684ace2be621a12563' -H 'Content-Type: application/json' -d '{"msgtype": "text","text": {"content": "peddy done"}}'
else
sh $basedir/trio_sentieon.sh $sample $SLURM_NPROCS $filtmode
[ -e ../../$sample/2_mapping/$sample.depth.sample_gene_summary ] && cp ../../$sample/2_mapping/$sample.depth.sample_gene_summary ../../annotation
cp segtrio/$sample.{ann.hg19_multianno.txt,CADD,maf,link} ../../annotation
cd ../../annotation
echo $sample >> list
python $basedir/score_re.py -sn $sample
python $basedir/annotation_filt.py -sn $sample -c_f_m c_f_m
echo "end `date`" >> ../trio/$sample/finished
curl 'https://oapi.dingtalk.com/robot/send?access_token=dc5bf30a873950e8c73bb37abf1851ebaf5d83f2792f08684ace2be621a12563' -H 'Content-Type: application/json' -d '{"msgtype": "text","text": {"content": "$sample annotation excel done"}}'
fi
