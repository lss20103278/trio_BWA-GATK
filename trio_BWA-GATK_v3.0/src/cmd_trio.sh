#!/bin/sh

#SBATCH -J jobname
#SBATCH -p BP10
#SBATCH -N 1
#SBATCH -n 4

sample=list

sh run2_sentieon.sh $sample $SLURM_NPROCS
#samtools depth -b /DATA/sslyu/Project/Genetics_children_hospital/all_ch_freq_gene.bed /BP12_share/sslyu/bam/$sample.sort.mkdup.bam > ../depth/$sample.depth
scp /BP12_share/sslyu/bam/$sample.sort.mkdup.bam* klyang@122.112.248.194:/media/bpdata/data/exon_temp2/
curl 'https://oapi.dingtalk.com/robot/send?access_token=dc5bf30a873950e8c73bb37abf1851ebaf5d83f2792f08684ace2be621a12563' -H 'Content-Type: application/json' -d '{"msgtype": "text","text": {"content": "$sample mapping done"}}'
