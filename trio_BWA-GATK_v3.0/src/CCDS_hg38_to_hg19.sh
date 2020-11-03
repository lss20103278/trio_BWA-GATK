#!/bin/sh

#SBATCH -J jobname
#SBATCH -p BP10
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=1573077420@qq.com

#########################################################################
# File Name: CCDS_hg38_to_hg19.sh
# Author: sslyu
# mail: 1573077420@qq.com
# Created Time:Thu 11 Apr 2019 01:31:28 PM CST
# Usage:
#########################################################################

cp /DATA/sslyu/ensembl_human/CCDS.20180614.txt CCDS.20180614.hg38.txt
cut -f 1,3,8,9 CCDS.20180614.hg38.txt |awk '{OFS="\t";print "chr"$1,$3,$4,$2}' |sed '1d' |grep -v -w "-" > CCDS.20180614.hg38.bed
CrossMap.py 
