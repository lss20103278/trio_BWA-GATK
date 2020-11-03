#########################################################################
# File Name: src/unmapped_duplicate.sh
# Author: sslyu
# mail: 1573077420@qq.com
# Created Time:Thu 11 Oct 2018 11:08:40 PM CST
#########################################################################
#!/bin/bash

list=$1

echo -e "sample\tUNMAPPED_READS\tUNPAIRED_READ_DUPLICATES\tREAD_PAIR_DUPLICATES\tPERCENT_DUPLICATION" > unmapped_duplicate
for i in `cat $list`
do
	sed -n '2p' $i/2_mapping/$i.metrics |awk -F'\t' '{OFS="\t";print "'$i'",$5,$6,$7,$9}' >> unmapped_duplicate
done
