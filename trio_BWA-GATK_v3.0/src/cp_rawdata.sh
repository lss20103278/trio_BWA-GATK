#!/bin/sh

#SBATCH -J 19p4126
#SBATCH -p BP10
#SBATCH -N 1
#SBATCH -n 4

sample=19p4126
ID=B190105109

dir=/DATA/rawdata/sslv/98_eryi_Panel/
md5=/DATA/rawdata/sslv/98_eryi_Panel/md5sum.txt
R1=R1_001.fastq.gz
R2=R2_001.fastq.gz
strategy=IDT-PANEL

mkdir -p /anaData/anaData004/children_hos_genetic/rawdata/2018/$strategy/$ID
awk '{if(match($2,"'$ID'")){print $0}}' $md5 >> /anaData/anaData004/children_hos_genetic/rawdata/2018/$strategy/$ID/$ID.md5
cp $dir/*$ID*gz /anaData/anaData004/children_hos_genetic/rawdata/2018/$strategy/$ID
cp $dir/*$ID*$R1 ../raw/$sample\_R1.fastq.gz
cp $dir/*$ID*$R2 ../raw/$sample\_R2.fastq.gz
for j in $R1 $R2; do md5sum $dir/*$ID*$j >> ../original_md5; done
for j in R1 R2; do md5sum ../raw/$sample\_$j.fastq.gz >> ../cp_md5; done
for j in $R1 $R2; do md5sum /anaData/anaData004/children_hos_genetic/rawdata/2018/$strategy/$ID/*$ID*$j >> ../backup_md5; done
for j in $R1 $R2; do a=`grep $ID ../original_md5 |grep $j |cut -d" " -f 1`; b=`grep $ID $md5 |grep $j |cut -d" " -f 1`; if [ "$a" != "$b" ]; then echo -e "original "$ID"_"$j" is problematic" >> ../check_md5.log; echo -e "original "$ID"_"$j" is problematic"; fi; done
pattern_original=($R1 $R2); pattern_copied=(R1 R2); for j in 0 1; do original=`grep $ID ../original_md5 |grep ${pattern_original[$j]} |cut -d" " -f 1`; copied=`grep $sample ../cp_md5 |grep ${pattern_copied[$j]} |cut -d" " -f 1`; if [ "$original" != "$copied" ]; then echo -e $ID"_"${pattern_original[$j]}" is not copied correctly" >> ../check_md5.log; echo -e $ID"_"${pattern_original[$j]}" is not copied correctly"; fi; done
for j in $R1 $R2; do original=`grep $ID ../original_md5 |grep $j |cut -d" " -f 1`; backup=`grep $ID ../backup_md5 |grep $j |cut -d" " -f 1`; if [ "$original" != "$backup" ]; then echo -e $ID"_"$j" is not backup correctly" >> ../check_md5.log; echo -e $ID"_"$j" is not backup correctly"; fi; done
curl 'https://oapi.dingtalk.com/robot/send?access_token=dc5bf30a873950e8c73bb37abf1851ebaf5d83f2792f08684ace2be621a12563' -H 'Content-Type: application/json' -d '{"msgtype": "text","text": {"content": "$sample copy rawdata done"}}'
