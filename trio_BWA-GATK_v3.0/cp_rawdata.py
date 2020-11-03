#!/usr/bin/env python
# -*- coding: utf-8 -*-
#########################################################################
# File Name: cp_rawdata.py
# Author: sslyu
# mail: 1573077420@qq.com
# Created Time:Tue 15 Jan 2019 06:26:43 PM CST
#########################################################################

import numpy as np
import pandas as pd
import sys
reload(sys)
sys.setdefaultencoding("utf8")
import os

def generate_cp_redo(sample,ID,rawdir,R1,R2):
    if len(os.listdir(rawdir+'/'+ID)) > 3:
        print "there are more than two fastq files of %s" % ID
        sys.exit()
    cp = open(sample+'/cp_rawdata.sh', 'w')
    cp.write('#!/bin/sh\n\n#SBATCH -J '+sample+'\n#SBATCH -p BP10\n#SBATCH -N 1\n#SBATCH -n 4\n#SBATCH --mail-type=FAIL\n#SBATCH --mail-user=1573077420@qq.com\n\n')
    cp.write('sample='+sample+'\nID='+ID+'\n\ndir='+rawdir+'\nR1='+R1+'\nR2='+R2+'\n\n')
    cp.write('cp $dir/$ID/*$ID*$R1 ../raw/$sample\_R1.fastq.gz\n')
    cp.write('cp $dir/$ID/*$ID*$R2 ../raw/$sample\_R2.fastq.gz\n')
    cp.write('cat $dir/$ID/$ID.md5 >> ../backup_md5\n')
    cp.write('for j in $R1 $R2; do md5sum $dir/$ID/*$ID*$j >> ../original_md5; done\n')
    cp.write('for j in R1 R2; do md5sum ../raw/$sample\_$j*.fastq.gz >> ../cp_md5; done\n')
    cp.write('for j in $R1 $R2; do original=`grep $ID ../original_md5 |grep $j |cut -d" " -f 1`; backup=`grep $ID ../backup_md5 |grep $j |cut -d" " -f 1`; if [ "$original" != "$backup" ]; then echo -e $ID"_"$j" is not backup correctly" >> ../check_md5.log; echo -e $ID"_"$j" is not backup correctly"; fi; done\n')
    cp.write('pattern_original=($R1 $R2); pattern_copied=(R1 R2); for j in 0 1; do original=`grep $ID ../original_md5 |grep ${pattern_original[$j]} |cut -d" " -f 1`; copied=`grep $sample ../cp_md5 |grep ${pattern_copied[$j]} |cut -d" " -f 1`; if [ "$original" != "$copied" ]; then echo -e $sample"_"$ID"_"${pattern_original[$j]}" is not copied correctly" >> ../check_md5.log; echo -e $sample"_"$ID"_"${pattern_original[$j]}" is not copied correctly"; fi; done\n')
    cp.write('curl \'https://oapi.dingtalk.com/robot/send?access_token=dc5bf30a873950e8c73bb37abf1851ebaf5d83f2792f08684ace2be621a12563\' -H \'Content-Type: application/json\' -d \'{"msgtype": "text","text": {"content": "'+sample+' copy rawdata done"}}\'\n')                                                                                                                              
    cp.close()                                

def generate_cp(sample,ID,rawdir,md5,R1,R2,panel,outdir):
    cp = open(sample+'/cp_rawdata.sh', 'w')                                                                                       
    cp.write('#!/bin/sh\n\n#SBATCH -J '+sample+'\n#SBATCH -p BP10\n#SBATCH -N 1\n#SBATCH -n 4\n#SBATCH --mail-type=FAIL\n#SBATCH --mail-user=1573077420@qq.com\n\n') #--mail-type=END,FAIL
    cp.write('sample='+sample+'\nID='+ID+'\n\ndir='+rawdir+'\nmd5='+md5+'\nR1='+R1+'\nR2='+R2+'\npanel='+panel+'\n\n')
    cp.write('outdir='+outdir+'\n')
    #if os.path.exists(outdir+'/'+ID):
    #    print "%s/%s exists ..." % (outdir, ID)
    #    sys.exit()
    cp.write('mkdir -p $outdir/$ID\n')
    cp.write('awk \'{if(match($2,"\'$ID\'")){print $0}}\' $md5 > $outdir/$ID/$ID.md5\n')
    cp.write('cp $dir/*$ID*$R1 $outdir/$ID\n')
    cp.write('cp $dir/*$ID*$R2 $outdir/$ID\n')
    cp.write('cp $dir/*$ID*$R1 ../raw/$sample\_R1.fastq.gz\n')
    cp.write('cp $dir/*$ID*$R2 ../raw/$sample\_R2.fastq.gz\n')

    cp.write('for j in $R1 $R2; do md5sum $dir/*$ID*$j >> ../original_md5; done\n')
    cp.write('for j in R1 R2; do md5sum ../raw/$sample\_$j*.fastq.gz >> ../cp_md5; done\n')
    cp.write('for j in $R1 $R2; do md5sum $outdir/$ID/*$ID*$j >> ../backup_md5; done\n')
    cp.write('for j in $R1 $R2; do a=`grep $ID ../original_md5 |grep $j |cut -d" " -f 1`; b=`grep $ID $md5 |grep $j |cut -d" " -f 1`; if [ "$a" != "$b" ]; then echo -e "original "$ID"_"$j" is problematic" >> ../check_md5.log; echo -e "original "$ID"_"$j" is problematic"; fi; done\n')
    cp.write('pattern_original=($R1 $R2); pattern_copied=(R1 R2); for j in 0 1; do original=`grep $ID ../original_md5 |grep ${pattern_original[$j]} |cut -d" " -f 1`; copied=`grep $sample ../cp_md5 |grep ${pattern_copied[$j]} |cut -d" " -f 1`; if [ "$original" != "$copied" ]; then echo -e $sample"_"$ID"_"${pattern_original[$j]}" is not copied correctly" >> ../check_md5.log; echo -e $sample"_"$ID"_"${pattern_original[$j]}" is not copied correctly"; fi; done\n')
    cp.write('for j in $R1 $R2; do original=`grep $ID ../original_md5 |grep $j |cut -d" " -f 1`; backup=`grep $ID ../backup_md5 |grep $j |cut -d" " -f 1`; if [ "$original" != "$backup" ]; then echo -e $ID"_"$j" is not backup correctly" >> ../check_md5.log; echo -e $ID"_"$j" is not backup correctly"; fi; done\n')
    cp.write('curl \'https://oapi.dingtalk.com/robot/send?access_token=dc5bf30a873950e8c73bb37abf1851ebaf5d83f2792f08684ace2be621a12563\' -H \'Content-Type: application/json\' -d \'{"msgtype": "text","text": {"content": "'+sample+' copy rawdata done"}}\'\n')  
    cp.close()

kwargs={'-sn':'', '-panel':'', '-filtmode':'', '-rawdir':'', '-md5':'', '-R1':'', '-R2':'', '-outdir':'', '-analysis_data':'', '-redo':''}

with open('config') as f:
    for l in f:
        l = l.strip('\n').split('\t')
        if len(l) != 0:
            for j in kwargs.keys():
                if l[0] == j:
                    kwargs.update({j:l[1]})

sn=kwargs.get('-sn')
panel=kwargs.get('-panel')
filtmode=kwargs.get('-filtmode')
rawdir = kwargs.get('-rawdir')
md5 = kwargs.get('-md5')
R1 = kwargs.get('-R1')
R2 = kwargs.get('-R2')
outdir = kwargs.get('-outdir')
analysis_data = kwargs.get('-analysis_data')
redo = kwargs.get('-redo')
workpath = os.getcwd().split('/')[-2:]
if u'分析任务单' in sn:
    sep = u'分析任务单'
else:
    sep = u'产品信息表'

if not os.path.exists('raw'):
    os.mkdir('raw')

if not os.path.exists('rename'):
    if sn.endswith('xlsx'):
        d = pd.read_excel(sn, dtype=str)
    elif sn.endswith('csv'):
        d = pd.read_csv(sn, dtype=str, encoding="utf-8", sep="\t")
    elif sn.endswith('txt'):
        d = pd.read_csv(sn, dtype=str, encoding="utf-8", sep="\t")
    else:
        print "the format of the input file is wrong ..."
        sys.exit()
    if 'sample' in d.columns and 'Sample_ID' in d.columns:
        d[['sample', 'Sample_ID']].to_csv('rename', index=False, encoding="utf-8", sep="\t")
    elif u'样本编号' in d.columns and u'原始样本ID' in d.columns:
        d.rename(columns={u'样本编号':'Sample_ID',u'原始样本ID':'sample'}, inplace=True)
        d[['sample', 'Sample_ID']].to_csv('rename', index=False, encoding="utf-8", sep="\t")
    else:
        print "lack sample column and Sample_ID column"
        sys.exit()

excel = pd.read_csv('rename', encoding="utf-8", sep="\t")  # rename header
list_ID = open('list', 'w')
for i in excel.index:
    if excel.isnull().loc[i,'sample']:
        ID = excel.loc[i]['Sample_ID']
        ID = str(ID)
        sample = ID
    else:
        sample = excel.loc[i]['sample']
        sample = str(sample)
        ID = excel.loc[i]['Sample_ID']
        ID = str(ID)
    list_ID.write(sample+'\n')
    if not os.path.exists(sample):
        os.mkdir(sample)
    if os.path.exists(sample+'/cmd.sh'):
        os.remove(sample+'/cmd.sh')
    if redo == 'yes':
        generate_cp_redo(sample, ID, rawdir, R1, R2)
    else:
        generate_cp(sample, ID, rawdir, md5, R1, R2, panel, outdir)
list_ID.close()

os.system('for i in `cat list`; do cd $i; sbatch cp_rawdata.sh; cd ..; done')        
