#!/usr/bin/env python
# -*- coding: utf-8 -*-
#########################################################################
# File Name: CNV.py
# Author: sslyu
# mail: 1573077420@qq.com
# Created Time:Tue 16 Apr 2019 07:52:56 PM CST
#########################################################################

import numpy as np
import pandas as pd
import sys
reload(sys)
sys.setdefaultencoding("utf8")
import os

maindir="/DATA/sslyu/trio_BWA-GATK_v3.0"

def CNV(panel):
    CNV_read = open('CNV/CNV_Read.R','w')
    CNV_read.write('.libPaths("/home/ana005/R/x86_64-pc-linux-gnu-library/3.4")\nlibrary(ExomeDepth)\n')
    print """
.libPaths("/home/ana005/R/x86_64-pc-linux-gnu-library/3.4")
library(ExomeDepth)
    """
    if 'IDT-PANEL' == panel:
        CNV_read.write('read.table("'+maindir+'/bed/IDT_gene.bed")->bed\n')
        print """
read.table("%s/bed/IDT_gene.bed")->bed
        """ % maindir
    elif 'WES' == panel:
        CNV_read.write('read.table("'+maindir+'/bed/wes_gene.bed")->bed\n')
        print """
read.table("%s/bed/wes_gene.bed")->bed
        """ % maindir
    elif 'T084V2' == panel:
        CNV_read.write('read.table("'+maindir+'/bed/T084V2.ig.cnv.bed")->bed\n')
        print """
read.table("%s/bed/T084V2.ig.cnv.bed")->bed
        """ % maindir
    elif 'T084V2_CNV' == panel:
        CNV_read.write('read.table("'+maindir+'/bed/T084V2_CNV.igt.cnv.bed")->bed\n')
        print """
read.table("%s/bed/T084V2_CNV.igt.cnv.bed")->bed
        """ % maindir
    elif 'T086V4' == panel:
        CNV_read.write('read.table("'+maindir+'/bed/T086V4.igt.cnv.bed")->bed\n')
        print """
read.table("%s/bed/T086V4.igt.cnv.bed")->bed
        """ % maindir
    elif 'Agilent-wes' == panel:
        CNV_read.write('read.table("'+maindir+'/doc/S06588914_Regions_cnv.bed")->bed\n')
        print """
read.table("%s/doc/S06588914_Regions_cnv.bed")->bed
        """ % maindir
    else:
        print 'bed is lost...'
        sys.exit()
    print """
colnames(bed)=c("chromosome","start","end","name")
read.table("list")->bamfile
unlist(bamfile)->bamfile
bamfile=as.character(bamfile)
bamfile <- as.vector(bamfile)
bam.counts <- getBamCounts(bed.frame=bed, bam.files=bamfile, referenceFasta="/SSD750/PB1/db1/Homo/refseq/hg19.fa")
save.image("bam.counts.Rdata")
q()
    """
    CNV_read.write('colnames(bed)=c("chromosome","start","end","name")\nread.table("list")->bamfile\nunlist(bamfile)->bamfile\nbamfile=as.character(bamfile)\nbamfile <- as.vector(bamfile)\nbam.counts <- getBamCounts(bed.frame=bed, bam.files=bamfile, referenceFasta="/SSD750/PB1/db1/Homo/refseq/hg19.fa")\nsave.image("bam.counts.Rdata")\nq()\n')
    CNV_read.close()
    CNV_bash = open('CNV/CNV.bash', 'w')
    CNV_bash.write('cd CNV/\n')
    CNV_bash.write('/DATA/ypliu/opt/R-3.4.3/bin/R # note: choose the right bed file\n')
    CNV_bash.write('/DATA/ypliu/opt/R-3.4.3/bin/R --slave --vanilla < '+maindir+'/src/CNV_Run.R > log 2>&1\n')
    CNV_bash.write('cd ..\n')
    CNV_bash.write('curl \'https://oapi.dingtalk.com/robot/send?access_token=dc5bf30a873950e8c73bb37abf1851ebaf5d83f2792f08684ace2be621a12563\' -H \'Content-Type: application/json\' -d \'{"msgtype": "text","text": {"content": "'+os.getcwd()+' cnv done"}}\'\n')
    CNV_bash.close()
    os.system('bash CNV/CNV.bash')

def cnv():
    d1 = pd.read_table('list', header=None)
    d1.rename(columns={0:'sample'}, inplace=True)    
    for i in d1['sample']:
        i = str(i)
        if not os.path.exists('/BP12_share/sslyu/bam/'+i+'.sort.mkdup.bam'):
            print '/BP12_share/sslyu/bam/'+i+'.sort.mkdup.bam doesn\'t exist'
            sys.exit()
    if not os.path.exists('CNV'):
        os.mkdir('CNV')
    dirname = os.getcwd().split('/')[-1]
    if not os.path.exists(dirname+'/CNV'):
        os.makedirs(dirname+'/CNV')
    if not os.path.exists('CNV/id') and not os.path.exists('CNV/list'):
        d1['sample'].to_csv('CNV/id', header=None, index=None)
        d1['CNV'] = '/BP12_share/sslyu/bam/'+d1['sample'].astype(str)+'.sort.mkdup.bam'
        d1['CNV'].to_csv('CNV/list', header=None, index=None)
    CNV(panel)

if not os.path.exists('config'):
    #shutil.copy(m_path+'/src/config', '.')
    print """
    Examples: 
    prepare config
    format: (delimiter \t)
    -sn excelfile 
    -panel   IDT-PANEL 
    -filtmode   hard 
    -rawdir    absolute_path_of_rawdata 
    -md5    absolute_path_of_rawdata_of_md5_file 
    -R1 suffix_of_raw_R1 
    -R2 suffix_of_raw_R2
    -outdir absolute_path_of_backup
    -analysis_data the dirname of the analysis data of the rawdata
    -redo   yes
    Note:
    make sure excel colnames contain 批次   收样日期    样本编号    原始样本ID    项目    性别    年龄    科室    浓度(ng/ul)    体积(ul)   总量(ug)    OD260/OD280 OD260/OD230 质检结果    上机日期    下机日期    Reads（M）   bases(Gb)  Q30
    the rawdata and md5 fileare in the same single dir 
    """
    sys.exit()

kwargs={'-sn':'', '-panel':'', '-filtmode':'', '-dir':'', '-md5':'', '-R1':'', '-R2':'', '-outdir':'', '-analysis_data':'', '-redo':''}

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
dir = kwargs.get('-dir')
md5 = kwargs.get('-md5')
R1 = kwargs.get('-R1')
R2 = kwargs.get('-R2')
outdir = kwargs.get('-outdir')
analysis_data = kwargs.get('-analysis_data')
redo = kwargs.get('-redo')

cnv()
os.system('grep -ci error CNV/log')
