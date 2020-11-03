#!/usr/bin/env python
# -*- coding: utf-8 -*-
#########################################################################
# File Name: /DATA/sslyu/trio_BWA-GATK_v3.0/create_email.py
# Author: sslyu
# mail: 1573077420@qq.com
# Created Time:Mon 22 Apr 2019 05:48:24 PM CST
#########################################################################

import numpy as np
import pandas as pd
import sys
reload(sys)
sys.setdefaultencoding("utf8")
import os
import smtplib
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from email.mime.base import MIMEBase
from email import encoders

if not os.path.exists('config'):
    print """
    Examples: 
    prepare config
    format: (delimiter \\t)
    -sn excelfile 
    -panel   IDT-PANEL 
    -filtmode   hard 
    -rawdir    absolute_path_of_rawdata 
    -md5    absolute_path_of_rawdata_of_md5_file 
    -R1 suffix_of_raw_R1 
    -R2 suffix_of_raw_R2
    -outdir absolute_path_of_backup
    -analysis_data the dirname of the analysis data of the rawdata
    -redo yes
    Note:
    make sure excel colnames contain 批次   收样日期    样本编号    原始样本ID    项目    性别    年龄    科室    浓度(ng/ul)    体积(ul)   总量(ug)    OD260/OD280 OD260/OD230 质检结果    上机日期    下机日期    Reads（M）   bases(Gb)  Q30
    the rawdata and md5 fileare in the same single dir 
    """
    sys.exit()

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

def create_email(n,panel,analysis_data):
    n = str(n)
    if os.getcwd().startswith('/DATA'):
        dirname = "/".join(os.getcwd().split('/')[-2:])
        dirname = "/DATA/sslyu/Project/Genetics_children_hospital/2019/"+dirname
    else:
        dirname = os.getcwd()
    
    if 'SR' in os.getcwd().split('/')[-1]:
        sn = '儿童遗传病科研全外分析'
    else:
        sn = os.getcwd().split('/')[-1].split('_')[1]
        batch = os.getcwd().split('/')[-1].split('_')[0]
        if len(os.getcwd().split('/')[-1].split('_')) == 3:
            #panel = os.getcwd().split('/')[-1].split('_')[2]
            if panel == 'IDT-PANEL':
                sn = sn+'儿医-'+n+'样-IDT4500基因检测分析（'+batch+'批）'
            if panel == 'WES':
                sn = sn+'儿医-'+n+'样-全外显子检测分析（'+batch+'批）'
            if panel == 'Agilent-wes':
                sn = sn+'儿医-'+n+'样-安捷伦亚全外显子检测分析（'+batch+'批）'
            if panel == 'T084V2':
                sn = sn+'儿医-'+n+'样-基因检测分析（'+batch+'批_'+panel+'）'
            if panel == 'T086V4':
                sn = sn+'儿医-'+n+'样-基因检测分析（'+batch+'批_'+panel+'）'
        elif len(os.getcwd().split('/')[-1].split('_')) == 4:
            if os.getcwd().split('/')[-1].split('_')[2] == "medical":
                panel = os.getcwd().split('/')[-1].split('_')[3]
                sn = sn+'儿医-'+n+'样-医学全外显子检测分析（'+batch+'批）'
            if panel == 'T084V2_CNV':
                sn = sn+'儿医-'+n+'样-基因检测分析（'+batch+'批_'+panel+'）'
    ofile = open('email.txt','w')
    ofile.write('以下为儿童遗传病分析结果：\n批次：'+sn+'结果\n\n位置\n'+dirname+'\n\n')
    if os.path.exists('peddy_check.log'):
        ofile.write('备注：\n')
        with open('peddy_check.log') as f:
            for l in f:
                ofile.write(l)
        ofile.write('\n')
    if os.path.exists('qcsum_check.log'):
        ofile.write('备注：\n')
        with open('qcsum_check.log') as f:
            for l in f:
                ofile.write(l)
        ofile.write('\n')                
    if os.path.exists('CNV.log'):
        ofile.write('备注：\n')
        with open('CNV.log') as f:
            for l in f:
                ofile.write(l)
        ofile.write('\n')                
    ofile.write('bed为'+panel+'\n')
    ofile.write('原始数据位于  /anaData/anaData004/children_hos_genetic/rawdata/2019/'+panel+'\n')
    date = os.getcwd().split('/')[-1].split('_')[1]
    batch = os.getcwd().split('/')[-1].split('_')[0]                             
    ofile.write('低通量原始数据位于 /anaData/anaData004/children_hos_genetic/rawdata/2019/LOWWGS/'+batch+'\n')
    ofile.write('原始数据的分析数据位于   '+analysis_data+'\n')
    ofile.write('bam文件上传地址    http://122.112.248.194/GMOD/BP_yun/?data=sample_data%2Fjson%2FBPdata%2Fexon_temp2%2Fdata&loc=chrX%3A70334659..70362518&tracks=DNA%2CUCSC_hg19.gff3&highlight=\n')
    ofile.write('祝好，\n吕珊珊\n\n')                                                                                           
    ofile.close()


if sn.endswith('xlsx'):
    excel = pd.read_excel(sn, dtype=str)
elif sn.endswith('csv'):
    excel = pd.read_csv(sn, dtype=str, encoding="utf-8", sep="\t")
elif sn.endswith('txt'):
    excel = pd.read_csv(sn, dtype=str, encoding="utf-8", sep="\t")
else:
    print "the format of the input file is wrong ..."
    sys.exit()

n = excel.shape[0]
create_email(n, panel,analysis_data)
os.system('cat email.txt')

