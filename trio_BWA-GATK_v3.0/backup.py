#!/usr/bin/env python
# -*- coding: utf-8 -*-
#########################################################################
# File Name: /DATA/sslyu/trio_BWA-GATK_v3.0/backup.py
# Author: sslyu
# mail: 1573077420@qq.com
# Created Time:Mon 22 Apr 2019 05:42:25 PM CST
#########################################################################

import numpy as np
import pandas as pd
import sys
reload(sys)
sys.setdefaultencoding("utf8")
import os

def backup(ofile, dirname, sn, sep, panel):
    #snfile = sn.split(sep)[0]+'交付信息表'+u'）'.join(sn.split(sep)[1].split(')'))
    ofile.write('cat CNV/log\n')
    ofile.write('grep -i error CNV/log\n')
    ofile.write('cp *交付信息表*xlsx '+dirname[-1]+'\n')
    ofile.write('cp CNV/*csv '+dirname[-1]+'/CNV\n')
    ofile.write('cp exon/*csv '+dirname[-1]+'\n')
    ofile.write('ls -R '+dirname[-1]+'\n\n')
    ofile.write('zip -r '+dirname[-1]+'.zip '+dirname[-1]+'\n')
    #ofile.write('if [ ! -e /DATA/BPshare/opm/遗传病-交大附属儿医/2018/'+dirname[-2]+' ]\n')
    ofile.write('if [ ! -e /DATA/BPshare/opm/遗传病-交大附属儿医/2019/'+dirname[-2]+' ]\n')
    ofile.write('then\n')
    #ofile.write('mkdir -p /DATA/BPshare/opm/遗传病-交大附属儿医/2018/'+dirname[-2]+'\n')
    ofile.write('mkdir -p /DATA/BPshare/opm/遗传病-交大附属儿医/2019/'+dirname[-2]+'\n')
    ofile.write('fi\n')
    #ofile.write('cp -r '+dirname[-1]+' /DATA/BPshare/opm/遗传病-交大附属儿医/2018/'+dirname[-2]+'\n')
    ofile.write('cp -r '+dirname[-1]+' /DATA/BPshare/opm/遗传病-交大附属儿医/2019/'+dirname[-2]+'\n')
    ofile.write('scp list klyang@122.112.248.194:/media/bpdata/data/exon_temp2/dir_list/list_'+dirname[-1]+'\n')
    ofile.write('curl \'https://oapi.dingtalk.com/robot/send?access_token=dc5bf30a873950e8c73bb37abf1851ebaf5d83f2792f08684ace2be621a12563\' -H \'Content-Type: application/json\' -d \'{"msgtype": "text","text": {"content": "'+dirname[-1]+' backup done"}}\'\n')
    ofile.write('ssh kly\n\n')
    return ofile

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
    -redo   yes
    Note:
    make sure excel colnames contain 批次   收样日期    样本编号    原始样本ID    项目    性别    年龄    科室    浓度(ng/ul)    体积(ul)   总量(ug)    OD260/OD280 OD260/OD230 质检结果    上机日期    下机日期    Reads（M）   bases(Gb)  Q30
    the rawdata and md5 fileare in the same single dir 
    """
    sys.exit()

kwargs={'-sn':'', '-panel':'', '-filtmode':'', '-rawdir':'', '-md5':'', '-R1':'', '-R2':'', '-outdir':'',  '-analysis_data':'', '-redo':''}

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
dirname = os.getcwd().split('/')[-2:]
if u'分析任务单' in sn:
    sep = u'分析任务单'
else:
    sep = u'产品信息表'

ofile = open('step5_backup.sh', 'w')
backup(ofile, dirname, sn, sep, panel)
ofile.close()
os.system('sh step5_backup.sh')
