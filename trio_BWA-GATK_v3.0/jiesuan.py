#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
version:3.0
by:lss
"""

import numpy as np
import pandas as pd
import sys
reload(sys)
sys.setdefaultencoding("utf8")
import os
import re
import time
import collections
import shutil
import subprocess

from lib.generate_trio_cmd_prepare import *
from lib.jiaofuxinxi_prepare import *
from lib.navicat_prepare import *
from lib.jiesuan_prepare import *

# requirements: qcsum,qcvcf,metrics,info of every sample

#def generate_navicat_jiesuan(navicat,table_map_all,panel,sn,sep):
#    table_map_all.index = table_map_all.index.astype(str)
#    table1 = navicat[[u'姓名', 'exp_reads', 'exp_bases(Mb)', 'exp_Q30', 'ana_coverage1X(%)', 'ana_coverage20X(%)', 'ana_avg_depth(X)', 'ana_indel(genomic)', 'ana_snp(genomic)', 'ana_ts/tv(genomic)', 'ana_snp(exonic)', 'ana_ts/tv(exonic)', u'process_送样日期', u'申请科室']]
#    table1.columns = ['Sample', 'Total Reads', 'Total Bases(Gb)', 'Q30', '1Xcov', '20Xcov', 'Ave Depth', 'Indel_genomic', 'Snp_genomic', 'ts/tv_genomic', 'snp_exonic', 'ts/tv_exonic',u'送检日期', u'送检科室']
#    index = []
#    for i in range(navicat.shape[0]):
#        if navicat[u'样本ID'].isnull()[i]:     
#            index.append(str(navicat[u'样本编号'][i]))
#        elif navicat[u'样本ID'][i] == "nan":
#            index.append(str(navicat[u'样本编号'][i]))
#        else:
#            index.append(str(navicat[u'样本ID'][i]))
#    navicat.index = index
#    table1.index = index
#    table2 = pd.concat([table1,table_map_all], axis=1)
#    if os.path.exists('Q30'):
#        Q30 = pd.read_table('Q30', index_col=0)
#        Q30.index = Q30.index.astype(str)
#        print Q30[Q30.columns[:]]
#    Total_Reads = Q30["readnum"].tolist()
#    table2.insert(table2.columns.tolist().index('Unmapped reads'), 'Valid reads', Total_Reads)
#    print table2
#    mapping_ratio = table2.apply(lambda x:round((int(x['Valid reads'])-int(x['Unmapped reads']))/float(x['Valid reads'])*100,2), axis=1)
#    table2.insert(table2.columns.tolist().index('Unmapped reads')+1, 'Mapping Ratio', [str(i)+'%' for i in mapping_ratio])
#    tmp_list = ['Valid reads','Unmapped reads','Mapping Ratio','UNPAIRED_READ_DUPLICATES','READ_PAIR_DUPLICATES','Percent duplication']        
#    tmp = table2[tmp_list]
#    table2 = table2.drop(tmp_list, axis=1)
#    for i in range(len(tmp_list)):
#        table2.insert(table2.columns.tolist().index('Q30')+1+i, tmp_list[i], tmp[tmp_list[i]])
#    table2.insert(table2.columns.tolist().index(u'送检日期'), u'序号', range(1,table2.shape[0]+1))
#    if panel == 'WES':
#        table2[u'检测项目'] = u'全外显子测序'
#        table2[u'金额/元'] = u'2970'
#    else:
#        table2[u'检测项目'] = panel
#        table2[u'金额/元'] = u'0'
#    ofile = sn.split(sep)[0]+'jiesuan'+u'）'.join(sn.split(sep)[1].split(')'))
#    if 'xlsx' in ofile:
#        ofile = ofile[:-4]+'txt'
#    table2.to_csv(ofile, index=None, sep="\t", encoding="utf-8")
#    #if 'txt' in ofile:
#    #    ofile = ofile[:-3]+'xlsx'
#    #wb = pd.ExcelWriter(ofile,engine='openpyxl')
#    #table2.to_excel(wb,index=False)
#    #wb.save()
#    table3 = pd.concat([navicat,table2], axis=1)
#    table3.index = table3['Sample']
#    total_table = pd.read_csv('/DATA/sslyu/Project/Genetics_children_hospital/navicat_jiesuan.csv', sep="\t", encoding="utf-8")
#    total_table.index = total_table['Sample']
#    idx = np.unique(total_table.index, return_index=True)[1] # important delete multiple rows with the same index
#    total_table = total_table.iloc[idx] # important
#    print table3.index
#    for i in table3.index:
#        if i.decode("utf-8") in total_table.index:
#            print "recurrent sample"
#            total_table = total_table.drop(i.decode("utf-8"))
#    navicat_jiesuan = total_table.append(table3)
#    navicat_jiesuan.to_csv('/DATA/sslyu/Project/Genetics_children_hospital/navicat_jiesuan.csv', index=False, sep="\t", encoding="utf-8")
#
#def Jiesuan(sn, sep, panel):
#    ofile = sn.split(sep)[0]+'navicat'+u'）'.join(sn.split(sep)[1].split(')'))
#    #if 'txt' in ofile:
#    #    ofile = ofile[:-3]+'xlsx'    
#    #try:
#    #    navicat = pd.read_excel(ofile)
#    if 'xlsx' in ofile:
#        ofile = ofile[:-4]+'txt'
#    try:
#        navicat = pd.read_csv(ofile, sep="\t", encoding="utf-8")
#    except:
#        print "%s is lost" % ofile
#    #navicat = Navicat(sn, sep, panel)
#    
#    if not os.path.exists('unmapped_duplicate'):
#        subprocess.call(["sh",m_path+"/src/unmapped_duplicate.sh","list"])
#    
#    #table_map_all = pd.read_table('unmapped_duplicate', header=None, index_col=0)    
#    table_map_all = pd.read_table('unmapped_duplicate', index_col=0)    
#    table_map_all.index = table_map_all.index.astype(str)
#    table_map_all.columns = ['Unmapped reads', 'UNPAIRED_READ_DUPLICATES', 'READ_PAIR_DUPLICATES', 'Percent duplication']
#    
#    generate_navicat_jiesuan(navicat,table_map_all,panel,sn,sep)

m_path = os.path.split(os.path.abspath(sys.argv[0]))[0]
print m_path
if len(sys.argv) > 1:
    cmd = sys.argv[1]

if not os.path.exists('config'):
    #shutil.copy(m_path+'/src/config', '.')
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
    
log = open('log', 'a')

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

Jiesuan(sn, sep, panel)
