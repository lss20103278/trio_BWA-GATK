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

from generate_trio_cmd_prepare import *
from jiaofuxinxi_prepare import *

script_path = "/DATA/sslyu/trio_BWA-GATK_v3.0"
# requirements: qcsum,qcvcf,metrics,info of every sample

def generate_qcvcf(d1):
    trio = [item for item, count in collections.Counter(d1['familyname'].tolist()).items() if count > 1]
    single = [item for item, count in collections.Counter(d1['familyname'].tolist()).items() if count == 1]
    qcvcf = pd.DataFrame(columns=['sample', 'indel_genomic', u'snp_genomic', u'ts/tv_genomic',u'snp_exonic', u'ts/tv_exonic'])
    qcvcf_list = []                                                                                                             
    d1.index = d1['sample']
    for j in d1['sample'].tolist():
        if j in single:
            qcvcf_list.append(j)
        elif j in trio:
            qcvcf_list.append(j)
        else:
            continue
    for j in qcvcf_list:
        if j in single:
            fn = j+'/3_variants/'+j+'.qcvcf'
        else:
            fn = 'trio/'+j+'/triotmp/'+j+'.qcvcf'
        df = pd.read_table(fn, dtype=str)
        qcvcf = pd.concat([qcvcf,df], ignore_index=True)
    qcvcf.index = qcvcf['sample']
    qcvcf.to_csv('all.qcvcf', sep="\t", index=None)
    return qcvcf

def generate_navicat(jiaofuxinxi_full, sn, sep, panel):
    if not os.path.exists('all.qcvcf'):
        qcvcf = generate_qcvcf(jiaofuxinxi_full)
    else:
        qcvcf = pd.read_csv('all.qcvcf', sep="\t")
        qcvcf.index = qcvcf['sample']
    if not os.path.exists('all.qcsum'):
        qcsum = generate_qcsum(jiaofuxinxi_full)
    else:
        qcsum = pd.read_csv('all.qcsum', sep="\t")
        qcsum.index = qcsum['sample']
    header = pd.read_excel(''+script_path+'/doc/navicat_header.xlsx')
    #total_table = pd.read_excel(''+script_path+'/bk/navicat_bak.xlsx')
    #total_table = pd.read_excel('/DATA/sslyu/Project/Genetics_children_hospital/navicat_bak.xlsx')                              
    total_table = pd.read_csv('/DATA/sslyu/Project/Genetics_children_hospital/navicat_bak.csv', sep="\t", encoding="utf-8")
    total_table.index = total_table[u'姓名']
    for i in jiaofuxinxi_full[u'姓名']:
        #print total_table.index.tolist().index(i.decode("utf-8"))
        if i.decode("utf-8") in total_table.index.tolist():
            print "recurrent sample"
            total_table = total_table.drop(i.decode("utf-8"))
    colname3 = header.columns
    sample_n = jiaofuxinxi_full.shape[0]
    start = int(total_table['order'].tolist()[-1])+1
    header[u'order'] = range(int(start), int(start)+sample_n)
    jiaofuxinxi_full.rename(columns={u'原始样本ID':u'样本ID', u'科室':u'申请科室'}, inplace=True)
    jiaofuxinxi_full_index = []
    for i in jiaofuxinxi_full.index:
        if jiaofuxinxi_full[u'样本ID'].isnull()[i]:
            jiaofuxinxi_full_index.append(jiaofuxinxi_full[u'样本编号'][i])
        else:
            jiaofuxinxi_full_index.append(jiaofuxinxi_full[u'样本ID'][i])
    jiaofuxinxi_full.index = jiaofuxinxi_full_index
    if not os.path.exists('Q30'):
        os.system('python /DATA/sslyu/trio_BWA-GATK_v3.0/src/fqstat_info.py -l list > Q30')
    #Q30 = pd.read_table('Q30', header=None, index_col=0)
    Q30 = pd.read_table('Q30', index_col=0)
    Q30.index = Q30.index.astype(str)
    print Q30[Q30.columns[:]]
    jiaofuxinxi_full_columns_ix = [0,2,3,-2,4,5,6,-12,7,8,9,10,11,12,13,13,-6,-1]
    jiaofuxinxi_full[jiaofuxinxi_full.columns[16:19]] = Q30[Q30.columns[:3]]
    print jiaofuxinxi_full[jiaofuxinxi_full.columns[16:19]]
    jiaofuxinxi_full.index = range(sample_n)
    header.index = range(sample_n)
    jiaofuxinxi_full_columns_ix = [0,2,3,-2,4,5,6,-12,7,16,17,18,8,9,10,11,12,13,13,-6,-1] # 16,17,18: reads,bases,Q30
    #jiaofuxinxi_full_columns_ix = [0,2,3,-2,4,5,6,-12,7,17,18,19,8,9,10,11,12,13,13,-6,-1] # 17,18,19: reads,bases,Q30
    jiaofuxinxi_full_columns = []
    for i in jiaofuxinxi_full_columns_ix:
        jiaofuxinxi_full_columns.append(jiaofuxinxi_full.columns.tolist()[i])
    print '\t'.join(jiaofuxinxi_full_columns)
    header[[u'批次', u'样本编号', u'样本ID', u'策略', u'姓名', u'性别', u'年龄',u'样本间关系', u'申请科室',u'exp_reads',u'exp_bases(Mb)', u'exp_Q30', u'exp_浓度(ng/ul)', u'exp_体积(ul)',u'exp_质量(ug)', u'exp_OD260/OD280', u'exp_OD260/OD230', u'exp_QC结论',u'ana_qc_comment',u'ana_pepline', u'ana_vcf']] = jiaofuxinxi_full[jiaofuxinxi_full_columns[:21]]
    #header[[u'批次', u'样本编号', u'样本ID', u'策略', u'姓名', u'性别', u'年龄',u'样本间关系', u'申请科室',u'exp_reads',u'exp_bases(Mb)', u'exp_Q30', u'exp_浓度(ng/ul)', u'exp_体积(ul)',u'exp_质量(ug)', u'exp_OD260/OD280', u'exp_OD260/OD230', u'exp_QC结>>论',u'ana_qc_comment',u'ana_pepline', u'ana_vcf']] = jiaofuxinxi_full[jiaofuxinxi_full_columns] ValueError: Wrong number of items passed 2, placement implies 1
    header.index = jiaofuxinxi_full.index # very important
    header[u'家系关系'] = None
    header[u'家系关系'] = jiaofuxinxi_full[u'家系关系']                                                                                       
    header[[u'process_送样日期', u'process_上机日期', u'process_下机日期']] = None
    header[[u'process_送样日期', u'process_上机日期', u'process_下机日期']] = jiaofuxinxi_full[[u'收样日期', u'上机日期',u'下机日期']]
    header.index = header[u'样本ID']
    header[[u'ana_coverage1X(%)',u'ana_coverage20X(%)', u'ana_avg_depth(X)', u'ana_duplication(%)']] = qcsum[[u'1Xcov', u'20Xcov', u'AvgDepth', u'duplication']]
    header[[u'ana_indel(genomic)', u'ana_snp(genomic)', u'ana_ts/tv(genomic)',u'ana_snp(exonic)', u'ana_ts/tv(exonic)']] = qcvcf[[u'indel_genomic', u'snp_genomic', u'ts/tv_genomic',u'snp_exonic', u'ts/tv_exonic']]
    annotation_ver = "3.0"
    header[u'ana_annotation'] = annotation_ver
    header[u'ana_reporter'] = u'吕珊珊'
    #bk_rawdata = u'/anaData/anaData004/children_hos_genetic/rawdata/2018/'+panel+'/'
    bk_rawdata = u'/anaData/anaData004/children_hos_genetic/rawdata/2019/'+panel+'/'
    for i in jiaofuxinxi_full[u'样本编号']:
        os.system('ls '+bk_rawdata+i)
    print bk_rawdata
    header[u'bk_rawdata'] = bk_rawdata
    #bk_reports = u'/DATA/BPshare/opm/遗传病-交大附属儿医/2018/'+'/'.join(os.getcwd().split('/')[-2:])
    bk_reports = u'/DATA/BPshare/opm/遗传病-交大附属儿医/2019/'+'/'.join(os.getcwd().split('/')[-2:])
    header[u'bk_reports'] = bk_reports
    print header
    #result = total_table.append(header)
    header.index = header[u'姓名']
    result = pd.concat([total_table,header])

    ofile = sn.split(sep)[0]+'navicat'+u'）'.join(sn.split(sep)[1].split(')'))
    if 'xlsx' in ofile:
        ofile = ofile[:-4]+'txt'
    header.to_csv(ofile, index=None, sep="\t", encoding="utf-8")
    #if 'txt' in ofile:
    #    ofile = ofile[:-3]+'xlsx'
    #wb=pd.ExcelWriter(ofile, engine='openpyxl')
    #header.to_excel(wb,index=False)                                                                                             
    #wb.save()
    database=raw_input('Insert into database? yes or no: ')
    if database == 'yes':
        ofile1 = script_path+'/doc/navicat.xlsx'
        wb1=pd.ExcelWriter(ofile1, engine='openpyxl')
        result.to_excel(wb1,index=False)
        wb1.save()
        #os.system('cp '+script_path+'/doc/navicat.xlsx '+script_path+'/bk/navicat_bak.xlsx')
        os.system('cp '+script_path+'/doc/navicat.xlsx /DATA/sslyu/Project/Genetics_children_hospital/navicat_bak.xlsx')
        result.to_csv('/DATA/sslyu/Project/Genetics_children_hospital/navicat_bak.csv', index=None, sep="\t", encoding="utf-8")
    return header

def Navicat(sn, sep, panel):
    jiaofuxinxi_full = Excel(sn,panel)
    navicat = generate_navicat(jiaofuxinxi_full, sn, sep, panel)
    return navicat

