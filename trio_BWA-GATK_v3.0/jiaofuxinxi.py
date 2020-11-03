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

from lib.jiaofuxinxi_prepare import *

# requirements: qcsum,qcvcf,metrics,info of every sample

#def time_convert(time_str):
#    time_obj = datetime.datetime.strptime(str(time_str), '%Y-%m-%d %H:%M:%S')
#    time_converted = time_obj.strftime('%Y-%m-%d')
#    return time_converted
#
#def append_fenxiliucheng(d1):
#    trio = [item for item, count in collections.Counter(d1['familyname'].tolist()).items() if count > 1]
#    single = [item for item, count in collections.Counter(d1['familyname'].tolist()).items() if count == 1]
#    ped_v = []
#    for j in d1['sample']:
#        if j in single:
#            ped_v.append(0)
#        else:
#            ped_v.append(1)
#    d1['pedigree'] = ped_v
#    d1[u'分析流程'] = d1['pedigree'].map(lambda x:'single' if x==0 or x=='nan' else 'trio')
#    return d1
#
#def flatten_list(nested):                                                                                                       
#    all_item = []
#    for sublist in nested:
#        for item in sublist:
#            all_item.append(item)
#    return all_item            
#    #if isinstance(nested, list):
#    #    for sublist in nested:                                                                                                  
#    #        for item in flatten_list(sublist):
#    #            yield item
#    #else:
#    #    yield nested
# 
#def generate_qcsum(d1):
#    qcsum = pd.DataFrame(columns=['sample', '1Xcov', '20Xcov', 'AvgDepth', 'duplication'])
#    for j in d1['sample']:
#        fn = j+'/2_mapping/'+j+'.qcsum'
#        df = pd.read_table(fn, dtype=str)
#        qcsum = pd.concat([qcsum,df], ignore_index=True)
#    qcsum.index = qcsum['sample']
#    qcsum.to_csv('all.qcsum', sep="\t", index=None)
#    return qcsum
#
#def generate_jiaofuxinxi(d1,merge_columns, qcsum,panel,sep,sn):
#    d1 = append_fenxiliucheng(d1)
#    merge_columns.append(u'分析流程')
#    d1_merge = d1[merge_columns]
#    qcsum_merge = qcsum.drop(u'duplication', axis=1)
#    jiaofuxinxi_table = pd.merge(d1_merge, qcsum_merge, on="sample") # notice the type of the value
#    jiaofuxinxi_full = pd.merge(d1, qcsum_merge, on="sample")
#    l2 = [u'策略', 'vcf', u'质检结果']
#    l3 = [panel, 'available', u'合格']
#    for j in range(len(l2)):
#        if l2[j] not in jiaofuxinxi_table.columns:
#            jiaofuxinxi_table.insert(len(jiaofuxinxi_table.columns.tolist()), l2[j], l3[j])
#            jiaofuxinxi_full.insert(len(jiaofuxinxi_full.columns.tolist()), l2[j], l3[j])
#    jiaofuxinxi_table = jiaofuxinxi_table.drop('sample', axis=1)
#    #print jiaofuxinxi_table.columns[16:19]
#    jiaofuxinxi_table_index = []
#    for i in jiaofuxinxi_table.index:
#        if jiaofuxinxi_table[u'原始样本ID'].isnull()[i]:
#            jiaofuxinxi_table_index.append(jiaofuxinxi_table[u'样本编号'][i])
#        else:
#            jiaofuxinxi_table_index.append(jiaofuxinxi_table[u'原始样本ID'][i])
#    jiaofuxinxi_table.index = jiaofuxinxi_table_index
#    if not os.path.exists('Q30'):
#        os.system('python /DATA/sslyu/trio_BWA-GATK_v3.0/src/fqstat_info.py -l list > Q30')
#    #Q30 = pd.read_table('Q30', header=None, index_col=0)
#    #Q30.index = Q30.index.astype(str)
#    Q30 = pd.read_table('Q30', index_col=0)
#    Q30.index = Q30.index.astype(str)
#    jiaofuxinxi_table[jiaofuxinxi_table.columns[16:19]] = Q30[Q30.columns[:3]]
#    #print jiaofuxinxi_table[jiaofuxinxi_table.columns[16:19]]
#    if os.path.exists('peddy/pedtrio/results/peddy.sex_check.csv') and os.path.exists('peddy/pedtrio/results/peddy.ped_check.csv'):
#        sex = pd.read_csv('peddy/pedtrio/results/peddy.sex_check.csv', dtype=str)
#        sex.index = sex['sample_id']
#        sex_check = []
#        for i in jiaofuxinxi_table.index:
#            if sex.loc[i]['error'] == 'True':
#                sex_check.append(u'错误')
#            else:
#                sex_check.append(u'正确')
#        jiaofuxinxi_table[u'性别检查'] = sex_check
#        ped = pd.read_csv('peddy/pedtrio/results/peddy.ped_check.csv', dtype=str)
#        ped_check = {}
#        for i in ped.index:
#            if ped.loc[i]['parent_error'] == 'True':
#                ped_check[ped.loc[i]['sample_a']] = ped.loc[i]['sample_b']
#        ped_check_excel = []
#        pedigree_dict = generate_pedigree(d1)
#        pedigree_dict_values = flatten_list(pedigree_dict.values())
#        for i in jiaofuxinxi_table.index:
#            if jiaofuxinxi_table.loc[i,u'姓名'] not in pedigree_dict_values:
#                #print jiaofuxinxi_table.loc[i,u'姓名'].decode('utf-8')
#                ped_check_excel.append(u'无家系')
#            else:
#                if i in ped_check.keys() or i in ped_check.values():
#                    if u'姐姐' in jiaofuxinxi_table.loc[i,u'姓名'] or u'妹妹' in jiaofuxinxi_table.loc[i,u'姓名'] or u'弟弟' in jiaofuxinxi_table.loc[i,u'姓名'] or u'哥哥' in jiaofuxinxi_table.loc[i,u'姓名']: #处理家系中存在先证者兄弟姐妹的情况，要求：兄弟姐妹的名子为先证者名字+姐姐、妹妹、弟弟或哥哥
#                        ped_check_excel.append(u'有亲子关系')
#                    else:
#                        ped_check_excel.append(u'亲子关系错误')
#                else:
#                    ped_check_excel.append(u'有亲子关系')
#        jiaofuxinxi_table[u'亲缘关系验证'] = ped_check_excel
#    else:
#        jiaofuxinxi_table[u'性别检查'] = 'unknown'
#        jiaofuxinxi_table[u'亲缘关系验证'] = 'unknown'
#
#    ofile = sn.split(sep)[0]+'交付信息表'+u'）'.join(sn.split(sep)[1].split(')'))
#    if 'xlsx' in ofile:
#        ofile = ofile[:-4]+'txt'
#    jiaofuxinxi_table.to_csv(ofile, index=None, sep="\t", encoding="utf-8")
#    if 'txt' in ofile:
#        ofile = ofile[:-3]+'xlsx'
#    wb=pd.ExcelWriter(ofile,engine='openpyxl')
#    jiaofuxinxi_table.to_excel(wb,index=False)
#    wb.save()
#    return jiaofuxinxi_full
#
#def qcsum_check(panel):
#    ofile = open('qcsum_check.log', 'w')
#    d2 = pd.read_table('all.qcsum', dtype=str)
#    d2.index = d2['sample']
#    if panel == 'PANEL' or panel == 'Agilent-PANEL' or panel == 'IDT-PANEL':
#        for i in d2.index:
#            if float(d2.loc[i]['20Xcov']) < 95:
#                ofile.write('the 20X coverage of '+i+' is less than 95%\n')
#    else:
#        for i in d2.index:
#            if float(d2.loc[i]['20Xcov']) < 90:
#                ofile.write('the 20X coverage of '+i+' is less than 90%\n')
#    ofile.close()
#    if os.path.getsize('qcsum_check.log') == 0:
#        os.remove('qcsum_check.log')
#    else:
#        os.system('cat qcsum_check.log')
#
#def Excel(sn, panel):    
#    if u'分析任务单' in sn:
#        sep = u'分析任务单'
#    else:
#        sep = u'产品信息表'
#    
#    if sn.endswith('xlsx'):
#        excel = pd.read_excel(sn, dtype=str)
#    elif sn.endswith('csv'):
#        excel = pd.read_csv(sn, dtype=str, encoding="utf-8", sep="\t")
#    elif sn.endswith('txt'):
#        excel = pd.read_csv(sn, dtype=str, encoding="utf-8", sep="\t")
#    else:
#        print "the format of the input file is wrong ..."
#        sys.exit()
#    for i in [u'批次', u'收样日期', u'样本编号', u'原始样本ID', u'姓名', u'性别', u'年龄', u'科室', u'浓度(ng/ul)', u'体积(ul)', u'质量(ug)', u'OD260/OD280', u'OD260/OD230', u'质检结果', u'上机日期', u'下机日期', u'Reads（M）', 'bases(Gb)', 'Q30']:
#        if i not in excel.columns:
#            print i
#            print "the header of the excel is not correct, exiting ..."
#            sys.exit()
#    excel[u'收样日期'] = excel[u'收样日期'].fillna('')
#    for i in excel.index:
#        if u'00:00:00' in excel.loc[i,u'收样日期']:
#            v = time_convert(excel.loc[i,u'收样日期'])
#            excel.loc[i,u'收样日期'] = v
#    if u'样本编号' not in excel.columns.tolist():
#        print """
#        The header line of the excel is not correct
#        """
#        sys.exit()
#    merge_columns = excel.columns.tolist()
#    merge_columns = [u'批次', u'收样日期', u'样本编号', u'原始样本ID', u'姓名', u'性别', u'年龄', u'科室', u'浓度(ng/ul)', u'体积(ul)', u'质量(ug)', u'OD260/OD280', u'OD260/OD230', u'质检结果', u'上机日期', u'下机日期', u'Reads（M）', 'bases(Gb)', 'Q30', 'sample']
#    if sn.endswith('xlsx'):
#        d1 = append_sample_excel(excel)
#    else:
#        d1 = append_sample_txt(excel)
#    d1 = append_gender(d1)
#    d1 = append_relation(d1)
#    d1 = append_pedigree(d1)
#    d1 = append_phenotype(d1)
#    d1 = append_father_mother(d1)
#    d1 = d1.replace('nan', '')
#
#    if not os.path.exists('all.qcsum'):
#        qcsum = generate_qcsum(d1)
#    else:
#        qcsum = pd.read_csv('all.qcsum', sep="\t")
#        qcsum['sample'] = qcsum['sample'].astype(str) # In case of number-only ID
#        qcsum.index = qcsum['sample']
#    jiaofuxinxi_full = generate_jiaofuxinxi(d1,merge_columns, qcsum,panel,sep,sn)
#    return jiaofuxinxi_full

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

Excel(sn,panel)
qcsum_check(panel)

