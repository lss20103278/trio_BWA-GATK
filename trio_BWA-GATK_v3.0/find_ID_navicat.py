#!/usr/bin/env python
# -*- coding: utf-8 -*-
#########################################################################
# File Name: find_ID_navicat.py
# Author: sslyu
# mail: 1573077420@qq.com
# Created Time:Wed 15 May 2019 10:50:56 AM CST
#########################################################################

import numpy as np
import pandas as pd
import sys
reload(sys)
sys.setdefaultencoding("utf8")
import os
from lib.generate_trio_cmd_prepare import *

def extract_ID_navicat(ID):
    global navicat
    if ID in navicat[u'样本ID'].tolist():
        sample = navicat[navicat[u'样本ID'] == ID]
    elif ID in navicat[u'样本编号'].tolist():
        sample = navicat[navicat[u'样本编号'] == ID]
    else:
        print "%s not in navicat" % ID
        sys.exit()
    if not sample[u'样本间关系'].isnull().any():
        famID = sample[u'样本间关系'].tolist()[0]
        sample = navicat[navicat[u'样本间关系'] == famID]
    pedigree = generate_pedigree(sample)
    for i in sample.index:
        for k in pedigree.keys():
            if u'爸爸' in sample.loc[i,u"姓名"]:
                sample.loc[i,u"姓名"] = k+u'父'
            elif u'妈妈' in sample.loc[i,u"姓名"]:
                sample.loc[i,u"姓名"] = k+u'母'
    return sample
def generate_file(sample):        
    global ID
    sample_header = [u"批次", u"收样日期", u"样本编号", u"原始样本ID", u"姓名", u"性别", u"年龄", u"科室", u"浓度(ng/ul)", u"体积(ul)", u"质量(ug)", u"OD260/OD280", u"OD260/OD230", u"质检结果", u"上机日期", u"下机日期", u"Reads（M）", u"bases(Gb)", u"Q30", u"exon"]
    sample.rename(columns={u'样本ID':u'原始样本ID', u'申请科室':u"科室", u'exp_reads':u"Reads（M）", u'exp_bases(Mb)':u"bases(Gb)", u'exp_Q30':u"Q30", u'exp_浓度(ng/ul)':u"浓度(ng/ul)", u'exp_体积(ul)':u"体积(ul)", u'exp_质量(ug)':u"质量(ug)", u'exp_OD260/OD280':u"OD260/OD280", u'exp_OD260/OD230':u"OD260/OD230", u'exp_QC结论':u"质检结果", u'process_送样日期':u"收样日期", u'process_上机日期':u"上机日期", u'process_下机日期':u"下机日期"}, inplace=True)
    sample['exon'] = ""
    sample[sample_header].to_csv(u"儿童遗传病（交大附属儿医）分析任务单_"+"".join(ID)+".txt", index=None, sep="\t", encoding="utf-8")

navicat = pd.read_csv('/DATA/sslyu/Project/Genetics_children_hospital/navicat_bak.csv', encoding="utf-8", sep="\t")

ID = sys.argv[1].split(',')
if len(ID) == 1:
    sample = extract_ID_navicat(ID[0])
else:
    sample = pd.DataFrame(columns = navicat.columns)
    for sn in ID:
        tmp = extract_ID_navicat(sn)
        sample = sample.append(tmp, ignore_index=True)
print sample
generate_file(sample)        
