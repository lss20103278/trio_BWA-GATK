#!/usr/bin/env python
# -*- coding: utf-8 -*-
#########################################################################
# File Name: generate_ped.py
# Author: sslyu
# mail: 1573077420@qq.com
# Created Time:Fri 28 Sep 2018 05:15:56 PM CST
#########################################################################

import numpy as np
import pandas as pd
import sys
reload(sys)
sys.setdefaultencoding("utf8")
import os
import re
import collections
import shutil

script_path = os.path.split(os.path.abspath(sys.argv[0]))[0]
print script_path

def append_sample_excel(df):
    ID = []
    if u'样本编号' in df.columns and u'原始样本ID' in df.columns:
        for i in range(len(df[u'原始样本ID'])):
            if df[u'原始样本ID'][i] == 'nan':
                ID.append(str(df.iloc[i][u'样本编号']))
            else:
                ID.append(str(df.iloc[i][u'原始样本ID']))
    df['sample'] = ID
    df['sample'] = df['sample'].apply(str)
    return df

def append_sample_txt(df):
    ID = []
    if u'样本编号' in df.columns and u'原始样本ID' in df.columns:
        print u'原始样本ID'
        for i in range(len(df[u'原始样本ID'].isnull())):
            if df[u'原始样本ID'].isnull()[i]:                
                ID.append(df.iloc[i][u'样本编号'])
            else:
                ID.append(df.iloc[i][u'原始样本ID'])
    df['sample'] = ID
    df['sample'] = df['sample'].apply(str)                                                                                        
    return df

def generate_pedigree(df):
    pedigree = {}
    if u'姓名' in df.columns.tolist():
        df.index = df[u'姓名']
    elif 'name' in df.columns.tolist():
        df.index = df['name']
    else:
        print "lacks name column"
        sys.exit()
    for i in df.index:
        pedigree[i] = []
    for k in pedigree:
        for i in df.index:
            if k in i: # the name of the child should be totally included in other family member's names
                pedigree[k].append(i)
    list_keys = pedigree.keys()
    for k in list_keys:
        if len(pedigree[k]) == 1:
            pedigree.pop(k)
    return pedigree

def flatten_list(nested):
    all_item = []
    for sublist in nested:
        for item in sublist:
            all_item.append(item)
    return all_item            
    #if isinstance(nested, list):
    #    for sublist in nested:                                                                                                  
    #        for item in flatten_list(sublist):
    #            yield item
    #else:
    #    yield nested    # flatten_list(pedigree_dict.values()) problem: only return the first list of pedigree_dict.values()

def generate_single(df):
    pedigree_dict = generate_pedigree(df)
    trio = pedigree_dict.keys()
    trio_m = []
    if u'姓名' in df.columns.tolist():
        for i in df[u'姓名']:
            for k in pedigree_dict:
                trio_m.extend(pedigree_dict[k])
        single = [i for i in df[u'姓名'] if i not in trio_m]
    elif 'name' in df.columns.tolist():
        for i in df['name']:
            for k in pedigree_dict:
                trio_m.extend(pedigree_dict[k])
        single = [i for i in df['name'] if i not in trio_m]
    return single

def generate_exon(sample,exon):
    global script_path
    single_cmd = open(sample+'/exon.sh', 'w')
    single_cmd.write('#!/bin/sh\n\n#SBATCH -J '+sample+'\n#SBATCH -p BP10\n#SBATCH -N 1\n#SBATCH -n 4\n#SBATCH --mail-type=FAIL\n#SBATCH --mail-user=1573077420@qq.com\n\n')
    single_cmd.write('sample='+sample+'\n\n')
    if len(exon) > 0 and exon[0] != "":
        for i in exon:
            if os.path.exists(script_path+'/bed/gene/'+i+'.bed'):
                single_cmd.write('samtools depth -b %s/bed/gene/%s.bed /BP12_share/sslyu/bam/%s.sort.mkdup.bam > ../exon/%s_%s.bed.depth\n' % (script_path,i,sample,sample,i))
            else:
                single_cmd.write('python %s/bed/gene/CDS/extract_CDS_region.py %s' % (script_path,i))
                single_cmd.write('python %s/bed/gene/Exon/extract_exon_region.py %s' % (script_path,i))
                single_cmd.write('python %s/bed/gene/merge_exon_cds.py %s' % (script_path,i))
                single_cmd.write('samtools depth -b %s/bed/gene/%s.bed /BP12_share/sslyu/bam/%s.sort.mkdup.bam > ../exon/%s_%s.bed.depth\n' % (script_path,i,sample,sample,i))
                #print "%s/bed/gene/%s.bed doesn't exist" % (script_path,i)
                #sys.exit()
        single_cmd.write('cd ../exon\n')
        for i in exon:
            if os.path.exists(script_path+'/bed/gene/CDS/'+i+'.bed'):
                single_cmd.write('python /DATA/sslyu/trio_BWA-GATK_v3.0/bed/gene/CDS/cds_depth_barplot.py -gene %s -sn %s\n' % (i, sample))
            if os.path.exists(script_path+'/bed/gene/Exon/'+i+'.bed'):
                single_cmd.write('python /DATA/sslyu/trio_BWA-GATK_v3.0/bed/gene/Exon/exon_depth_barplot.py -gene %s -sn %s\n' % (i, sample))
        single_cmd.write('cd ../%s\n' % sample)
    single_cmd.close()

if not os.path.exists('config'):
    #shutil.copy(script_path+'/src/config', '.')
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
    make sure excel/csv colnames contain    样本编号    原始样本ID    姓名    性别   exon   cmd
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
workpath = os.getcwd().split('/')[-2:]
if u'分析任务单' in sn:
    sep = u'分析任务单'
else:                                                                                            
    sep = u'产品信息表'
if sn.endswith('xlsx'):
    excel = pd.read_excel(sn, dtype=str)
elif sn.endswith('csv'):
    excel = pd.read_csv(sn, dtype=str, encoding="utf-8", sep="\t")
elif sn.endswith('txt'):
    excel = pd.read_csv(sn, dtype=str, encoding="utf-8", sep="\t")
else:
    print "the format of the input file is wrong ..."
    sys.exit()

#for i in [u'样本编号', u'原始样本ID', u'姓名', u'性别']:
#    if i not in excel.columns:
#        print i
#        print "the header of the excel/csv is not correct, exiting ..."
#        sys.exit()
if 'sample' not in excel.columns:
    if sn.endswith('xlsx'):
        excel = append_sample_excel(excel)
    else:
        excel = append_sample_txt(excel)
d1 = excel
d1 = d1.replace("nan", "")  # important
d1 = d1.fillna('')  # important

if not os.path.exists('exon'):
    os.mkdir('exon')

pedigree_dict = generate_pedigree(d1) # dictionary of pedigree samples, keys:names of patient samples, values:names of the whole family members including the patient sample
trio = pedigree_dict.keys() 
single = generate_single(d1)
trio_all = [i for i in d1.index if i not in single]
if 'exon' not in d1.columns:
    d1['exon'] = ''

for i in d1['sample']:                            
    if not os.path.exists(str(i)):
        os.mkdir(str(i))
list_exon_sample = open('list_exon_sample', 'w')
if len(single) > 0:
    for i in single:
        sample = d1.loc[i]['sample'] # index:name
        print 'single '+sample
        if os.path.exists(sample+'/exon.sh'):
            os.remove(sample+'/exon.sh')
        try:
            exon = d1.loc[i]['exon'].split(',')
            print exon
        except:
            print "lacks exon column"
            sys.exit()
        if len(exon) > 0 and exon[0] != "":
            generate_exon(sample, exon)
            list_exon_sample.write(sample+'\n')
if len(trio) > 0:
    for k in pedigree_dict:
        for i in pedigree_dict[k]:
            sample = d1.loc[i]['sample']
            print 'trio '+sample
            try:
                exon = d1.loc[i]['exon'].split(',')
                print exon
            except:
                print "lacks exon column"
                sys.exit()
            if len(exon) > 0 and exon[0] != "":
                generate_exon(sample,exon)
                list_exon_sample.write(sample+'\n')
list_exon_sample.close()
