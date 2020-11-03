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
from lib.generate_cmd_prepare import *

script_path = os.path.split(os.path.abspath(sys.argv[0]))[0]
print script_path

#def append_sample_excel(df):
#    ID = []
#    if u'样本编号' in df.columns and u'原始样本ID' in df.columns:
#        for i in range(len(df[u'原始样本ID'])):
#            if df[u'原始样本ID'][i] == 'nan':
#                ID.append(str(df.iloc[i][u'样本编号']))
#            else:
#                ID.append(str(df.iloc[i][u'原始样本ID']))
#    df['sample'] = ID
#    df['sample'] = df['sample'].apply(str)
#    return df
#
#def append_sample_txt(df):
#    ID = []
#    if u'样本编号' in df.columns and u'原始样本ID' in df.columns:
#        print u'原始样本ID'
#        for i in range(len(df[u'原始样本ID'].isnull())):
#            if df[u'原始样本ID'].isnull()[i]:                
#                ID.append(df.iloc[i][u'样本编号'])
#            else:
#                ID.append(df.iloc[i][u'原始样本ID'])
#    df['sample'] = ID
#    df['sample'] = df['sample'].apply(str)                                                                                        
#    return df
#
#def generate_pedigree(df):
#    pedigree = {}
#    if u'姓名' in df.columns.tolist():
#        df.index = df[u'姓名']
#    elif 'name' in df.columns.tolist():
#        df.index = df['name']
#    else:
#        print "lacks name column"
#        sys.exit()
#    for i in df.index:
#        pedigree[i] = []
#    for k in pedigree:
#        for i in df.index:
#            if k in i: # the name of the child should be totally included in other family member's names
#                pedigree[k].append(i)
#    list_keys = pedigree.keys()
#    for k in list_keys:
#        if len(pedigree[k]) == 1:
#            pedigree.pop(k)
#    return pedigree
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
#    #    yield nested    # flatten_list(pedigree_dict.values()) problem: only return the first list of pedigree_dict.values()
#
#def generate_single(df):
#    pedigree_dict = generate_pedigree(df)
#    trio = pedigree_dict.keys()
#    trio_m = []
#    if u'姓名' in df.columns.tolist():
#        for i in df[u'姓名']:
#            for k in pedigree_dict:
#                trio_m.extend(pedigree_dict[k])
#        single = [i for i in df[u'姓名'] if i not in trio_m]
#    elif 'name' in df.columns.tolist():
#        for i in df['name']:
#            for k in pedigree_dict:
#                trio_m.extend(pedigree_dict[k])
#        single = [i for i in df['name'] if i not in trio_m]
#    return single
#
#def generate_run2_sentieon(path,str):
#    global script_path
#    run2_sentieon = open(path+'/run2_sentieon.sh', 'w')
#    n = 0
#    with open(script_path+'/src/run2_sentieon.sh') as f:
#        for l in f:
#            n = n+1
#            if n>0 and n<= 16:
#                run2_sentieon.write(l)
#            elif n>17 and n<25:
#                run2_sentieon.write(l)
#            elif n>= 26:
#                run2_sentieon.write(l)
#            elif n==25:
#                run2_sentieon.write('cmd="'+str+'"\n')
#            elif n==17:
#                run2_sentieon.write('basedir='+script_path+'/src\n')
#            else:
#                continue
#    run2_sentieon.close()                
#
#def generate_run2_gatk4(path,str):
#    global script_path
#    run2_gatk4 = open(path+'/run2_gatk4.sh', 'w')
#    n = 0
#    with open(script_path+'/src/run2_gatk4.sh') as f:
#        for l in f:
#            n = n+1
#            if n>0 and n<= 17:
#                run2_gatk4.write(l)
#            elif n==18:
#                run2_gatk4.write('basedir='+script_path+'/src\n')
#            elif n==24:
#                run2_gatk4.write('cmd="'+str+'"\n')
#            else:
#                run2_gatk4.write(l)
#    run2_gatk4.close()
#
#def generate_cmd_gatk4_part1(sample,exon):
#    global script_path
#    single_cmd = open(sample+'/cmd.sh', 'w')
#    single_cmd.write('#!/bin/sh\n\n#SBATCH -J '+sample+'\n#SBATCH -p BP10\n#SBATCH -N 1\n#SBATCH -n 4\n#SBATCH --mail-type=FAIL\n#SBATCH --mail-user=1573077420@qq.com\n\n')
#    single_cmd.write('sample='+sample+'\n\n')
#    single_cmd.write('sh run2_gatk4.sh $sample $SLURM_NPROCS\n')
#    single_cmd.write('/home/ana005/anaconda2/bin/iTools Fqtools stat -InFq ../raw/'+sample+'_R1.fastq.gz -InFq ../raw/'+sample+'_R2.fastq.gz -OutStat '+sample+'.info\n')
#    if len(exon) > 0 and exon[0] != "":
#        for i in exon:
#            if os.path.exists(script_path+'/bed/gene/'+i+'.bed'):
#                single_cmd.write('samtools depth -b %s/bed/gene/%s.bed /BP12_share/sslyu/bam/%s.sort.mkdup.bam > ../exon/%s_%s.bed.depth\n' % (script_path,i,sample,sample,i))
#            else:
#                print "%s/bed/gene/%s.bed doesn't exist" % (script_path,i)
#                sys.exit()
#        single_cmd.write('cd ../exon\n')
#        for i in exon:
#            single_cmd.write('python /DATA/sslyu/trio_BWA-GATK_v3.0/src/exon_depth_barplot.py -gene %s -sn %s\n' % (i, sample))
#        single_cmd.write('cd ../%s\n' % sample)
#    single_cmd.write('scp /BP12_share/sslyu/bam/'+sample+'.sort.mkdup.bam* klyang@122.112.248.194:/media/bpdata/data/exon_temp2/\n')                                                                                                                                
#    single_cmd.write('if [ -e ../gvcf/%s.g.vcf ]; then ' % sample)
#    single_cmd.write('curl \'https://oapi.dingtalk.com/robot/send?access_token=dc5bf30a873950e8c73bb37abf1851ebaf5d83f2792f08684ace2be621a12563\' -H \'Content-Type: application/json\' -d \'{"msgtype": "text","text": {"content": "'+sample+' mapping done"}}\'\n')
#    single_cmd.write('else curl \'https://oapi.dingtalk.com/robot/send?access_token=dc5bf30a873950e8c73bb37abf1851ebaf5d83f2792f08684ace2be621a12563\' -H \'Content-Type: application/json\' -d \'{"msgtype": "text","text": {"content": "'+sample+' run2_gatk4.sh error"}}\'\n')
#    single_cmd.write('fi\n')
#    single_cmd.close()
#
#def generate_cmd_sentieon_part1(sample,exon):
#    global script_path
#    single_cmd = open(sample+'/cmd.sh', 'w')
#    single_cmd.write('#!/bin/sh\n\n#SBATCH -J '+sample+'\n#SBATCH -p BP10\n#SBATCH -N 1\n#SBATCH -n 4\n#SBATCH --mail-type=FAIL\n#SBATCH --mail-user=1573077420@qq.com\n\n')
#    single_cmd.write('sample='+sample+'\n\n')
#    single_cmd.write('sh run2_sentieon.sh $sample $SLURM_NPROCS\n')
#    single_cmd.write('/home/ana005/anaconda2/bin/iTools Fqtools stat -InFq ../raw/'+sample+'_R1.fastq.gz -InFq ../raw/'+sample+'_R2.fastq.gz -OutStat '+sample+'.info\n')
#    if len(exon) > 0 and exon[0] != "":
#        for i in exon:
#            if os.path.exists(script_path+'/bed/gene/'+i+'.bed'):
#                single_cmd.write('samtools depth -b %s/bed/gene/%s.bed /BP12_share/sslyu/bam/%s.sort.mkdup.bam > ../exon/%s_%s.bed.depth\n' % (script_path,i,sample,sample,i))
#            else:
#                print "%s/bed/gene/%s.bed doesn't exist" % (script_path,i)
#                sys.exit()
#        single_cmd.write('cd ../exon\n')
#        for i in exon:
#            single_cmd.write('python /DATA/sslyu/trio_BWA-GATK_v3.0/src/exon_depth_barplot.py -gene %s -sn %s\n' % (i, sample))
#        single_cmd.write('cd ../%s\n' % sample)                                                                                   
#    single_cmd.write('scp /BP12_share/sslyu/bam/'+sample+'.sort.mkdup.bam* klyang@122.112.248.194:/media/bpdata/data/exon_temp2/\n')
#    single_cmd.write('if [ -e ../gvcf/%s.g.vcf ]; then ' % sample)
#    single_cmd.write('curl \'https://oapi.dingtalk.com/robot/send?access_token=dc5bf30a873950e8c73bb37abf1851ebaf5d83f2792f08684ace2be621a12563\' -H \'Content-Type: application/json\' -d \'{"msgtype": "text","text": {"content": "'+sample+' mapping done"}}\'\n')
#    single_cmd.write('else curl \'https://oapi.dingtalk.com/robot/send?access_token=dc5bf30a873950e8c73bb37abf1851ebaf5d83f2792f08684ace2be621a12563\' -H \'Content-Type: application/json\' -d \'{"msgtype": "text","text": {"content": "'+sample+' run2_gatk4.sh error"}}\'\n')
#    single_cmd.write('fi\n')
#    single_cmd.close()
#
#def generate_cmd_part2(sample, dirname):
#    single_cmd = open(sample+'/cmd.sh', 'a')
#    single_cmd.write('cp 2_mapping/$sample.depth.sample_gene_summary ../annotation\n')
#    single_cmd.write('cp 3_variants/$sample.{ann.hg19_multianno.txt,CADD,maf,link} ../annotation\n')
#    single_cmd.write('cd ../annotation\n')
#    single_cmd.write('echo $sample >> list\n')
#    single_cmd.write('python '+script_path+'/src/score_re.py -sn $sample\n')
#    single_cmd.write('python '+script_path+'/src/annotation_filt.py -sn $sample -c_f_m c\n')
#    single_cmd.write('mv '+sample+'_children_hospital*xlsx ../'+dirname[-1]+'\n')
#    single_cmd.write('cp ../'+sample+'/3_variants/'+sample+'.vcf ../'+dirname[-1]+'/vcf\n')
#    single_cmd.write('echo "end `date`" >> ../$sample/finished\n')
#    single_cmd.write('samplexlsx=`ls ../%s/%s_children_hospital*xlsx`\n' % (dirname[-1], sample))
#    single_cmd.write('if [ -e $samplexlsx ]; then ')
#    single_cmd.write('curl \'https://oapi.dingtalk.com/robot/send?access_token=dc5bf30a873950e8c73bb37abf1851ebaf5d83f2792f08684ace2be621a12563\' -H \'Content-Type: application/json\' -d \'{"msgtype": "text","text": {"content": "'+sample+' annotation excel done"}}\'\n')
#    single_cmd.write('else curl \'https://oapi.dingtalk.com/robot/send?access_token=dc5bf30a873950e8c73bb37abf1851ebaf5d83f2792f08684ace2be621a12563\' -H \'Content-Type: application/json\' -d \'{"msgtype": "text","text": {"content": "'+sample+' annotation_filt.py error"}}\'; fi\n')
#    single_cmd.close()
#    ofile = open('annotation/annotation.sh', 'a')
#    ofile.write('python '+script_path+'/src/annotation_filt.py -sn '+sample+' -c_f_m c\n')                                        
#    ofile.close()
#
#def generate_dbevn(panel):
#    dbevn = open('dbevn.sh', 'w')
#    if panel != 'none':
#        with open(script_path+'/src/dbevn.sh') as f:
#            for l in f:
#                if l.startswith('#panel'):
#                    if panel in l:
#                        dbevn.write(l[1:])
#                    else:
#                        dbevn.write(l)
#                else:
#                    dbevn.write(l)
#    else:
#        with open(script_path+'/src/dbevn.sh') as f:
#            for l in f:
#                dbevn.write(l)
#        panel = raw_input('please enter the absolute path of the bed:')                                                           
#        dbevn.write('panel=\"'+panel+'\"\n')
#    dbevn.close()

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
    -software   sentieon/gatk4
    Note:
    make sure excel/csv colnames contain    样本编号    原始样本ID    姓名    性别   exon   cmd
    """
    sys.exit()

kwargs={'-sn':'', '-panel':'', '-filtmode':'', '-rawdir':'', '-md5':'', '-R1':'', '-R2':'', '-outdir':'', '-analysis_data':'', '-redo':'', '-software':''}
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
software = kwargs.get('-software')
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

d1['sample'].to_csv('list', header=None, index=None)
report_dirname = os.getcwd().split('/')[-1]
for i in ['gvcf', 'annotation', report_dirname, 'exon']:
    if not os.path.exists(i):
        os.mkdir(i)
if not os.path.exists(report_dirname+'/vcf'):
    os.makedirs(report_dirname+'/vcf')
generate_dbevn(panel)   # important

pedigree_dict = generate_pedigree(d1) # dictionary of pedigree samples, keys:names of patient samples, values:names of the whole family members including the patient sample
trio = pedigree_dict.keys() 
single = generate_single(d1)
trio_all = [i for i in d1.index if i not in single]
if 'exon' not in d1.columns:
    d1['exon'] = ''
if 'cmd' not in d1.columns:
    cmd = []
    for i in d1.index:
        if i in single:
            #cmd.append("clean mapping sorting precalling gvcf index_single gtgvcf hardfilt left-normalize depth qcsum qcvcf annotation")  # gatk4 cmd
            cmd.append("clean mapping precalling gvcf index_single gtgvcf hardfilt left-normalize depth qcsum qcvcf annotation")   # sentieon cmd
        elif i in trio_all:
            #cmd.append("clean mapping sorting precalling gvcf index_trio depth qcsum")  # gatk4 cmd
            cmd.append("clean mapping sorting precalling gvcf index_trio depth qcsum")  # sentieon cmd
    d1['cmd'] = cmd
print d1['cmd']

for i in d1['sample']:                            
    if not os.path.exists(str(i)):
        os.mkdir(str(i))
if len(single) > 0:
    list_single = open('list_single', 'w')
    for i in single:
        sample = d1.loc[i]['sample'] # index:name
        print 'single '+sample
        list_single.write(sample+'\n')
        if os.path.exists(sample+'/cmd.sh'):
            os.remove(sample+'/cmd.sh')
        try:
            exon = d1.loc[i]['exon'].split(',')
            print exon
        except:
            print "lacks exon column"
            sys.exit()
        try:
            cmd_single = d1.loc[i]['cmd']
        except:
            print "lacks cmd columns"
            sys.exit()
        if software == 'sentieon':
            generate_cmd_sentieon_part1(sample, exon)
            generate_run2_sentieon(sample,cmd_single)
        else:
            generate_cmd_gatk4_part1(sample,exon)
            generate_run2_gatk4(sample,cmd_single)
        generate_cmd_part2(sample, workpath)
    list_single.close()
if len(trio) > 0:
    list_trio = open('list_trio', 'w')
    for k in pedigree_dict:
        for i in pedigree_dict[k]:
            sample = d1.loc[i]['sample']
            print 'trio '+sample
            list_trio.write(sample+'\n')
            if os.path.exists(sample+'/cmd.sh'):
                os.remove(sample+'/cmd.sh')
            try:
                exon = d1.loc[i]['exon'].split(',')
                print exon
            except:
                print "lacks exon column"
                sys.exit()
            try:
                cmd_trio = d1.loc[i]['cmd']
            except:
                print "lacks cmd columns"
                sys.exit()
            if software == 'sentieon':
                generate_cmd_sentieon_part1(sample, exon)
                generate_run2_sentieon(sample,cmd_trio)
            else:
                generate_cmd_gatk4_part1(sample,exon)
                generate_run2_gatk4(sample,cmd_trio)
    list_trio.close()

if os.path.exists("annotation/annotation.sh"):
    os.system('sort annotation/annotation.sh |uniq > annotation/tmp; mv annotation/tmp annotation/annotation.sh')

#os.system('for i in `cat list`; do cd $i; sbatch cmd.sh; cd ..; done')    
