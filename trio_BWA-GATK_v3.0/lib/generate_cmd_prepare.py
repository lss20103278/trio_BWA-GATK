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
script_path = "/DATA/sslyu/trio_BWA-GATK_v3.0/"
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

def generate_run2_sentieon(path,str):
    global script_path
    run2_sentieon = open(path+'/run2_sentieon.sh', 'w')
    n = 0
    with open(script_path+'/src/run2_sentieon.sh') as f:
        for l in f:
            n = n+1
            if n>0 and n<= 16:
                run2_sentieon.write(l)
            elif n>17 and n<25:
                run2_sentieon.write(l)
            elif n>= 26:
                run2_sentieon.write(l)
            elif n==25:
                run2_sentieon.write('cmd="'+str+'"\n')
            elif n==17:
                run2_sentieon.write('basedir='+script_path+'/src\n')
            else:
                continue
    run2_sentieon.close()                

def generate_run2_gatk4(path,str):
    global script_path
    run2_gatk4 = open(path+'/run2_gatk4.sh', 'w')
    n = 0
    with open(script_path+'/src/run2_gatk4.sh') as f:
        for l in f:
            n = n+1
            if n>0 and n<= 17:
                run2_gatk4.write(l)
            elif n==18:
                run2_gatk4.write('basedir='+script_path+'/src\n')
            elif n==24:
                run2_gatk4.write('cmd="'+str+'"\n')
            else:
                run2_gatk4.write(l)
    run2_gatk4.close()

def generate_cmd_gatk4_part1(sample,exon):
    global script_path
    single_cmd = open(sample+'/cmd.sh', 'w')
    single_cmd.write('#!/bin/sh\n\n#SBATCH -J '+sample+'\n#SBATCH -p BP10\n#SBATCH -N 1\n#SBATCH -n 4\n#SBATCH --mail-type=FAIL\n#SBATCH --mail-user=1573077420@qq.com\n\n')
    single_cmd.write('sample='+sample+'\n\n')
    single_cmd.write('sh run2_gatk4.sh $sample $SLURM_NPROCS\n')
    single_cmd.write('/home/ana005/anaconda2/bin/iTools Fqtools stat -InFq ../raw/'+sample+'_R1.fastq.gz -InFq ../raw/'+sample+'_R2.fastq.gz -OutStat '+sample+'.info\n')
    #if len(exon) > 0 and exon[0] != "":
    #    for i in exon:
    #        if os.path.exists(script_path+'/bed/gene/'+i+'.bed'):
    #            single_cmd.write('samtools depth -b %s/bed/gene/%s.bed /BP12_share/sslyu/bam/%s.sort.mkdup.bam > ../exon/%s_%s.bed.depth\n' % (script_path,i,sample,sample,i))
    #        else:
    #            print "%s/bed/gene/%s.bed doesn't exist" % (script_path,i)
    #            sys.exit()
    #    single_cmd.write('cd ../exon\n')
    #    for i in exon:
    #        single_cmd.write('python /DATA/sslyu/trio_BWA-GATK_v3.0/src/exon_depth_barplot.py -gene %s -sn %s\n' % (i, sample))
    #    single_cmd.write('cd ../%s\n' % sample)
    single_cmd.write('scp /BP12_share/sslyu/bam/'+sample+'.sort.mkdup.bam* klyang@122.112.248.194:/media/bpdata/data/exon_temp2/\n')                                                                                                                                
    single_cmd.write('if [ -e ../gvcf/%s.g.vcf ]; then ' % sample)
    single_cmd.write('curl \'https://oapi.dingtalk.com/robot/send?access_token=dc5bf30a873950e8c73bb37abf1851ebaf5d83f2792f08684ace2be621a12563\' -H \'Content-Type: application/json\' -d \'{"msgtype": "text","text": {"content": "'+sample+' mapping done"}}\'\n')
    single_cmd.write('else curl \'https://oapi.dingtalk.com/robot/send?access_token=dc5bf30a873950e8c73bb37abf1851ebaf5d83f2792f08684ace2be621a12563\' -H \'Content-Type: application/json\' -d \'{"msgtype": "text","text": {"content": "'+sample+' run2_gatk4.sh error"}}\'\n')
    single_cmd.write('fi\n')
    single_cmd.close()

def generate_cmd_sentieon_part1(sample,exon):
    global script_path
    single_cmd = open(sample+'/cmd.sh', 'w')
    single_cmd.write('#!/bin/sh\n\n#SBATCH -J '+sample+'\n#SBATCH -p BP10\n#SBATCH -N 1\n#SBATCH -n 4\n#SBATCH --mail-type=FAIL\n#SBATCH --mail-user=1573077420@qq.com\n\n')
    single_cmd.write('sample='+sample+'\n\n')
    single_cmd.write('sh run2_sentieon.sh $sample $SLURM_NPROCS\n')
    single_cmd.write('/home/ana005/anaconda2/bin/iTools Fqtools stat -InFq ../raw/'+sample+'_R1.fastq.gz -InFq ../raw/'+sample+'_R2.fastq.gz -OutStat '+sample+'.info\n')
    #if len(exon) > 0 and exon[0] != "":
    #    for i in exon:
    #        if os.path.exists(script_path+'/bed/gene/'+i+'.bed'):
    #            single_cmd.write('samtools depth -b %s/bed/gene/%s.bed /BP12_share/sslyu/bam/%s.sort.mkdup.bam > ../exon/%s_%s.bed.depth\n' % (script_path,i,sample,sample,i))
    #        else:
    #            print "%s/bed/gene/%s.bed doesn't exist" % (script_path,i)
    #            sys.exit()
    #    single_cmd.write('cd ../exon\n')
    #    for i in exon:
    #        single_cmd.write('python /DATA/sslyu/trio_BWA-GATK_v3.0/src/exon_depth_barplot.py -gene %s -sn %s\n' % (i, sample))
    #    single_cmd.write('cd ../%s\n' % sample)                                                                                   
    single_cmd.write('scp /BP12_share/sslyu/bam/'+sample+'.sort.mkdup.bam* klyang@122.112.248.194:/media/bpdata/data/exon_temp2/\n')
    single_cmd.write('if [ -e ../gvcf/%s.g.vcf ]; then ' % sample)
    single_cmd.write('curl \'https://oapi.dingtalk.com/robot/send?access_token=dc5bf30a873950e8c73bb37abf1851ebaf5d83f2792f08684ace2be621a12563\' -H \'Content-Type: application/json\' -d \'{"msgtype": "text","text": {"content": "'+sample+' mapping done"}}\'\n')
    single_cmd.write('else curl \'https://oapi.dingtalk.com/robot/send?access_token=dc5bf30a873950e8c73bb37abf1851ebaf5d83f2792f08684ace2be621a12563\' -H \'Content-Type: application/json\' -d \'{"msgtype": "text","text": {"content": "'+sample+' run2_sentieon.sh error"}}\'\n')
    single_cmd.write('fi\n')
    single_cmd.close()

def generate_cmd_part2(sample, dirname):
    single_cmd = open(sample+'/cmd.sh', 'a')
    single_cmd.write('cp 2_mapping/$sample.depth.sample_gene_summary ../annotation\n')
    single_cmd.write('cp 3_variants/$sample.{ann.hg19_multianno.txt,CADD,maf,link} ../annotation\n')
    single_cmd.write('cd ../annotation\n')
    single_cmd.write('echo $sample >> list\n')
    single_cmd.write('python '+script_path+'/src/score_re.py -sn $sample\n')
    single_cmd.write('python '+script_path+'/src/annotation_filt.py -sn $sample -c_f_m c\n')
    single_cmd.write('mv '+sample+'_children_hospital*xlsx ../'+dirname[-1]+'\n')
    single_cmd.write('cp ../'+sample+'/3_variants/'+sample+'.vcf ../'+dirname[-1]+'/vcf\n')
    single_cmd.write('echo "end `date`" >> ../$sample/finished\n')
    single_cmd.write('samplexlsx=`ls ../%s/%s_children_hospital*xlsx`\n' % (dirname[-1], sample))
    single_cmd.write('if [ -e $samplexlsx ]; then ')
    single_cmd.write('curl \'https://oapi.dingtalk.com/robot/send?access_token=dc5bf30a873950e8c73bb37abf1851ebaf5d83f2792f08684ace2be621a12563\' -H \'Content-Type: application/json\' -d \'{"msgtype": "text","text": {"content": "'+sample+' annotation excel done"}}\'\n')
    single_cmd.write('else curl \'https://oapi.dingtalk.com/robot/send?access_token=dc5bf30a873950e8c73bb37abf1851ebaf5d83f2792f08684ace2be621a12563\' -H \'Content-Type: application/json\' -d \'{"msgtype": "text","text": {"content": "'+sample+' annotation_filt.py error"}}\'; fi\n')
    single_cmd.close()
    ofile = open('annotation/annotation.sh', 'a')
    ofile.write('python '+script_path+'/src/annotation_filt.py -sn '+sample+' -c_f_m c\n')                                        
    ofile.close()

def generate_dbevn(panel):
    dbevn = open('dbevn.sh', 'w')
    if panel != 'none':
        with open(script_path+'/src/dbevn.sh') as f:
            for l in f:
                if l.startswith('#panel'):
                    if panel in l:
                        dbevn.write(l[1:])
                    else:
                        dbevn.write(l)
                else:
                    dbevn.write(l)
    else:
        with open(script_path+'/src/dbevn.sh') as f:
            for l in f:
                dbevn.write(l)
        panel = raw_input('please enter the absolute path of the bed:')                                                           
        dbevn.write('panel=\"'+panel+'\"\n')
    dbevn.close()

#if not os.path.exists('config'):
#    #shutil.copy(script_path+'/src/config', '.')
#    print """
#    Examples: 
#    prepare config
#    format: (delimiter \\t)
#    -sn excelfile 
#    -panel   IDT-PANEL 
#    -filtmode   hard 
#    -rawdir    absolute_path_of_rawdata 
