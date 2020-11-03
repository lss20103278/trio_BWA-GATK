#!/usr/bin/env python
# -*- coding: utf-8 -*-
#########################################################################
# File Name: /DATA/sslyu/trio_BWA-GATK_v3.0/generate_trio_cmd.py
# Author: sslyu
# mail: 1573077420@qq.com
# Created Time:Wed 08 May 2019 06:35:11 PM CST
#########################################################################

import numpy as np
import pandas as pd
import os
import sys
reload(sys)
sys.setdefaultencoding("utf8")
import shutil
from lib.generate_trio_cmd_prepare import *

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
#def append_gender(df):
#    if u'性别' in df.columns:
#        df['gender'] = df[u'性别'].apply(lambda x:'1' if '男' in x else '2' if '女' in x else 'unknown')
#    else:
#        df['gender'] = 'unknown'
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
#def append_relation(df):
#    if u'姓名' in df.columns.tolist():
#        df.index = df[u'姓名']                                                                  
#    elif 'name' in df.columns.tolist():
#        df.index = df['name']
#    else:
#        print "lacks name column"
#        sys.exit()
#    if 'relationship' not in df.columns:
#        df['relationship'] = None
#    pedigree_dict = generate_pedigree(df)
#    trio = pedigree_dict.keys()
#    for i in df.index:
#        for j in trio:
#            if j in i:
#                if j == i:
#                    df.loc[i,'relationship'] = u'子'
#                elif u'父' in i:
#                    df.loc[i,'relationship'] = u'父'  #there is no '父' in i
#                elif u'母' in i:
#                    df.loc[i,'relationship'] = u'母'  #there is no '母' in i
#                elif u'爷' in i:
#                    df.loc[i,'relationship'] = u'爷'
#                elif u'奶' in i:
#                    df.loc[i,'relationship'] = u'奶'
#                elif u'外公' in i:
#                    df.loc[i,'relationship'] = u'外公'
#                elif u'外婆' in i:
#                    df.loc[i,'relationship'] = u'外婆'
#                elif u'姐' in i:
#                    df.loc[i,'relationship'] = u'姐'
#                elif u'妹' in i:
#                    df.loc[i,'relationship'] = u'妹'
#                elif u'哥' in i:
#                    df.loc[i,'relationship'] = u'哥'
#                elif u'弟' in i:
#                    df.loc[i,'relationship'] = u'弟'
#                else:
#                    df.loc[i,'relationship'] = 'other'
#    df[u'家系关系'] = None
#    df[u'家系关系'] = df['relationship']
#    return df
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
#def append_phenotype(df):
#    pedigree_dict = generate_pedigree(df)
#    trio = pedigree_dict.keys()
#    if u'姓名' in df.columns.tolist():
#        name = u'姓名'
#    elif 'name' in df.columns.tolist():
#        name = 'name'
#    else:
#        print "lacks name column"
#        sys.exit()
#    try:
#        single = [i for i in df['name'] if i not in flatten_list(pedigree_dict.values())] # dict.values()
#    except:
#        single = []
#        #pass
#    df['phenotype1'] = '0'
#    df['phenotype2'] = None
#    if len(trio) > 0 and len(single) > 0:
#        for i in df.index: 
#            if i in trio:
#                df.loc[i,'phenotype2'] = '2'
#            elif i in single:
#                df.loc[i,'phenotype2'] = '2'
#            else:
#                df.loc[i,'phenotype2'] = '1'
#    elif len(trio) > 0:
#        for i in df.index:
#            if i in trio:
#                df.loc[i,'phenotype2'] = '2'
#            else:
#                df.loc[i,'phenotype2'] = '1'
#    else:
#        for i in df.index:
#            df.loc[i,'phenotype2'] = '2'                                                                                        
#    return df
#
#def generate_single(df):
#    pedigree_dict = generate_pedigree(df)
#    trio = pedigree_dict.keys()
#    trio_m = []
#    if u'姓名' in df.columns.tolist():
#        name = u'姓名'
#    elif 'name' in df.columns.tolist():
#        name = 'name'
#    else:
#        print "lacks name column"
#        sys.exit()
#    for i in df[name]:
#        for k in pedigree_dict:
#            trio_m.extend(pedigree_dict[k])
#    single = [i for i in df[name] if i not in trio_m]
#    return single
#
#def append_father_mother(df):
#    df['father'] = '0'
#    df['mother'] = '0'
#    pedigree_dict = generate_pedigree(df)
#    trio = pedigree_dict.keys()
#    if u'姓名' in df.columns.tolist():
#        name = u'姓名'
#    elif 'name' in df.columns.tolist():
#        name = 'name'
#    else:
#        print "lacks name column"
#        sys.exit()
#    try:
#        single = [i for i in df[name] if i not in pedigree_dict.values()[0]]
#    except:                                                                                      
#        single = df[name].tolist()
#    df.index = df[name]
#    for i in trio: # i 患病孩子
#        relation = ''
#        for j in pedigree_dict[i]: 
#            relation = relation+df.loc[j]['relationship']
#        for j in df.index:
#            if j in pedigree_dict[i]: 
#                if u'父' in relation and u'母' in relation: #父母都有
#                    if j == i: #j为患病孩子
#                        for k in pedigree_dict[i]:
#                            if u'父' in k:
#                                father_name = k 
#                            if u'母' in k:
#                                mother_name = k
#                        df.loc[j,'father'] = df.loc[father_name]['sample']
#                        df.loc[j,'mother'] = df.loc[mother_name]['sample']
#                    else: # j为患病孩子家系中其他人
#                        if len(pedigree_dict[i]) == 3: #家系结构为子父母，父母的父母为0
#                            df.loc[j,'father'] = '0' # j为父或母
#                        else: #家系结构为子父母和其他家人
#                            if u'姐' in j or u'妹' in j or u'哥' in j or u'弟' in j: # 家系结构为 孩子，孩子的兄弟姐妹，父母，j为兄>弟姐妹                                                                                                                          
#                                for k in pedigree_dict[i]:
#                                    if u'父' in k:
#                                        father_name = k
#                                    if u'母' in k:
#                                        mother_name = k
#                                df.loc[j,'father'] = df.loc[father_name]['sample']
#                                df.loc[j,'mother'] = df.loc[mother_name]['sample']
#                            elif u'爷' not in relation and u'奶' not in relation and u'外婆' not in relation and u'外公' not in relation: 
#                                if u'父' in j:
#                                    df.loc[j,'father'] = '0'
#                                    df.loc[j,'mother'] = '0'
#                                if u'母' in j:
#                                    df.loc[j,'father'] = '0'
#                                    df.loc[j,'mother'] = '0'
#                            else:
#                                father = raw_input('father ID of '+j+'(if not exist, please enter 0): ')
#                                df.loc[j,'father'] = father
#                                mother = raw_input('mother ID of '+j+'(if not exist, please enter 0):')
#                                df.loc[j,'mother'] = mother
#                else: # 父母一方没有
#                    if j == i: # j为患病孩子
#                        if u'父' in relation: # 只有父亲
#                            for k in pedigree_dict[i]:
#                                if u'父' in k:
#                                    father_name = k
#                            df.loc[j,'father'] = df.loc[father_name]['sample']
#                        if u'母' in relation: # 只有母亲
#                            for k in pedigree_dict[i]:
#                                if u'母' in k:
#                                    mother_name = k
#                            df.loc[j,'mother'] = df.loc[mother_name]['sample']
#                    elif u'父' in j: # j为父亲
#                        if u'爷' not in relation and u'奶' not in relation and u'外婆' not in relation and u'外公' not in relation: 
#                            df.loc[j,'father'] = '0'
#                            df.loc[j,'mother'] = '0'
#                    elif u'母' in j: # j为母亲
#                        if u'爷' not in relation and u'奶' not in relation and u'外婆' not in relation and u'外公' not in relation: 
#                            df.loc[j,'father'] = '0'
#                            df.loc[j,'mother'] = '0'
#                    else: # j为除父母外的其他家人
#                        father = raw_input('father ID of '+j+'(if not exist, please enter 0):')
#                        df.loc[j,'father'] = father
#                        mother = raw_input('mother ID of '+j+'(if not exist, please enter 0):')
#                        df.loc[j,'mother'] = mother
#    return df                        
#
#def generate_trio_cmd_gatk4(k, filtmode, pedigree_dict, d1, dirname):
#    sample = d1.loc[k]['sample']
#    path = 'trio/'+sample
#    if not os.path.exists(path):
#        os.makedirs(path)
#    trio_cmd = open(path+'/trio_cmd_gatk4.sh', 'w')
#    #trio_cmd.write('#!/bin/sh\n\n#SBATCH -J '+sample+'\n#SBATCH -p BP10\n#SBATCH -N 1\n#SBATCH -n 4\n#SBATCH --mail-type=END,FAIL\n#SBATCH --mail-user=1573077420@qq.com\n\nsample='+sample+'\n\n')
#    trio_cmd.write('#!/bin/sh\n\n#SBATCH -J '+sample+'\n#SBATCH -p BP10\n#SBATCH -N 1\n#SBATCH -n 4\n#SBATCH --mail-type=FAIL\n#SBATCH --mail-user=1573077420@qq.com\n\nsample='+sample+'\n\n')
#    trio_cmd.write('filtmode=\''+filtmode+'\'\npath=`pwd`\nvar=`echo $path | awk -F \'/\' \'{print $NF}\'`\n')
#    trio_cmd.write('if [ $var = \'peddy\' ]\nthen\nsh '+script_path+'/src/trio_gatk4.sh $sample $SLURM_NPROCS\n')
#    trio_cmd.write('if [ -e pedtrio/results/peddy.ped_check.csv ]; then ')
#    trio_cmd.write('curl \'https://oapi.dingtalk.com/robot/send?access_token=dc5bf30a873950e8c73bb37abf1851ebaf5d83f2792f08684ace2be621a12563\' -H \'Content-Type: application/json\' -d \'{"msgtype": "text","text": {"content": "peddy done"}}\'\n')
#    trio_cmd.write('else curl \'https://oapi.dingtalk.com/robot/send?access_token=dc5bf30a873950e8c73bb37abf1851ebaf5d83f2792f08684ace2be621a12563\' -H \'Content-Type: application/json\' -d \'{"msgtype": "text","text": {"content": "peddy error"}}\'; fi\n')
#    trio_cmd.write('else\nsh '+script_path+'/src/trio_gatk4.sh $sample $SLURM_NPROCS $filtmode\n')
#    trio_cmd.write('[ -e ../../$sample/2_mapping/$sample.depth.sample_gene_summary ] && cp ../../$sample/2_mapping/$sample.depth.sample_gene_summary ../../annotation\n')
#    trio_cmd.write('cp segtrio/$sample.{ann.hg19_multianno.txt,CADD,maf,link} ../../annotation\n')
#    trio_cmd.write('cd ../../annotation\n')
#    trio_cmd.write('echo $sample >> list\n')
#    trio_cmd.write('python '+script_path+'/src/score_re.py -sn $sample\n')
#    relation = []
#    for j in pedigree_dict[k]:
#        relation.append(d1.loc[j]['relationship'])
#    a = ' '.join(relation)
#    print a
#    annotationfile = open('annotation/annotation.sh', 'a')
#    if u'父' in a and u'母' in a:
#        trio_cmd.write('python '+script_path+'/src/annotation_filt.py -sn $sample -c_f_m c_f_m\n')
#        annotationfile.write('python '+script_path+'/src/annotation_filt.py -sn '+sample+' -c_f_m c_f_m\n')
#    else:
#        if u'父' in a:
#            trio_cmd.write('python '+script_path+'/src/annotation_filt.py -sn $sample -c_f_m c_f\n')
#            annotationfile.write('python '+script_path+'/src/annotation_filt.py -sn '+sample+' -c_f_m c_f\n')
#        else:
#            trio_cmd.write('python '+script_path+'/src/annotation_filt.py -sn $sample -c_f_m c_m\n')
#            annotationfile.write('python '+script_path+'/src/annotation_filt.py -sn '+sample+' -c_f_m c_m\n')
#    trio_cmd.write('mv '+sample+'_children_hospital*xlsx ../'+dirname[-1]+'\n')
#    trio_cmd.write('cp ../trio/'+sample+'/sep/*.vcf ../'+dirname[-1]+'/vcf\n')
#    trio_cmd.write('echo "end `date`" >> ../$sample/finished\n')
#    trio_cmd.write('samplexlsx=`ls ../%s/%s_children_hospital*xlsx`\n' % (dirname[-1], sample))
#    trio_cmd.write('if [ -e $samplexlsx ]; then ')
#    trio_cmd.write('curl \'https://oapi.dingtalk.com/robot/send?access_token=dc5bf30a873950e8c73bb37abf1851ebaf5d83f2792f08684ace2be621a12563\' -H \'Content-Type: application/json\' -d \'{"msgtype": "text","text": {"content": "'+sample+' annotation excel done"}}\'\n')
#    trio_cmd.write('else curl \'https://oapi.dingtalk.com/robot/send?access_token=dc5bf30a873950e8c73bb37abf1851ebaf5d83f2792f08684ace2be621a12563\' -H \'Content-Type: application/json\' -d \'{"msgtype": "text","text": {"content": "'+sample+' trio_cmd_gatk4.sh"}}\'; fi\n')
#    trio_cmd.write('fi\n')
#    trio_cmd.close()
#    annotationfile.close()
#
#def generate_trio_cmd_sentieon(k, filtmode, pedigree_dict, d1, dirname):              
#    sample = d1.loc[k]['sample']
#    path = 'trio/'+sample
#    if not os.path.exists(path):
#        os.makedirs(path)
#    trio_cmd = open(path+'/trio_cmd_sentieon.sh', 'w')
#    #trio_cmd.write('#!/bin/sh\n\n#SBATCH -J '+sample+'\n#SBATCH -p BP10\n#SBATCH -N 1\n#SBATCH -n 4\n#SBATCH --mail-type=END,FAIL\n#SBATCH --mail-user=1573077420@qq.com\n\nsample='+sample+'\n\n')
#    trio_cmd.write('#!/bin/sh\n\n#SBATCH -J '+sample+'\n#SBATCH -p BP10\n#SBATCH -N 1\n#SBATCH -n 4\n#SBATCH --mail-type=FAIL\n#SBATCH --mail-user=1573077420@qq.com\n\nsample='+sample+'\n\n')
#    trio_cmd.write('filtmode=\''+filtmode+'\'\npath=`pwd`\nvar=`echo $path | awk -F \'/\' \'{print $NF}\'`\n')
#    trio_cmd.write('if [ $var = \'peddy\' ]\nthen\nsh '+script_path+'/src/trio_sentieon.sh $sample $SLURM_NPROCS\n')
#    trio_cmd.write('if [ -e pedtrio/results/peddy.ped_check.csv ]; then ')
#    trio_cmd.write('curl \'https://oapi.dingtalk.com/robot/send?access_token=dc5bf30a873950e8c73bb37abf1851ebaf5d83f2792f08684ace2be621a12563\' -H \'Content-Type: application/json\' -d \'{"msgtype": "text","text": {"content": "peddy done"}}\'\n')
#    trio_cmd.write('else curl \'https://oapi.dingtalk.com/robot/send?access_token=dc5bf30a873950e8c73bb37abf1851ebaf5d83f2792f08684ace2be621a12563\' -H \'Content-Type: application/json\' -d \'{"msgtype": "text","text": {"content": "peddy error"}}\'; fi\n')
#    trio_cmd.write('else\nsh '+script_path+'/src/trio_sentieon.sh $sample $SLURM_NPROCS $filtmode\n')
#    trio_cmd.write('[ -e ../../$sample/2_mapping/$sample.depth.sample_gene_summary ] && cp ../../$sample/2_mapping/$sample.depth.sample_gene_summary ../../annotation\n')
#    trio_cmd.write('cp segtrio/$sample.{ann.hg19_multianno.txt,CADD,maf,link} ../../annotation\n')
#    trio_cmd.write('cd ../../annotation\n')
#    trio_cmd.write('echo $sample >> list\n')
#    trio_cmd.write('python '+script_path+'/src/score_re.py -sn $sample\n')
#    relation = []
#    for j in pedigree_dict[k]:
#        relation.append(d1.loc[j]['relationship'])
#    a = ' '.join(relation)
#    print a
#    annotationfile = open('annotation/annotation.sh', 'a')
#    if u'父' in a and u'母' in a:                                                                                               
#        trio_cmd.write('python '+script_path+'/src/annotation_filt.py -sn $sample -c_f_m c_f_m\n')
#        annotationfile.write('python '+script_path+'/src/annotation_filt.py -sn '+sample+' -c_f_m c_f_m\n')
#    else:
#        if u'父' in a:
#            trio_cmd.write('python '+script_path+'/src/annotation_filt.py -sn $sample -c_f_m c_f\n')
#            annotationfile.write('python '+script_path+'/src/annotation_filt.py -sn '+sample+' -c_f_m c_f\n')
#        else:
#            trio_cmd.write('python '+script_path+'/src/annotation_filt.py -sn $sample -c_f_m c_m\n')
#            annotationfile.write('python '+script_path+'/src/annotation_filt.py -sn '+sample+' -c_f_m c_m\n')
#    trio_cmd.write('mv '+sample+'_children_hospital*xlsx ../'+dirname[-1]+'\n')
#    trio_cmd.write('cp ../trio/'+sample+'/sep/*.vcf ../'+dirname[-1]+'/vcf\n')
#    trio_cmd.write('echo "end `date`" >> ../$sample/finished\n')
#    trio_cmd.write('samplexlsx=`ls ../%s/%s_children_hospital*xlsx`\n' % (dirname[-1], sample))
#    trio_cmd.write('if [ -e $samplexlsx ]; then ')
#    trio_cmd.write('curl \'https://oapi.dingtalk.com/robot/send?access_token=dc5bf30a873950e8c73bb37abf1851ebaf5d83f2792f08684ace2be621a12563\' -H \'Content-Type: application/json\' -d \'{"msgtype": "text","text": {"content": "'+sample+' annotation excel done"}}\'\n')
#    trio_cmd.write('else curl \'https://oapi.dingtalk.com/robot/send?access_token=dc5bf30a873950e8c73bb37abf1851ebaf5d83f2792f08684ace2be621a12563\' -H \'Content-Type: application/json\' -d \'{"msgtype": "text","text": {"content": "'+sample+' trio_cmd_sentieon.sh"}}\'; fi\n')
#    trio_cmd.write('fi\n')
#    trio_cmd.close()
#    annotationfile.close()
#
#def generate_single_ped(k,d):
#    d.index = d['sample']
#    ped_mendel = d.loc[k][['sample', 'father', 'mother', 'gender', 'phenotype2']]
#    ped_mendel = ped_mendel.apply(str)
#    fname = d.loc[k]['sample']+'/'+d.loc[k]['sample']+'.mendel.ped'
#    fname = 'ped/'+d.loc[k]['sample']+'.mendel.ped'
#    ped = open(fname, 'w')
#    ped.write(str(d.loc[k]['sample'])+'\t'+'\t'.join(ped_mendel)+'\n')
#    ped.close()
#    if u'姓名' in d.columns.tolist():
#        name = u'姓名'
#    elif 'name' in d.columns.tolist():
#        name = 'name'
#    else:
#        print "lacks name column"
#        sys.exit()
#    d.index = d[name]
#
#def generate_ped(k,d):
#    pedigree_dict = generate_pedigree(d)
#    ped_f = d.loc[pedigree_dict[k]][['sample', 'father', 'mother', 'gender', 'phenotype1']]
#    ped_f.index = ped_f['sample']
#    ped_mendel = d.loc[pedigree_dict[k]][['sample', 'father', 'mother', 'gender', 'phenotype2']]
#    ped_mendel.index = ped_mendel['sample']
#    f1name = 'ped/'+d.loc[k]['sample']+'.ped'
#    f2name = 'ped/'+d.loc[k]['sample']+'.mendel.ped'
#    ped1 = open(f1name, 'w')
#    ped2 = open(f2name, 'w')
#    relation = []
#    for i in pedigree_dict[k]:
#        relation.append(d.loc[i]['relationship'])
#    a = ' '.join(relation)
#    key1 = d.loc[k]['sample']
#    ped1.write(key1+'\t'+'\t'.join(ped_f.loc[key1])+'\n')
#    ped2.write(key1+'\t'+'\t'.join(ped_mendel.loc[key1])+'\n')
#    if '父' in a and '母' in a:
#        for j in ['父', '母']:
#            for i in pedigree_dict[k]:
#                if i != k and j in i:
#                    key2 = d.loc[i]['sample']
#                    ped1.write(key1+'\t'+'\t'.join(ped_f.loc[key2])+'\n')
#                    ped2.write(key1+'\t'+'\t'.join(ped_mendel.loc[key2])+'\n')
#        n_c_f_m = [i for i in pedigree_dict[k] if '父' not in i and '母' not in i and i != k]
#        for i in n_c_f_m:
#            key2 = d.loc[i]['sample']
#            ped1.write(key1+'\t'+'\t'.join(ped_f.loc[key2])+'\n')                                                               
#            ped2.write(key1+'\t'+'\t'.join(ped_mendel.loc[key2])+'\n')
#    elif '父' in a or '母' in a:
#        parent = ''
#        if '父' in a:
#            parent = '父'
#        else:
#            parent = '母'
#        for i in pedigree_dict[k]:
#            if parent in i and i != k:
#                key2 = d.loc[i]['sample']
#                ped1.write(key1+'\t'+'\t'.join(ped_f.loc[key2])+'\n')
#                ped2.write(key1+'\t'+'\t'.join(ped_mendel.loc[key2])+'\n')
#        n_c_f_m = [i for i in pedigree_dict[k] if parent not in i and i != k]
#        for i in n_c_f_m:
#            key2 = d.loc[i]['sample']
#            ped1.write(key1+'\t'+'\t'.join(ped_f.loc[key2])+'\n')
#            ped2.write(key1+'\t'+'\t'.join(ped_mendel.loc[key2])+'\n')
#    else:
#        n_c = [i for i in pedigree_dict[k] if i != k]
#        for i in n_c:
#            key2 = d.loc[i]['sample']
#            ped1.write(key1+'\t'+'\t'.join(ped_f.loc[key2])+'\n')
#            ped2.write(key1+'\t'+'\t'.join(ped_mendel.loc[key2])+'\n')
#    ped1.close()
#    ped2.close()

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

if 'sample' not in excel.columns:
    if sn.endswith('xlsx'):
        excel = append_sample_excel(excel)
    else:
        excel = append_sample_txt(excel)

#for i in [u'样本编号', u'原始样本ID', u'姓名', u'性别']:
#    if i not in excel.columns:
#        print i
#        print "the header of the excel/csv is not correct, exiting ..."
#        sys.exit()

d1 = append_gender(excel)
d1 = append_relation(d1)
d1 = append_phenotype(d1)
d1 = append_father_mother(d1)
d1 = d1.replace("nan", "")
d1 = d1.fillna('')

pedigree_dict = generate_pedigree(d1) # dictionary of pedigree samples, keys:names of patient samples, values:names of the whole family members including the patient sample
trio = pedigree_dict.keys() 
single = generate_single(d1)

report_dirname = os.getcwd().split('/')[-1]
for i in ['ped', 'peddy', 'annotation', report_dirname]:
    if not os.path.exists(i):
        os.mkdir(i)
if not os.path.exists(report_dirname+'/vcf'):
    os.makedirs(report_dirname+'/vcf')

if len(single) > 0:
    for i in single:
        sample = d1.loc[i]['sample'] # index:姓名
        generate_single_ped(sample,d1)
if len(trio) > 0:
    if not os.path.exists('trio'):
        os.mkdir('trio')
    ofile = open('trio/list', 'w')
    for k in pedigree_dict:
        generate_ped(k,d1)
        sample = d1.loc[k]['sample']
        ofile.write(sample+'\n')
        if software == 'sentieon':
            generate_trio_cmd_sentieon(k, filtmode, pedigree_dict, d1, workpath)
        else:
            generate_trio_cmd_gatk4(k, filtmode, pedigree_dict, d1, workpath)
    ofile.close()

os.system('for i in `cat trio/list`; do cd trio/$i; sbatch trio_cmd_sentieon.sh; cd ..; done')
