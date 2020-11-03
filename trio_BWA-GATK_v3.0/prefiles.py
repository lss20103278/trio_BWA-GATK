#!/usr/bin/env python
# -*- coding: utf-8 -*-
#########################################################################
# File Name: prefiles.py
# Author: sslyu
# mail: 1573077420@qq.com
# Created Time:Fri 19 Apr 2019 10:26:34 AM CST
#########################################################################

import numpy as np
import pandas as pd
import os
import sys
reload(sys)
sys.setdefaultencoding("utf8")
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

def append_gender(df):
    df['gender'] = df[u'性别'].apply(lambda x:'1' if '男' in x else '2' if '女' in x else 'unknown')
    return df                                                                                                                   

def generate_pedigree(df):
    pedigree = {}
    df.index = df[u'姓名']
    for i in df[u'姓名']:
        pedigree[i] = []
    for k in pedigree:
        for i in df[u'姓名']:
            if k in i: # the name of the child should be totally included in other family member's names
                pedigree[k].append(i)
    list_keys = pedigree.keys()
    for k in list_keys:
        if len(pedigree[k]) == 1:
            pedigree.pop(k)
    df.index = df[u'姓名']
    return pedigree

def append_relation(df):
    df.index = df[u'姓名']                                                                                                       
    if 'relationship' not in df.columns:
        df['relationship'] = None
    pedigree_dict = generate_pedigree(df)
    trio = pedigree_dict.keys()
    for i in df.index:
        for j in trio:
            if j in i:
                if j == i:
                    df.loc[i,'relationship'] = u'子'
                elif u'父' in i:
                    df.loc[i,'relationship'] = u'父'  #there is no '父' in i
                elif u'母' in i:
                    df.loc[i,'relationship'] = u'母'  #there is no '母' in i
                elif u'爷' in i:
                    df.loc[i,'relationship'] = u'爷'
                elif u'奶' in i:
                    df.loc[i,'relationship'] = u'奶'
                elif u'外公' in i:
                    df.loc[i,'relationship'] = u'外公'
                elif u'外婆' in i:
                    df.loc[i,'relationship'] = u'外婆'
                elif u'姐' in i:
                    df.loc[i,'relationship'] = u'姐'
                elif u'妹' in i:
                    df.loc[i,'relationship'] = u'妹'
                elif u'哥' in i:
                    df.loc[i,'relationship'] = u'哥'
                elif u'弟' in i:
                    df.loc[i,'relationship'] = u'弟'
                else:
                    df.loc[i,'relationship'] = 'other'
    df[u'家系关系'] = None
    df[u'家系关系'] = df['relationship']
    return df

def append_pedigree(df):
    pedigree_dict = generate_pedigree(df)
    df.index = df[u'姓名']
    df['familyname'] = df['sample']
    for i in df.index:
        for k in pedigree_dict:
            if k in i:  # the name of the child should be totally included in other family member's names
                df.loc[i,'familyname'] = df.loc[k]['sample']
    df[u'样本间关系'] = None
    for i in df.index:
        for k in pedigree_dict:
            if k in i:
                df.loc[i,u'样本间关系'] = 'fam'+df.loc[k]['sample']
    return df                    

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

def append_phenotype(df):
    pedigree_dict = generate_pedigree(df)
    trio = pedigree_dict.keys()
    try:
        single = [i for i in df[u'姓名'] if i not in flatten_list(pedigree_dict.values())] # dict.values()
    except:
        pass
    df['phenotype1'] = '0'
    df['phenotype2'] = None
    if len(trio) > 0 and len(single) > 0:
        for i in df.index: 
            if i in trio:
                df.loc[i,'phenotype2'] = '2'
            elif i in single:
                df.loc[i,'phenotype2'] = '2'
            else:
                df.loc[i,'phenotype2'] = '1'
    elif len(trio) > 0:
        for i in df.index:
            if i in trio:
                df.loc[i,'phenotype2'] = '2'
            else:
                df.loc[i,'phenotype2'] = '1'
    else:
        for i in df.index:
            df.loc[i,'phenotype2'] = '2'                                                                                        
    return df

def generate_single(df):
    pedigree_dict = generate_pedigree(df)
    trio = pedigree_dict.keys()
    trio_m = []
    for i in df[u'姓名']:
        for k in pedigree_dict:
            trio_m.extend(pedigree_dict[k])
    single = [i for i in df[u'姓名'] if i not in trio_m]
    return single

def append_father_mother(df):
    df['father'] = '0'
    df['mother'] = '0'
    pedigree_dict = generate_pedigree(df)
    trio = pedigree_dict.keys()
    try:
        single = [i for i in df[u'姓名'] if i not in pedigree_dict.values()[0]]
    except:                                                                                                                     
        single = df[u'姓名'].tolist()
    df.index = df[u'姓名']
    for i in trio: # i 患病孩子
        relation = ''
        for j in pedigree_dict[i]: 
            relation = relation+df.loc[j]['relationship']
        for j in df.index:
            if j in pedigree_dict[i]: 
                if u'父' in relation and u'母' in relation: #父母都有
                    if j == i: #j为患病孩子
                        for k in pedigree_dict[i]:
                            if u'父' in k:
                                father_name = k 
                            if u'母' in k:
                                mother_name = k
                        df.loc[j,'father'] = df.loc[father_name]['sample']
                        df.loc[j,'mother'] = df.loc[mother_name]['sample']
                    else: # j为患病孩子家系中其他人
                        if len(pedigree_dict[i]) == 3: #家系结构为子父母，父母的父母为0
                            df.loc[j,'father'] = '0' # j为父或母
                        else: #家系结构为子父母和其他家人
                            if u'姐' in j or u'妹' in j or u'哥' in j or u'弟' in j: # 家系结构为 孩子，孩子的兄弟姐妹，父母，j为兄>弟姐妹                                                                                                                          
                                for k in pedigree_dict[i]:
                                    if u'父' in k:
                                        father_name = k
                                    if u'母' in k:
                                        mother_name = k
                                df.loc[j,'father'] = df.loc[father_name]['sample']
                                df.loc[j,'mother'] = df.loc[mother_name]['sample']
                            elif u'爷' not in relation and u'奶' not in relation and u'外婆' not in relation and u'外公' not in relation: 
                                if u'父' in j:
                                    df.loc[j,'father'] = '0'
                                    df.loc[j,'mother'] = '0'
                                if u'母' in j:
                                    df.loc[j,'father'] = '0'
                                    df.loc[j,'mother'] = '0'
                            else:
                                father = raw_input('father ID of '+j+'(if not exist, please enter 0): ')
                                df.loc[j,'father'] = father
                                mother = raw_input('mother ID of '+j+'(if not exist, please enter 0):')
                                df.loc[j,'mother'] = mother
                else: # 父母一方没有
                    if j == i: # j为患病孩子
                        if u'父' in relation: # 只有父亲
                            for k in pedigree_dict[i]:
                                if u'父' in k:
                                    father_name = k
                            df.loc[j,'father'] = df.loc[father_name]['sample']
                        if u'母' in relation: # 只有母亲
                            for k in pedigree_dict[i]:
                                if u'母' in k:
                                    mother_name = k
                            df.loc[j,'mother'] = df.loc[mother_name]['sample']
                    elif u'父' in j: # j为父亲
                        if u'爷' not in relation and u'奶' not in relation and u'外婆' not in relation and u'外公' not in relation: 
                            df.loc[j,'father'] = '0'
                            df.loc[j,'mother'] = '0'
                    elif u'母' in j: # j为母亲
                        if u'爷' not in relation and u'奶' not in relation and u'外婆' not in relation and u'外公' not in relation: 
                            df.loc[j,'father'] = '0'
                            df.loc[j,'mother'] = '0'
                    else: # j为除父母外的其他家人
                        father = raw_input('father ID of '+j+'(if not exist, please enter 0):')
                        df.loc[j,'father'] = father
                        mother = raw_input('mother ID of '+j+'(if not exist, please enter 0):')
                        df.loc[j,'mother'] = mother
    return df                        

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

def generate_cp_redo(sample,ID,rawdir,R1,R2):
    if len(os.listdir(rawdir+'/'+ID)) > 3:
        print "there are more than two fastq files of %s" % ID
        sys.exit()
    cp = open(sample+'/cp_rawdata.sh', 'w')
    cp.write('#!/bin/sh\n\n#SBATCH -J '+sample+'\n#SBATCH -p BP10\n#SBATCH -N 1\n#SBATCH -n 4\n#SBATCH --mail-type=FAIL\n#SBATCH --mail-user=1573077420@qq.com\n\n')
    cp.write('sample='+sample+'\nID='+ID+'\n\ndir='+rawdir+'\nR1='+R1+'\nR2='+R2+'\n\n')
    cp.write('cp $dir/$ID/*$ID*$R1 ../raw/$sample\_R1.fastq.gz\n')
    cp.write('cp $dir/$ID/*$ID*$R2 ../raw/$sample\_R2.fastq.gz\n')
    cp.write('cat $dir/$ID/$ID.md5 >> ../backup_md5\n')
    cp.write('for j in $R1 $R2; do md5sum $dir/$ID/*$ID*$j >> ../original_md5; done\n')
    cp.write('for j in R1 R2; do md5sum ../raw/$sample\_$j*.fastq.gz >> ../cp_md5; done\n')
    cp.write('for j in $R1 $R2; do original=`grep $ID ../original_md5 |grep $j |cut -d" " -f 1`; backup=`grep $ID ../backup_md5 |grep $j |cut -d" " -f 1`; if [ "$original" != "$backup" ]; then echo -e $ID"_"$j" is not backup correctly" >> ../check_md5.log; echo -e $ID"_"$j" is not backup correctly"; fi; done\n')
    cp.write('pattern_original=($R1 $R2); pattern_copied=(R1 R2); for j in 0 1; do original=`grep $ID ../original_md5 |grep ${pattern_original[$j]} |cut -d" " -f 1`; copied=`grep $sample ../cp_md5 |grep ${pattern_copied[$j]} |cut -d" " -f 1`; if [ "$original" != "$copied" ]; then echo -e $sample"_"$ID"_"${pattern_original[$j]}" is not copied correctly" >> ../check_md5.log; echo -e $sample"_"$ID"_"${pattern_original[$j]}" is not copied correctly"; fi; done\n')
    cp.write('curl \'https://oapi.dingtalk.com/robot/send?access_token=dc5bf30a873950e8c73bb37abf1851ebaf5d83f2792f08684ace2be621a12563\' -H \'Content-Type: application/json\' -d \'{"msgtype": "text","text": {"content": "'+sample+' copy rawdata done"}}\'\n')                                                                                                                              
    cp.close()

def generate_cp(sample,ID,rawdir,md5,R1,R2,panel,outdir):                                                                       
    cp = open(sample+'/cp_rawdata.sh', 'w')
    cp.write('#!/bin/sh\n\n#SBATCH -J '+sample+'\n#SBATCH -p BP10\n#SBATCH -N 1\n#SBATCH -n 4\n#SBATCH --mail-type=FAIL\n#SBATCH --mail-user=1573077420@qq.com\n\n') #--mail-type=END,FAIL
    cp.write('sample='+sample+'\nID='+ID+'\n\ndir='+rawdir+'\nmd5='+md5+'\nR1='+R1+'\nR2='+R2+'\npanel='+panel+'\n\n')
    cp.write('outdir='+outdir+'\n')
    #if os.path.exists(outdir+'/'+ID):
    #    print "%s/%s exists ..." % (outdir, ID)
    #    sys.exit()
    cp.write('mkdir -p $outdir/$ID\n')
    cp.write('awk \'{if(match($2,"\'$ID\'")){print $0}}\' $md5 > $outdir/$ID/$ID.md5\n')
    cp.write('cp $dir/*$ID*$R1 $outdir/$ID\n')
    cp.write('cp $dir/*$ID*$R2 $outdir/$ID\n')
    cp.write('cp $dir/*$ID*$R1 ../raw/$sample\_R1.fastq.gz\n')
    cp.write('cp $dir/*$ID*$R2 ../raw/$sample\_R2.fastq.gz\n')

    cp.write('for j in $R1 $R2; do md5sum $dir/*$ID*$j >> ../original_md5; done\n')
    cp.write('for j in R1 R2; do md5sum ../raw/$sample\_$j*.fastq.gz >> ../cp_md5; done\n')
    cp.write('for j in $R1 $R2; do md5sum $outdir/$ID/*$ID*$j >> ../backup_md5; done\n')
    cp.write('for j in $R1 $R2; do a=`grep $ID ../original_md5 |grep $j |cut -d" " -f 1`; b=`grep $ID $md5 |grep $j |cut -d" " -f 1`; if [ "$a" != "$b" ]; then echo -e "original "$ID"_"$j" is problematic" >> ../check_md5.log; echo -e "original "$ID"_"$j" is problematic"; fi; done\n')
    cp.write('pattern_original=($R1 $R2); pattern_copied=(R1 R2); for j in 0 1; do original=`grep $ID ../original_md5 |grep ${pattern_original[$j]} |cut -d" " -f 1`; copied=`grep $sample ../cp_md5 |grep ${pattern_copied[$j]} |cut -d" " -f 1`; if [ "$original" != "$copied" ]; then echo -e $sample"_"$ID"_"${pattern_original[$j]}" is not copied correctly" >> ../check_md5.log; echo -e $sample"_"$ID"_"${pattern_original[$j]}" is not copied correctly"; fi; done\n')
    cp.write('for j in $R1 $R2; do original=`grep $ID ../original_md5 |grep $j |cut -d" " -f 1`; backup=`grep $ID ../backup_md5 |grep $j |cut -d" " -f 1`; if [ "$original" != "$backup" ]; then echo -e $ID"_"$j" is not backup correctly" >> ../check_md5.log; echo -e $ID"_"$j" is not backup correctly"; fi; done\n')
    cp.write('curl \'https://oapi.dingtalk.com/robot/send?access_token=dc5bf30a873950e8c73bb37abf1851ebaf5d83f2792f08684ace2be621a12563\' -H \'Content-Type: application/json\' -d \'{"msgtype": "text","text": {"content": "'+sample+' copy rawdata done"}}\'\n')
    cp.close()
    
def generate_cmd_gatk4_part1(sample,exon):
    global script_path
    single_cmd = open(sample+'/cmd.sh', 'w')
    single_cmd.write('#!/bin/sh\n\n#SBATCH -J '+sample+'\n#SBATCH -p BP10\n#SBATCH -N 1\n#SBATCH -n 4\n#SBATCH --mail-type=FAIL\n#SBATCH --mail-user=1573077420@qq.com\n\n')
    single_cmd.write('sample='+sample+'\n\n')
    single_cmd.write('sh run2_gatk4.sh $sample $SLURM_NPROCS\n')
    single_cmd.write('/home/ana005/anaconda2/bin/iTools Fqtools stat -InFq ../raw/'+sample+'_R1.fastq.gz -InFq ../raw/'+sample+'_R2.fastq.gz -OutStat '+sample+'.info\n')
    if len(exon) > 0 and exon[0] != "":
        for i in exon:
            if os.path.exists(script_path+'/bed/gene/'+i+'.bed'):
                single_cmd.write('samtools depth -b %s/bed/gene/%s.bed /BP12_share/sslyu/bam/%s.sort.mkdup.bam > ../exon/%s_%s.bed.depth\n' % (script_path,i,sample,sample,i))
            else:
                print "%s/bed/gene/%s.bed doesn't exist" % (script_path,i)
                sys.exit()
        single_cmd.write('cd ../exon\n')
        for i in exon:
            single_cmd.write('python /DATA/sslyu/trio_BWA-GATK_v3.0/src/exon_depth_barplot.py -gene %s -sn %s\n' % (i, sample))
        single_cmd.write('cd ../%s\n' % sample)
    single_cmd.write('scp /BP12_share/sslyu/bam/'+sample+'.sort.mkdup.bam* klyang@122.112.248.194:/media/bpdata/data/exon_temp2/\n')
    single_cmd.write('curl \'https://oapi.dingtalk.com/robot/send?access_token=dc5bf30a873950e8c73bb37abf1851ebaf5d83f2792f08684ace2be621a12563\' -H \'Content-Type: application/json\' -d \'{"msgtype": "text","text": {"content": "'+sample+' mapping done"}}\'\n')
    single_cmd.close()

def generate_cmd_sentieon_part1(sample,exon):
    global script_path
    single_cmd = open(sample+'/cmd.sh', 'w')
    single_cmd.write('#!/bin/sh\n\n#SBATCH -J '+sample+'\n#SBATCH -p BP10\n#SBATCH -N 1\n#SBATCH -n 4\n#SBATCH --mail-type=FAIL\n#SBATCH --mail-user=1573077420@qq.com\n\n')
    single_cmd.write('sample='+sample+'\n\n')
    single_cmd.write('sh run2_sentieon.sh $sample $SLURM_NPROCS\n')
    single_cmd.write('/home/ana005/anaconda2/bin/iTools Fqtools stat -InFq ../raw/'+sample+'_R1.fastq.gz -InFq ../raw/'+sample+'_R2.fastq.gz -OutStat '+sample+'.info\n')
    if len(exon) > 0 and exon[0] != "":
        for i in exon:
            if os.path.exists(script_path+'/bed/gene/'+i+'.bed'):
                single_cmd.write('samtools depth -b %s/bed/gene/%s.bed /BP12_share/sslyu/bam/%s.sort.mkdup.bam > ../exon/%s_%s.bed.depth\n' % (script_path,i,sample,sample,i))
            else:
                print "%s/bed/gene/%s.bed doesn't exist" % (script_path,i)
                sys.exit()
        single_cmd.write('cd ../exon\n')
        for i in exon:
            single_cmd.write('python /DATA/sslyu/trio_BWA-GATK_v3.0/src/exon_depth_barplot.py -gene %s -sn %s\n' % (i, sample))
        single_cmd.write('cd ../%s\n' % sample)
    single_cmd.write('scp /BP12_share/sslyu/bam/'+sample+'.sort.mkdup.bam* klyang@122.112.248.194:/media/bpdata/data/exon_temp2/\n')
    single_cmd.write('curl \'https://oapi.dingtalk.com/robot/send?access_token=dc5bf30a873950e8c73bb37abf1851ebaf5d83f2792f08684ace2be621a12563\' -H \'Content-Type: application/json\' -d \'{"msgtype": "text","text": {"content": "'+sample+' mapping done"}}\'\n')
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
    single_cmd.write('curl \'https://oapi.dingtalk.com/robot/send?access_token=dc5bf30a873950e8c73bb37abf1851ebaf5d83f2792f08684ace2be621a12563\' -H \'Content-Type: application/json\' -d \'{"msgtype": "text","text": {"content": "'+sample+' annotation excel done"}}\'\n')
    single_cmd.close()
    ofile = open('annotation/annotation.sh', 'a')
    ofile.write('python '+script_path+'/src/annotation_filt.py -sn '+sample+' -c_f_m c\n')
    ofile.close()

def generate_trio_cmd_gatk4(k, filtmode, pedigree_dict, d1, dirname):
    sample = d1.loc[k]['sample']
    path = 'trio/'+sample
    if not os.path.exists(path):
        os.makedirs(path)
    trio_cmd = open(path+'/trio_cmd_gatk4.sh', 'w')
    #trio_cmd.write('#!/bin/sh\n\n#SBATCH -J '+sample+'\n#SBATCH -p BP10\n#SBATCH -N 1\n#SBATCH -n 4\n#SBATCH --mail-type=END,FAIL\n#SBATCH --mail-user=1573077420@qq.com\n\nsample='+sample+'\n\n')
    trio_cmd.write('#!/bin/sh\n\n#SBATCH -J '+sample+'\n#SBATCH -p BP10\n#SBATCH -N 1\n#SBATCH -n 4\n#SBATCH --mail-type=FAIL\n#SBATCH --mail-user=1573077420@qq.com\n\nsample='+sample+'\n\n')
    trio_cmd.write('filtmode=\''+filtmode+'\'\npath=`pwd`\nvar=`echo $path | awk -F \'/\' \'{print $NF}\'`\n')
    trio_cmd.write('if [ $var = \'peddy\' ]\nthen\nsh '+script_path+'/src/trio_gatk4.sh $sample $SLURM_NPROCS\n')
    trio_cmd.write('curl \'https://oapi.dingtalk.com/robot/send?access_token=dc5bf30a873950e8c73bb37abf1851ebaf5d83f2792f08684ace2be621a12563\' -H \'Content-Type: application/json\' -d \'{"msgtype": "text","text": {"content": "peddy done"}}\'\n')
    trio_cmd.write('else\nsh '+script_path+'/src/trio_gatk4.sh $sample $SLURM_NPROCS $filtmode\n')
    trio_cmd.write('[ -e ../../$sample/2_mapping/$sample.depth.sample_gene_summary ] && cp ../../$sample/2_mapping/$sample.depth.sample_gene_summary ../../annotation\n')
    trio_cmd.write('cp segtrio/$sample.{ann.hg19_multianno.txt,CADD,maf,link} ../../annotation\n')
    trio_cmd.write('cd ../../annotation\n')
    trio_cmd.write('echo $sample >> list\n')
    trio_cmd.write('python '+script_path+'/src/score_re.py -sn $sample\n')
    relation = []
    for j in pedigree_dict[k]:
        relation.append(d1.loc[j]['relationship'])
    a = ' '.join(relation)
    print a
    annotationfile = open('annotation/annotation.sh', 'a')
    if u'父' in a and u'母' in a:
        trio_cmd.write('python '+script_path+'/src/annotation_filt.py -sn $sample -c_f_m c_f_m\n')
        annotationfile.write('python '+script_path+'/src/annotation_filt.py -sn '+sample+' -c_f_m c_f_m\n')
    else:
        if u'父' in a:
            trio_cmd.write('python '+script_path+'/src/annotation_filt.py -sn $sample -c_f_m c_f\n')
            annotationfile.write('python '+script_path+'/src/annotation_filt.py -sn '+sample+' -c_f_m c_f\n')
        else:
            trio_cmd.write('python '+script_path+'/src/annotation_filt.py -sn $sample -c_f_m c_m\n')
            annotationfile.write('python '+script_path+'/src/annotation_filt.py -sn '+sample+' -c_f_m c_m\n')
    trio_cmd.write('mv '+sample+'_children_hospital*xlsx ../'+dirname[-1]+'\n')
    trio_cmd.write('cp ../trio/'+sample+'/sep/*.vcf ../'+dirname[-1]+'/vcf\n')
    trio_cmd.write('echo "end `date`" >> ../$sample/finished\n')
    trio_cmd.write('curl \'https://oapi.dingtalk.com/robot/send?access_token=dc5bf30a873950e8c73bb37abf1851ebaf5d83f2792f08684ace2be621a12563\' -H \'Content-Type: application/json\' -d \'{"msgtype": "text","text": {"content": "'+sample+' annotation excel done"}}\'\n')
    trio_cmd.write('fi\n')
    trio_cmd.close()
    annotationfile.close()

def generate_trio_cmd_sentieon(k, filtmode, pedigree_dict, d1, dirname):              
    sample = d1.loc[k]['sample']
    path = 'trio/'+sample
    if not os.path.exists(path):
        os.makedirs(path)
    trio_cmd = open(path+'/trio_cmd_sentieon.sh', 'w')
    #trio_cmd.write('#!/bin/sh\n\n#SBATCH -J '+sample+'\n#SBATCH -p BP10\n#SBATCH -N 1\n#SBATCH -n 4\n#SBATCH --mail-type=END,FAIL\n#SBATCH --mail-user=1573077420@qq.com\n\nsample='+sample+'\n\n')
    trio_cmd.write('#!/bin/sh\n\n#SBATCH -J '+sample+'\n#SBATCH -p BP10\n#SBATCH -N 1\n#SBATCH -n 4\n#SBATCH --mail-type=FAIL\n#SBATCH --mail-user=1573077420@qq.com\n\nsample='+sample+'\n\n')
    trio_cmd.write('filtmode=\''+filtmode+'\'\npath=`pwd`\nvar=`echo $path | awk -F \'/\' \'{print $NF}\'`\n')
    trio_cmd.write('if [ $var = \'peddy\' ]\nthen\nsh '+script_path+'/src/trio_sentieon.sh $sample $SLURM_NPROCS\n')
    trio_cmd.write('curl \'https://oapi.dingtalk.com/robot/send?access_token=dc5bf30a873950e8c73bb37abf1851ebaf5d83f2792f08684ace2be621a12563\' -H \'Content-Type: application/json\' -d \'{"msgtype": "text","text": {"content": "peddy done"}}\'\n')
    trio_cmd.write('else\nsh '+script_path+'/src/trio_sentieon.sh $sample $SLURM_NPROCS $filtmode\n')
    trio_cmd.write('[ -e ../../$sample/2_mapping/$sample.depth.sample_gene_summary ] && cp ../../$sample/2_mapping/$sample.depth.sample_gene_summary ../../annotation\n')
    trio_cmd.write('cp segtrio/$sample.{ann.hg19_multianno.txt,CADD,maf,link} ../../annotation\n')
    trio_cmd.write('cd ../../annotation\n')
    trio_cmd.write('echo $sample >> list\n')
    trio_cmd.write('python '+script_path+'/src/score_re.py -sn $sample\n')
    relation = []
    for j in pedigree_dict[k]:
        relation.append(d1.loc[j]['relationship'])
    a = ' '.join(relation)
    print a
    annotationfile = open('annotation/annotation.sh', 'a')
    if u'父' in a and u'母' in a:                                                                                               
        trio_cmd.write('python '+script_path+'/src/annotation_filt.py -sn $sample -c_f_m c_f_m\n')
        annotationfile.write('python '+script_path+'/src/annotation_filt.py -sn '+sample+' -c_f_m c_f_m\n')
    else:
        if u'父' in a:
            trio_cmd.write('python '+script_path+'/src/annotation_filt.py -sn $sample -c_f_m c_f\n')
            annotationfile.write('python '+script_path+'/src/annotation_filt.py -sn '+sample+' -c_f_m c_f\n')
        else:
            trio_cmd.write('python '+script_path+'/src/annotation_filt.py -sn $sample -c_f_m c_m\n')
            annotationfile.write('python '+script_path+'/src/annotation_filt.py -sn '+sample+' -c_f_m c_m\n')
    trio_cmd.write('mv '+sample+'_children_hospital*xlsx ../'+dirname[-1]+'\n')
    trio_cmd.write('cp ../trio/'+sample+'/sep/*.vcf ../'+dirname[-1]+'/vcf\n')
    trio_cmd.write('echo "end `date`" >> ../$sample/finished\n')
    trio_cmd.write('curl \'https://oapi.dingtalk.com/robot/send?access_token=dc5bf30a873950e8c73bb37abf1851ebaf5d83f2792f08684ace2be621a12563\' -H \'Content-Type: application/json\' -d \'{"msgtype": "text","text": {"content": "'+sample+' annotation excel done"}}\'\n')
    trio_cmd.write('fi\n')
    trio_cmd.close()
    annotationfile.close()

def generate_single_ped(k,d):
    d.index = d['sample']
    ped_mendel = d.loc[k][['sample', 'father', 'mother', 'gender', 'phenotype2']]
    ped_mendel = ped_mendel.apply(str)
    fname = d.loc[k]['sample']+'/'+d.loc[k]['sample']+'.mendel.ped'
    fname = 'ped/'+d.loc[k]['sample']+'.mendel.ped'
    ped = open(fname, 'w')
    ped.write(str(d.loc[k]['sample'])+'\t'+'\t'.join(ped_mendel)+'\n')
    ped.close()
    d.index = d[u'姓名']

def generate_ped(k,d):
    pedigree_dict = generate_pedigree(d)
    ped_f = d.loc[pedigree_dict[k]][['sample', 'father', 'mother', 'gender', 'phenotype1']]
    ped_f.index = ped_f['sample']
    ped_mendel = d.loc[pedigree_dict[k]][['sample', 'father', 'mother', 'gender', 'phenotype2']]
    ped_mendel.index = ped_mendel['sample']
    f1name = 'ped/'+d.loc[k]['sample']+'.ped'
    f2name = 'ped/'+d.loc[k]['sample']+'.mendel.ped'
    ped1 = open(f1name, 'w')
    ped2 = open(f2name, 'w')
    relation = []
    for i in pedigree_dict[k]:
        relation.append(d.loc[i]['relationship'])
    a = ' '.join(relation)
    key1 = d.loc[k]['sample']
    ped1.write(key1+'\t'+'\t'.join(ped_f.loc[key1])+'\n')
    ped2.write(key1+'\t'+'\t'.join(ped_mendel.loc[key1])+'\n')
    if '父' in a and '母' in a:
        for j in ['父', '母']:
            for i in pedigree_dict[k]:
                if i != k and j in i:
                    key2 = d.loc[i]['sample']
                    ped1.write(key1+'\t'+'\t'.join(ped_f.loc[key2])+'\n')
                    ped2.write(key1+'\t'+'\t'.join(ped_mendel.loc[key2])+'\n')
        n_c_f_m = [i for i in pedigree_dict[k] if '父' not in i and '母' not in i and i != k]
        for i in n_c_f_m:
            key2 = d.loc[i]['sample']
            ped1.write(key1+'\t'+'\t'.join(ped_f.loc[key2])+'\n')                                                               
            ped2.write(key1+'\t'+'\t'.join(ped_mendel.loc[key2])+'\n')
    elif '父' in a or '母' in a:
        parent = ''
        if '父' in a:
            parent = '父'
        else:
            parent = '母'
        for i in pedigree_dict[k]:
            if parent in i and i != k:
                key2 = d.loc[i]['sample']
                ped1.write(key1+'\t'+'\t'.join(ped_f.loc[key2])+'\n')
                ped2.write(key1+'\t'+'\t'.join(ped_mendel.loc[key2])+'\n')
        n_c_f_m = [i for i in pedigree_dict[k] if parent not in i and i != k]
        for i in n_c_f_m:
            key2 = d.loc[i]['sample']
            ped1.write(key1+'\t'+'\t'.join(ped_f.loc[key2])+'\n')
            ped2.write(key1+'\t'+'\t'.join(ped_mendel.loc[key2])+'\n')
    else:
        n_c = [i for i in pedigree_dict[k] if i != k]
        for i in n_c:
            key2 = d.loc[i]['sample']
            ped1.write(key1+'\t'+'\t'.join(ped_f.loc[key2])+'\n')
            ped2.write(key1+'\t'+'\t'.join(ped_mendel.loc[key2])+'\n')
    ped1.close()
    ped2.close()

if len(sys.argv) > 1:
    cmd = sys.argv[1]

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

if sn.endswith('xlsx'):
    excel = append_sample_excel(excel)
else:
    excel = append_sample_txt(excel)
for i in [u'样本编号', u'原始样本ID', u'姓名', u'性别']:
    if i not in excel.columns:
        print i
        print "the header of the excel/csv is not correct, exiting ..."
        sys.exit()

d1 = append_gender(excel)
d1 = append_relation(d1)
d1 = append_pedigree(d1)
d1 = append_phenotype(d1)
d1 = append_father_mother(d1)
d1 = d1.replace("nan", "")
d1 = d1.fillna('')

pedigree_dict = generate_pedigree(d1) # dictionary of pedigree samples, keys:names of patient samples, values:names of the whole family members including the patient sample
trio = pedigree_dict.keys() 
single = generate_single(d1)

generate_dbevn(panel)

report_dirname = os.getcwd().split('/')[-1]
for i in ['ped', 'peddy', 'gvcf', 'annotation', 'raw', 'CNV', report_dirname, 'exon']:
    if not os.path.exists(i):
        os.mkdir(i)
if not os.path.exists(report_dirname+'/vcf'):
    os.makedirs(report_dirname+'/vcf')
if d1.shape[1] > 1:
    if not os.path.exists(report_dirname+'/CNV'):
        os.makedirs(report_dirname+'/CNV')

if u'样本编号' in d1.columns and u'原始样本ID' in d1.columns:
    d1_tmp = d1[[u'样本编号', u'原始样本ID']]
d1_tmp.to_csv('rename', index=None, header=None, sep="\t")

d1['sample'].to_csv('list', header=None, index=None)
d1['sample'].to_csv('CNV/id', header=None, index=None)
d1['CNV'] = '/BP12_share/sslyu/bam/'+d1['sample']+'.sort.mkdup.bam'
d1['CNV'].to_csv('CNV/list', header=None, index=None)                                                

for i in d1['sample']:                                                                                                      
    if not os.path.exists(str(i)):
        os.mkdir(str(i))
for i in ['run2_sentieon.sh', 'cp_rawdata.sh']:
    shutil.copy(script_path+'/src/'+i,os.getcwd())

run2_sentieon = open(script_path+'/src/run2_sentieon.sh').readlines()
for i in range(28,37):
    print run2_sentieon[i].strip('\n')

if os.path.exists('annotation/annotation.sh'):
    os.remove('annotation/annotation.sh')
if len(single) > 0:
    shutil.copy(script_path+'/src/cmd_single_sentieon.sh',os.getcwd())
    list_single = open('list_single', 'w')
    if panel == "WES":
        cmd_single = "clean mapping precalling gvcf index_single gtgvcf vqsr left-normalize depth qcsum qcvcf annotation"
    else:
        cmd_single = "clean mapping precalling gvcf index_single gtgvcf hardfilt left-normalize depth qcsum qcvcf annotation"
    for i in single:
        sample = d1.loc[i]['sample'] # index:姓名
        print 'single '+sample
        ID = d1.loc[i][u'样本编号']
        list_single.write(sample+'\n')
        if os.path.exists(sample+'/cmd.sh'):
            os.remove(sample+'/cmd.sh')
        if redo == 'yes':
            generate_cp_redo(sample, ID, rawdir, R1, R2)
        else:
            generate_cp(sample, ID, rawdir, md5, R1, R2, panel, outdir)
        exon = d1.loc[i]['exon'].split(',')
        print exon
        cmd_single = d1.loc[i]['cmd']
        #generate_cmd_sentieon_part1(sample, exon)
        #generate_run2_sentieon(sample,cmd_single)
        generate_cmd_gatk4_part1(sample,exon)
        generate_run2_gatk4(sample,cmd_single)
        generate_cmd_part2(sample, workpath)
        generate_single_ped(sample,d1)
    list_single.close()
if len(trio) > 0:
    shutil.copy(script_path+'/src/cmd_trio_sentieon.sh',os.getcwd())
    shutil.copy(script_path+'/src/trio_cmd_sentieon.sh', os.getcwd())
    if not os.path.exists('trio'):
        os.mkdir('trio')
    cmd_trio = "clean mapping precalling gvcf index_trio depth qcsum"
    list_trio = open('list_trio', 'w')
    ofile = open('trio/list', 'w')
    for k in pedigree_dict:
        generate_ped(k,d1)
        sample = d1.loc[k]['sample']
        ofile.write(sample+'\n')
        #generate_trio_cmd_sentieon(k, filtmode, pedigree_dict, d1, workpath)
        generate_trio_cmd_gatk4(k, filtmode, pedigree_dict, d1, workpath)
        for i in pedigree_dict[k]:
            sample = d1.loc[i]['sample']
            print 'trio '+sample
            ID = str(d1.loc[i][u'样本编号'])
            list_trio.write(sample+'\n')
            if redo == 'yes':
                generate_cp_redo(sample, ID, rawdir, R1, R2)
            else:
                generate_cp(sample, ID, rawdir, md5, R1, R2, panel, outdir)
            exon = d1.loc[i]['exon'].split(',')
            cmd_trio = d1.loc[i]['cmd']
            print cmd_trio
            #generate_cmd_sentieon_part1(sample,exon)
            #generate_run2_sentieon(sample,cmd_trio)
            generate_cmd_gatk4_part1(sample,exon)
            generate_run2_gatk4(sample,cmd_trio)
    list_trio.close()
    ofile.close()
if not os.path.exists("annotation/annotation.sh"):
    print """
    annotation/annotation.sh is missing
    execute -prepare_sentieon or -prepare_gatk4 first
    """
    sys.exit()
for i in single:
    if not os.path.exists("ped/"+d1.loc[i]['sample']+".mendel.ped"):
        print """
        ped/%s.mendel.ped is missing
        execute -prepare_sentieon or -prepare_gatk4 first
        """ % (d1.loc[i]['sample'])
        sys.exit()
for i in trio:
    if not os.path.exists("ped/"+d1.loc[i]['sample']+".mendel.ped"):
        print """
        ped/%s.mendel.ped is missing
        execute -prepare_sentieon or -prepare_gatk4 first
        """ % (d1.loc[i]['sample'])
        sys.exit()
os.system('sort annotation/annotation.sh |uniq > annotation/tmp; mv annotation/tmp annotation/annotation.sh')        
if os.path.exists('list_single'):
    os.system('for i in `cat list_single`; do cat ped/$i.mendel.ped >> peddy_tmp; done')
if os.path.exists('trio/list'):
    os.system('for i in `cat trio/list`; do cat ped/$i.mendel.ped >> peddy_tmp; done')
os.system('sort peddy_tmp |uniq |awk \'{OFS="\t"; print \'0\',$2,$3,$4,$5,$6}\' > ped/peddy.ped; rm peddy_tmp')
with open('ped/peddy.ped') as f:
    for l in f:                                                                                                             
        print l.strip('\n')

