#!/usr/bin/env python
# -*- coding: utf-8 -*-
#########################################################################
# File Name: cmdfile.py
# Author: sslyu
# mail: 1573077420@qq.com
# Created Time:Fri 19 Apr 2019 03:55:09 PM CST
#########################################################################

import numpy as np
import pandas as pd
import os
import sys
reload(sys)
sys.setdefaultencoding("utf8")
import shutil

m_path = os.path.split(os.path.abspath(sys.argv[0]))[0]
print m_path

def generate_dbevn(strategy):
    dbevn = open('dbevn.sh', 'w')
    if strategy != 'none':
        with open(m_path+'/src/dbevn.sh') as f:
            for l in f:
                if l.startswith('#panel'):
                    if strategy in l:                                                                                           
                        dbevn.write(l[1:])
                    else:
                        dbevn.write(l)
                else:
                    dbevn.write(l)
    else:
        with open(m_path+'/src/dbevn.sh') as f:
            for l in f:
                dbevn.write(l)
        panel = raw_input('please enter the absolute path of the bed:')
        dbevn.write('panel=\"'+panel+'\"\n')
    dbevn.close()

def generate_run2_sentieon(path,str):
    global m_path
    run2_sentieon = open(path+'/run2_sentieon.sh', 'w')
    n = 0
    with open(m_path+'/src/run2_sentieon.sh') as f:
        for l in f:
            n = n+1
            if n>0 and n<= 16:
                run2_sentieon.write(l)
            elif n>17 and n<24:
                run2_sentieon.write(l)
            elif n>= 25:
                run2_sentieon.write(l)
            elif n==24:
                run2_sentieon.write('cmd="'+str+'"\n')
            elif n==17:
                run2_sentieon.write('basedir='+m_path+'/src\n')
            else:
                continue
    run2_sentieon.close()                                                                                                       

def generate_cp(sample,ID,dir,md5,R1,R2,strategy,outdir):
    cp = open(sample+'/cp_rawdata.sh', 'w')
    #cp.write('#!/bin/sh\n\n#SBATCH -J '+sample+'\n#SBATCH -p BP10\n#SBATCH -N 1\n#SBATCH -n 4\n#SBATCH --mail-type=END,FAIL\n#SBATCH --mail-user=1573077420@qq.com\n\n')
    cp.write('#!/bin/sh\n\n#SBATCH -J '+sample+'\n#SBATCH -p BP10\n#SBATCH -N 1\n#SBATCH -n 4\n#SBATCH --mail-type=FAIL\n#SBATCH --mail-user=1573077420@qq.com\n\n')                                                                                            
    cp.write('sample='+sample+'\nID='+ID+'\n\ndir='+dir+'\nmd5='+md5+'\nR1='+R1+'\nR2='+R2+'\nstrategy='+strategy+'\n\n')
    #cp.write('mkdir -p /anaData/anaData004/children_hos_genetic/rawdata/2018/$strategy/$ID\n')
    #cp.write('awk \'{if(match($2,"\'$ID\'")){print $0}}\' $md5 >> /anaData/anaData004/children_hos_genetic/rawdata/2018/$strategy/$ID/$ID.md5\n')
    #cp.write('cp $dir/*$ID*$R1 /anaData/anaData004/children_hos_genetic/rawdata/2018/$strategy/$ID\n')
    #cp.write('cp $dir/*$ID*$R2 /anaData/anaData004/children_hos_genetic/rawdata/2018/$strategy/$ID\n')
    cp.write('outdir='+outdir+'\n')
    cp.write('mkdir -p $outdir/$ID\n')
    cp.write('awk \'{if(match($2,"\'$ID\'")){print $0}}\' $md5 >> $outdir/$ID/$ID.md5\n')
    cp.write('cp $dir/*$ID*$R1 $outdir/$ID\n')
    cp.write('cp $dir/*$ID*$R2 $outdir/$ID\n')
    cp.write('cp $dir/*$ID*$R1 ../raw/$sample\_R1.fastq.gz\n')
    cp.write('cp $dir/*$ID*$R2 ../raw/$sample\_R2.fastq.gz\n')

    cp.write('for j in $R1 $R2; do md5sum $dir/*$ID*$j >> ../original_md5; done\n')
    cp.write('for j in R1 R2; do md5sum ../raw/$sample\_$j*.fastq.gz >> ../cp_md5; done\n')
    #cp.write('for j in $R1 $R2; do md5sum /anaData/anaData004/children_hos_genetic/rawdata/2018/$strategy/$ID/*$ID*$j >> ../backup_md5; done\n')
    cp.write('for j in $R1 $R2; do md5sum $outdir/$ID/*$ID*$j >> ../backup_md5; done\n')
    cp.write('for j in $R1 $R2; do a=`grep $ID ../original_md5 |grep $j |cut -d" " -f 1`; b=`grep $ID $md5 |grep $j |cut -d" " -f 1`; if [ "$a" != "$b" ]; then echo -e "original "$ID"_"$j" is problematic" >> ../check_md5.log; echo -e "original "$ID"_"$j" is problematic"; fi; done\n')
    cp.write('pattern_original=($R1 $R2); pattern_copied=(R1 R2); for j in 0 1; do original=`grep $ID ../original_md5 |grep ${pattern_original[$j]} |cut -d" " -f 1`; copied=`grep $sample ../cp_md5 |grep ${pattern_copied[$j]} |cut -d" " -f 1`; if [ "$original" != "$copied" ]; then echo -e $sample"_"$ID"_"${pattern_original[$j]}" is not copied correctly" >> ../check_md5.log; echo -e $sample"_"$ID"_"${pattern_original[$j]}" is not copied correctly"; fi; done\n')
    cp.write('for j in $R1 $R2; do original=`grep $ID ../original_md5 |grep $j |cut -d" " -f 1`; backup=`grep $ID ../backup_md5 |grep $j |cut -d" " -f 1`; if [ "$original" != "$backup" ]; then echo -e $ID"_"$j" is not backup correctly" >> ../check_md5.log; echo -e $ID"_"$j" is not backup correctly"; fi; done\n')
    cp.write('curl \'https://oapi.dingtalk.com/robot/send?access_token=dc5bf30a873950e8c73bb37abf1851ebaf5d83f2792f08684ace2be621a12563\' -H \'Content-Type: application/json\' -d \'{"msgtype": "text","text": {"content": "'+sample+' copy rawdata done"}}\'\n')
    cp.close()
    
def generate_cmd_sentieon_part1(sample):
    single_cmd = open(sample+'/cmd.sh', 'w')
    #single_cmd.write('#!/bin/sh\n\n#SBATCH -J '+sample+'\n#SBATCH -p BP10\n#SBATCH -N 1\n#SBATCH -n 4\n#SBATCH --mail-type=END,FAIL\n#SBATCH --mail-user=1573077420@qq.com\n\n')
    single_cmd.write('#!/bin/sh\n\n#SBATCH -J '+sample+'\n#SBATCH -p BP10\n#SBATCH -N 1\n#SBATCH -n 4\n#SBATCH --mail-type=FAIL\n#SBATCH --mail-user=1573077420@qq.com\n\n')
    single_cmd.write('sample='+sample+'\n\n')
    single_cmd.write('sh run2_sentieon.sh $sample $SLURM_NPROCS\n')
    #single_cmd.write('samtools depth -b /DATA/sslyu/Project/Genetics_children_hospital/all_ch_freq_gene.bed /BP12_share/sslyu/bam/'+sample+'.sort.mkdup.bam > ../depth/'+sample+'.depth\n')
    single_cmd.write('scp /BP12_share/sslyu/bam/'+sample+'.sort.mkdup.bam* klyang@122.112.248.194:/media/bpdata/data/exon_temp2/\n')                                                                                                                            
    single_cmd.write('curl \'https://oapi.dingtalk.com/robot/send?access_token=dc5bf30a873950e8c73bb37abf1851ebaf5d83f2792f08684ace2be621a12563\' -H \'Content-Type: application/json\' -d \'{"msgtype": "text","text": {"content": "'+sample+' mapping done"}}\'\n')
    single_cmd.close()

def generate_cmd_part2(sample, dirname):
    single_cmd = open(sample+'/cmd.sh', 'a')
    single_cmd.write('cp 2_mapping/$sample.depth.sample_gene_summary ../annotation\n')
    single_cmd.write('cp 3_variants/$sample.{ann.hg19_multianno.txt,CADD,maf,link} ../annotation\n')
    single_cmd.write('cd ../annotation\n')
    single_cmd.write('echo $sample >> list\n')
    single_cmd.write('python '+m_path+'/src/score_re.py -sn $sample\n')
    single_cmd.write('python '+m_path+'/src/annotation_filt.py -sn $sample -c_f_m c\n')
    single_cmd.write('mv '+sample+'_children_hospital*xlsx ../'+dirname[-1]+'\n')
    single_cmd.write('cp ../'+sample+'/3_variants/'+sample+'.vcf ../'+dirname[-1]+'/vcf\n')
    single_cmd.write('echo "end `date`" >> ../$sample/finished\n')
    single_cmd.write('curl \'https://oapi.dingtalk.com/robot/send?access_token=dc5bf30a873950e8c73bb37abf1851ebaf5d83f2792f08684ace2be621a12563\' -H \'Content-Type: application/json\' -d \'{"msgtype": "text","text": {"content": "'+sample+' annotation excel done"}}\'\n')
    single_cmd.close()
    ofile = open('annotation/annotation.sh', 'a')
    ofile.write('python '+m_path+'/src/annotation_filt.py -sn '+sample+' -c_f_m c\n')
    ofile.close()

def generate_trio_cmd_sentieon(k, filtmode, pedigree_dict, d1, dirname):              
    sample = d1.loc[k]['sample']
    path = 'trio/'+sample
    if not os.path.exists(path):
        os.makedirs(path)                                                                                                       
    trio_cmd = open(path+'/trio_cmd_sentieon.sh', 'w')
    #trio_cmd.write('#!/bin/sh\n\n#SBATCH -J '+sample+'\n#SBATCH -p BP10\n#SBATCH -N 1\n#SBATCH -n 4\n#SBATCH --mail-type=END,FAIL\n#SBATCH --mail-user=1573077420@qq.com\n\nsample='+sample+'\n\n')
    trio_cmd.write('#!/bin/sh\n\n#SBATCH -J '+sample+'\n#SBATCH -p BP10\n#SBATCH -N 1\n#SBATCH -n 4\n#SBATCH --mail-type=FAIL\n#SBATCH --mail-user=1573077420@qq.com\n\nsample='+sample+'\n\n')
    trio_cmd.write('filtmode=\''+filtmode+'\'\npath=`pwd`\nvar=`echo $path | awk -F \'/\' \'{print $NF}\'`\n')
    trio_cmd.write('if [ $var = \'peddy\' ]\nthen\nsh '+m_path+'/src/trio_sentieon.sh $sample $SLURM_NPROCS\n')
    trio_cmd.write('curl \'https://oapi.dingtalk.com/robot/send?access_token=dc5bf30a873950e8c73bb37abf1851ebaf5d83f2792f08684ace2be621a12563\' -H \'Content-Type: application/json\' -d \'{"msgtype": "text","text": {"content": "peddy done"}}\'\n')
    trio_cmd.write('else\nsh '+m_path+'/src/trio_sentieon.sh $sample $SLURM_NPROCS $filtmode\n')
    trio_cmd.write('[ -e ../../$sample/2_mapping/$sample.depth.sample_gene_summary ] && cp ../../$sample/2_mapping/$sample.depth.sample_gene_summary ../../annotation\n')
    trio_cmd.write('cp segtrio/$sample.{ann.hg19_multianno.txt,CADD,maf,link} ../../annotation\n')
    trio_cmd.write('cd ../../annotation\n')
    trio_cmd.write('echo $sample >> list\n')
    trio_cmd.write('python '+m_path+'/src/score_re.py -sn $sample\n')
    relation = []
    for j in pedigree_dict[k]:
        relation.append(d1.loc[j]['relationship'])
    a = ' '.join(relation)
    print a
    annotationfile = open('annotation/annotation.sh', 'a')
    if u'父' in a and u'母' in a:                                                                                               
        trio_cmd.write('python '+m_path+'/src/annotation_filt.py -sn $sample -c_f_m c_f_m\n')
        annotationfile.write('python '+m_path+'/src/annotation_filt.py -sn '+sample+' -c_f_m c_f_m\n')
    else:                                                                                                                       
        if u'父' in a:
            trio_cmd.write('python '+m_path+'/src/annotation_filt.py -sn $sample -c_f_m c_f\n')
            annotationfile.write('python '+m_path+'/src/annotation_filt.py -sn '+sample+' -c_f_m c_f\n')
        else:
            trio_cmd.write('python '+m_path+'/src/annotation_filt.py -sn $sample -c_f_m c_m\n')
            annotationfile.write('python '+m_path+'/src/annotation_filt.py -sn '+sample+' -c_f_m c_m\n')
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


