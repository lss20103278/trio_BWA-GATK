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

sys.path.append("/DATA/sslyu/trio_BWA-GATK_3.0/")
from lib.prepare import *
from lib.GAF import *
from lib.xlsx import *

kwargs_raw = sys.argv[1:]
kwargs={'-sn':'', '-l':'', '-strategy':'', '-filtmode':'vqsr'}
for i in range(len(kwargs_raw)):
	for j in kwargs.keys():
		if(kwargs_raw[i]==j):
			kwargs.update({j:kwargs_raw[i+1]})

#read sample list
if(kwargs.get('-sn')!=''):
    list=[]
    list.append(kwargs.get('-sn'))
    print "\nNote:sigle sample mode"
else:
    if kwargs.get('-l')!='':
        list=pd.read_table(kwargs.get('-l'),header=None)[0].tolist()
        print "\nNote:list samples mode"
strategy=kwargs.get('-strategy')
filtmode=kwargs.get('-filtmode')

if len(list) == 0:
    print """
    Examples: 
    python /DATA/sslyu/trio_BWA-GATK_3.0/src/step1_prepare.py -sn excelfile -strategy WES -filtmode vqsr
    python /DATA/sslyu/trio_BWA-GATK_3.0/src/step1_prepare.py -l \'list of names of excels\' -strategy IDT-PANEL -filtmode hard 
    """
    sys.exit()
else:    
    for sn in list:
        generate_dbevn(strategy)
        #if os.path.exists('dbevn.sh'):
        #    with open('dbevn.sh') as f:
        #        for l in f:
        #            if l.startswith('panel'):
        #                print l.strip('\n')
        #    truth = raw_input('Is dbevn.sh correct?')
        #    if truth == 'no':
        #        strategy = raw_input('strategy(WES PANEL Agilent-PANEL IDT-PANEL Agilent-wes brain blood-disease, if not exist, please enter none):')
        #        generate_dbevn(strategy)
        #    else:
        #        strategy = os.getcwd().split('/')[-1].split('_')[-1]
        #else:
        #    strategy = raw_input('strategy(WES PANEL Agilent-PANEL IDT-PANEL Agilent-wes brain blood-disease, if not exist, please enter none):')
        #    generate_dbevn(strategy)
        #    with open('dbevn.sh') as f:
        #        for l in f:
        #            if l.startswith('panel'):
        #                print l.strip('\n')
        #if strategy == 'IDT-PANEL':
        #    if not os.path.exists('CNV'):
        #        os.system('mkdir CNV')

        #filtmode = raw_input('-filtmode:(hard,vqsr,...) ')
        print '-filtmode: '+filtmode

        for i in ['ped', 'peddy', 'gvcf', 'annotation', 'raw']:
            if not os.path.exists(i):
                os.mkdir(i)
        excel = pd.read_excel(sn)
        #excel = append_sample_excel(excel)
        excel = append_sample_txt(excel)
        print excel[[u'姓名', u'性别', 'sample']]
        d1 = append_gender(excel)
        d1 = append_relation(d1)
        d1 = append_pedigree(d1)
        d1 = append_phenotype(d1)
        d1 = append_father_mother(d1)
        print d1

        #### step2: calculate Reads(M) bases(Mb) Q30
        #manually rectify the colnames of reads, bases and Q30 to 'Reads','bases(Mb)','Q30'
        #print d1.columns
        #wrong_header = raw_input('Please enter the colnames of reads, bases and Q30: ')
        #wrong_header = wrong_header.split(', ')
        #print d1[wrong_header]
        #right_header = ['Reads','bases(Mb)','Q30']
        #for i in range(len(wrong_header)):
        #    d1.rename(columns={wrong_header[i]:right_header[i]}, inplace=True)
        ##print unicode(' '.join(d1.columns.tolist()), "utf-8")

        #header_truth1 = raw_input('Do Reads bases(Mb) Q30 need to be calculated?')
        #if header_truth1 != 'no':
        #    sys.exit()
        #d1.to_csv(sn+'.tmp', index=None, sep="\t", encoding='utf-8') #Maybe the encoding should be considered?

        ##### step3: generate cmd files
        pedigree_dict = generate_pedigree(d1)
        trio = generate_trio(pedigree_dict)
        
        for i in d1['sample']:
        	if not os.path.exists(str(i)):
        		os.mkdir(str(i))
        truth = "yes"
        if truth != 'no':
            n = 0
            with open('/DATA/sslyu/trio_BWA-GATK_3.0/src/run2_sentieon.sh') as f:
            	for l in f:
            		n = n+1
            		if n>22 and n<=30:
            			print l
    
            single = generate_single(d1)
            d1.index = d1[u'姓名']
            if len(single) > 0:
                single_ID = [d1.loc[i]['sample'] for i in single]
                list_single = open('list_single', 'w')
                for i in single_ID:
                    list_single.write(str(i)+'\n')
                list_single.close()
                if strategy == "WES":
                    cmd_single = "clean mapping precalling gvcf index_single gtgvcf vqsr left-normalize depth qcsum qcvcf annotation"
                else:
                    cmd_single = "clean mapping precalling gvcf index_single gtgvcf hardfilt left-normalize depth qcsum qcvcf annotation"
                for i in single_ID:
                    generate_single_cmd(str(i))
                    generate_sinle_ped(str(i),d1)
                    generate_run2_sentieon(str(i),cmd_single)
            if len(trio) > 0:
                if not os.path.exists('trio'):
                    os.mkdir('trio')
                generate_trio_list(trio, d1)
                cmd_trio = "clean mapping precalling gvcf index_trio depth qcsum"
                list_trio = open('list_trio', 'w')
                for k in pedigree_dict:
                    generate_ped(k,d1)
                    ID = d1.loc[k]['sample']
                    path = 'trio/'+ID
                    if not os.path.exists(path):
                        os.makedirs(path)
                    #os.system('cp /DATA/sslyu/trio_BWA-GATK_3.0/src/{mendel_to_annovar.py,sort_sample.py,vep-mendelscan.sh,vfilt.sh} '+path)
                    trio_cmd = generate_trio_cmd(path, filtmode)
                    trio_cmd.write('[ -e ../../$sample/2_mapping/$sample.depth.sample_gene_summary ] && cp ../../$sample/2_mapping/$sample.depth.sample_gene_summary ../../annotation\n')
                    trio_cmd.write('cp segtrio/$sample.{ann.hg19_multianno.txt,CADD,maf,link} ../../annotation\n')
                    if not os.path.exists('annotation'):
                        os.makedirs('annotation')
                    trio_cmd.write('cd ../../annotation\n')
                    trio_cmd.write('echo $sample >> list\n')
                    trio_cmd.write('python /DATA/sslyu/trio_BWA-GATK_3.0/src/score_re.py -sn $sample\n')
                    relation = []
                    for j in pedigree_dict[k]:
                        relation.append(d1.loc[j]['relationship'])
                    a = ' '.join(relation)
                    annotationfile = open('annotation/annotation.sh', 'a')
                    if u'父' in a and u'母' in a:
                        trio_cmd.write('python /DATA/sslyu/trio_BWA-GATK_3.0/src/annotation_filt_ver3.0.py -sn $sample -c_f_m c_f_m\n')
                        annotationfile.write('python /DATA/sslyu/trio_BWA-GATK_3.0/src/annotation_filt_ver3.0.py -sn '+ID+' -c_f_m c_f_m\n')
                    else:
                        if u'父' in a:
                            trio_cmd.write('python /DATA/sslyu/trio_BWA-GATK_3.0/src/annotation_filt_ver3.0.py -sn $sample -c_f_m c_f\n')
                            annotationfile.write('python /DATA/sslyu/trio_BWA-GATK_3.0/src/annotation_filt_ver3.0.py -sn '+ID+' -c_f_m c_f\n')
                        else:
                            trio_cmd.write('python /DATA/sslyu/trio_BWA-GATK_3.0/src/annotation_filt_ver3.0.py -sn $sample -c_f_m c_m\n')
                            annotationfile.write('python /DATA/sslyu/trio_BWA-GATK_3.0/src/annotation_filt_ver3.0.py -sn '+ID+' -c_f_m c_m\n')
                    trio_cmd.write('echo "end `date`" >> ../trio/$sample/finished\n')
                    trio_cmd.write('fi\n')
                    trio_cmd.close()
                    annotationfile.close()

                    trio_k_ID = []
                    for i in pedigree_dict[k]:
                        trio_k_ID.append(d1.loc[i]['sample'])
                    for i in trio_k_ID:
                        list_trio.write(str(i)+'\n')
                        generate_cmd(str(i))
                        generate_run2_sentieon(str(i),cmd_trio)
                list_trio.close()
##############################need to wait the last sub_sep.sh to complete############################                       
######## step4:trio_ped_analysis                    
#                generate_ped(k,d1)
#                generate_peddy(trio,d1)
#                k_ID = d1.loc[k]['sample']
#                mkdir -p trio/k_ID
#                for i in trio_K_ID:
#                    os.system('cp '+i+'/3_variants/gvcf/* trio/gvcf')
#                os.system('sh trio/sub_sep.sh '+k) 
#                os.system('sh trio/sub_sep.sh peddy')
##################################################################################################

        ###### step5: generate related files and dirs
        if u'样本编号' in d1.columns and u'原始样本ID' in d1.columns:
            d1_tmp = d1[[u'样本编号', u'原始样本ID']]
        if u'样本ID' in d1.columns and u'原始编码' in d1.columns:
            d1_tmp = d1[[u'样本ID', u'原始编码']]
        d1_tmp.to_csv('rename'+str(list.index(sn)), index=None, header=None, sep="\t")
        d1['sample'].to_csv('list'+str(list.index(sn)), header=None, index=None)
        if os.path.exists('CNV'):
            d1['CNV'] = '/BP12_share/sslyu/bam/'+d1['sample']+'.sort.mkdup.bam'
            d1['CNV'].to_csv('CNV/list'+str(list.index(sn)), header=None, index=None)
            os.system('cat CNV/list'+str(list.index(sn))+' >> CNV/list')
            os.system('rm CNV/list'+str(list.index(sn)))
        os.system('cat list'+str(list.index(sn))+' >> list')
        os.system('rm list'+str(list.index(sn)))

os.system('sort list |uniq > tmp; mv tmp list')                
if os.path.exists('CNV'):
    os.system('cp list CNV/id')
    os.system('sort CNV/list |uniq > tmp; mv tmp CNV/list')
if os.path.exists('trio'):
    os.system('sort trio/list |uniq > tmp; mv tmp trio/list')
if not os.path.exists('peddy'):
    os.mkdir('peddy')
os.system('cp /DATA/sslyu/trio_BWA-GATK_3.0/src/trio_cmd.sh peddy')
os.system('echo peddy > peddy/list')
if os.path.exists('list_single'):
    os.system('for i in `cat list_single`; do cat ped/$i.mendel.ped >> peddy_tmp; done')
if os.path.exists('trio/list'):
    os.system('for i in `cat trio/list`; do cat ped/$i.mendel.ped >> peddy_tmp; done')
os.system('sort peddy_tmp |uniq |awk \'{OFS="\t"; print \'0\',$2,$3,$4,$5,$6}\' > ped/peddy.ped')
with open('ped/peddy.ped') as f:
    for l in f:
        print l.strip('\n')
os.system('sort annotation/annotation.sh |uniq > annotation/tmp; mv annotation/tmp annotation/annotation.sh')        
report_dirname = os.getcwd().split('/')[-1]
if not os.path.exists(report_dirname):
    os.system('mkdir '+report_dirname)
#####################################################################################
#import shutil
#shutil.copyfile('/DATA/sslyu/trio_BWA-GATK_2.7/sub.sh', 'sub.sh')

#        ####### step1: copy raw date    
#        ofile = open('step1_prepare.sh', 'w')
#        cp_judgement = raw_input('Do you need to copy the raw data? If not, please enter no:')
#        if cp_judgement != 'no':
#            if not os.path.exists('raw'):
#                os.system('mkdir raw')
#            ID = d1['sample']
#            ofile.write('mkdir raw\n')
#            
#            redo_judement = raw_input('Is this a redo work? If not, please enter no:  ')
#            if redo_judement == 'no':
#                single_judgement = raw_input('Is the raw data in a single dir or in a set of subdirs? If single, please enter single, if set, please enter set:')
#                if single_judgement == 'single':
#                    rawdata_dir = raw_input('Please enter the absolute path of the raw data:  ')
#                    if not rawdata_dir.startswith('/DATA'):
#                        rawdata_dir = os.getcwd()+'/'+rawdata_dir
#                    os.system('ls '+rawdata_dir)
#                    origin = raw_input('please enter the pattern of R1 and R2:  ')
#                    ofile.write('a=('+' '.join(d1[u'样本编号'])+')\n')
#                    ofile.write('c=('+' '.join(origin.split(' '))+')\n')
#                    n = str(len(d1[u'样本编号']))
#                    md5 = raw_input('please enter the absolute path of the md5 file:  ')
#                    if not md5.startswith('/DATA'):
#                        md5 = os.getcwd()+'/'+md5
#                    ofile.write('cp '+md5+' raw/origin_md5\n')
#                    d1.index = d1['sample']
#                    for i in d1['sample']:
#                        if not os.path.exists(str(i)):
#                            os.mkdir(str(i))
#                        cp_rawdata_file = open(str(i)+'/cp_rawdata.sh','w')
#                        cp_rawdata_file.write('#!/bin/sh\n\n#SBATCH -J '+i+'\n#SBATCH -p BP10\n#SBATCH -N 1\n#SBATCH -n 4\n\nsample='+i+'\n\nmkdir -p ../raw/bk/'+d1.loc[i][u'样本编号']+'\nawk \'{if(match($2,\"\'$list\'\")){print $0}}\' '+md5+' > ../raw/bk/'+d1.loc[i][u'样本编号']+'/'+d1.loc[i][u'样本编号']+'.md5\ncp '+rawdata_dir+'/*'+d1.loc[i][u'样本编号']+'*'+origin.split(' ')[0]+'* ../raw/bk/'+d1.loc[i][u'样本编号']+'\ncp '+rawdata_dir+'/*'+d1.loc[i][u'样本编号']+'*'+origin.split(' ')[1]+'* ../raw/bk/'+d1.loc[i][u'样本编号']+'\ncp '+rawdata_dir+'/*'+d1[u'样本编号'][i]+'*'+origin.split(' ')[0]+'* ../raw/'+i+'_R1.fastq.gz\ncp '+rawdata_dir+'/*'+d1[u'样本编号'][i]+'*'+origin.split(' ')[1]+'* ../raw/'+i+'_R2.fastq.gz\n')
#                        cp_rawdata_file.close()
#        
#                    ofile.write('for i in `seq '+n+'`; do mkdir -p raw/bk/${a[$i-1]}; awk \'{if(match($2,\"\'${a[$i-1]}\'\")){print $0}}\' '+md5+' > raw/bk/${a[$i-1]}/${a[$i-1]}.md5; for j in `seq 2`; do cp '+rawdata_dir+'/*${a[$i-1]}*${c[$j-1]}* raw/bk/${a[$i-1]}; done; done\n')
#                    for i in range(len(d1[u'样本编号'])):
#                    	for j in range(2):
#                    		ofile.write('cp '+rawdata_dir+'/*'+d1[u'样本编号'][i]+'*'+origin.split(' ')[j]+'* raw/'+ID[i]+'_R'+str(j+1)+'.fastq.gz;')
#                    	ofile.write('\n')
#                else:
#                    parent_raw_dir = raw_input('please enter the parent directory:  ')
#                    for i in d1[u'样本编号']:
#                        ##rawdata_dir = '/anaData/anaData004/children_hos_genetic/rawdata/2018/'+strategy+'/'
#                        #raw_dir = parent_raw_dir+strategy+'/'
#                        rawdata_dir = parent_raw_dir+'/*'+i+'*'
#                        os.system('ls '+rawdata_dir)
#                    origin = raw_input('please enter the pattern of R1 and R2:  ')
#                    md5 = raw_input('please enter the pattern of md5:  ')
#                    d1.index = d1['sample']
#                    for i in d1['sample']:
#                        if not os.path.exists(str(i)):
#                            os.mkdir(str(i))
#                        cp_rawdata_file = open(str(i)+'/cp_rawdata.sh','w')
#                        #cp_rawdata_file.write('#!/bin/sh\n\n#SBATCH -J '+i+'\n#SBATCH -p BP10\n#SBATCH -N 1\n#SBATCH -n 4\n\nsample='+i+'\n\nmkdir -p ../raw/bk/'+d1.loc[i][u'样本编号']+'\ncp -r '+parent_raw_dir+'/*'+d1.loc[i][u'样本编号']+'*/* ../raw/bk/'+d1.loc[i][u'样本编号']+'\ncp '+parent_raw_dir+'/*'+d1.loc[i][u'样本编号']+'*/*'+d1.loc[i][u'样本编号']+'*'+origin.split(' ')[0]+'* ../raw/'+i+'_R1.fastq.gz\ncp '+parent_raw_dir+'/*'+d1.loc[i][u'样本编号']+'*/*'+d1.loc[i][u'样本编号']+'*'+origin.split(' ')[1]+'* ../raw/'+i+'_R2.fastq.gz\n')
#                        cp_rawdata_file.write('#!/bin/sh\n\n#SBATCH -J '+i+'\n#SBATCH -p BP10\n#SBATCH -N 1\n#SBATCH -n 4\n\nsample='+i+'\n\nmkdir -p ../raw/bk/'+d1.loc[i][u'样本编号']+'\ncp '+parent_raw_dir+'/*'+d1.loc[i][u'样本编号']+'*/*'+md5+'* ../raw/bk/'+d1.loc[i][u'样本编号']+'/'+d1.loc[i][u'样本编号']+'.md5\ncp '+parent_raw_dir+'/*'+d1.loc[i][u'样本编号']+'*/*gz ../raw/bk/'+d1.loc[i][u'样本编号']+'/\ncp '+parent_raw_dir+'/*'+d1.loc[i][u'样本编号']+'*/*'+d1.loc[i][u'样本编号']+'*'+origin.split(' ')[0]+'* ../raw/'+i+'_R1.fastq.gz\ncp '+parent_raw_dir+'/*'+d1.loc[i][u'样本编号']+'*/*'+d1.loc[i][u'样本编号']+'*'+origin.split(' ')[1]+'* ../raw/'+i+'_R2.fastq.gz\n')
#                        cp_rawdata_file.close()
#                	for i in range(len(d1[u'样本编号'])):
#                		for j in range(2):
#                			ofile.write('cp '+parent_raw_dir+'/*'+d1[u'样本编号'][i]+'*/*'+d1[u'样本编号'][i]+'*'+origin.split(' ')[j]+'* raw/'+ID[i]+'_R'+str(j+1)+'.fastq.gz;')
#                		ofile.write('\n')
#            else:
#                parent_raw_dir = raw_input('please enter the parent directory:  ')
#                for i in d1[u'样本编号']:
#                    rawdata_dir = parent_raw_dir+'/'+i
#                    os.system('ls '+rawdata_dir)
#                origin = raw_input('please enter the pattern of R1 and R2:  ')
#                d1.index = d1['sample']
#                for i in d1['sample']:
#                    if not os.path.exists(str(i)):
#                        os.mkdir(str(i))
#                    cp_rawdata_file = open(str(i)+'/cp_rawdata.sh','w')
#                    cp_rawdata_file.write('#!/bin/sh\n\n#SBATCH -J '+i+'\n#SBATCH -p BP10\n#SBATCH -N 1\n#SBATCH -n 4\n\nsample='+i+'\n\ncp '+parent_raw_dir+'/'+d1.loc[i][u'样本编号']+'/*'+d1.loc[i][u'样本编号']+'*'+origin.split(' ')[0]+'* ../raw/'+i+'_R1.fastq.gz\ncp '+parent_raw_dir+'/'+d1.loc[i][u'样本编号']+'/*'+d1.loc[i][u'样本编号']+'*'+origin.split(' ')[1]+'* ../raw/'+i+'_R2.fastq.gz\n')
#                    cp_rawdata_file.close()
#                for i in range(len(d1[u'样本编号'])):
#                    for j in range(2):
#                        ofile.write('cp '+parent_raw_dir+'/'+d1[u'样本编号'][i]+'/*'+d1[u'样本编号'][i]+'*'+origin.split(' ')[j]+'* raw/'+ID[i]+'_R'+str(j+1)+'.fastq.gz;')
#                    ofile.write('\n')
#        ofile.close()

#if os.path.getsize('step1_prepare.sh') == 0:
#    os.system('rm step1_prepare.sh')
#if os.path.exists('step1_prepare.sh'):    
#    with open('step1_prepare.sh') as f:
#    	for l in f:
#    		print l.strip('\n')
