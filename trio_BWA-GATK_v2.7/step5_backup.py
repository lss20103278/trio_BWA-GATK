#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import pandas as pd

ofile = open('step5_backup.sh', 'w')
dirname = os.getcwd().split('/')[-2:]
filelist = os.listdir('.')
infile = ''
for i in filelist:
	if '分析任务单' in i:
		seperator = '分析任务单'
		break
	if '产品信息' in i:
		seperator = '产品信息'
		break
infile = i

#ofile.write('mkdir -p raw/bk\n')
#
#d1 = pd.read_excel(infile)
#sample_ID = d1[u'样本编号']
#ofile.write('for i in '+' '.join(sample_ID)+'; do mkdir -p raw/bk/$i; cp raw/*$i*gz raw/bk/; awk \'{if(match($2,\"\'$1\'\")){print $0}}\' raw/md5sum.txt > raw/bk/$i/$i.md5; done\n')
strategy = os.getcwd().split('/')[-1].split('_')[-1]
print strategy
truth = raw_input('Is the strategy correct? If not, please enter the correct strategy(WES PANEL Agilent_wes): ')
if truth != '':
	if truth != 'yes':
		strategy = truth
if strategy == 'WES' or strategy == 'wes':
	ofile.write('mv raw/bk/* /anaData/anaData004/children_hos_genetic/rawdata/2018/WES/\n')
elif strategy == 'Agilent_wes':
	ofile.write('mv raw/bk/* /anaData/anaData004/children_hos_genetic/rawdata/2018/安捷伦全外/\n')
else:
	ofile.write('mv raw/bk/* /anaData/anaData004/children_hos_genetic/rawdata/2018/PANEL/\n')
ofile.write('rm -r raw/bk\n')
ofile.write('mkdir -p /DATA/BPshare/opm/遗传病-交大附属儿医/2018/'+'/'.join(dirname)+'/vcf\n')
if os.path.exists('annotation/list_all'):
	ofile.write('for i in `cat annotation/list_all`; do cp trio/$i/sep/*vcf /DATA/BPshare/opm/遗传病-交大附属儿医/2018/'+'/'.join(dirname)+'/vcf; done\n')
else:
	ofile.write('for i in `cat annotation/list`; do cp trio/$i/sep/*vcf /DATA/BPshare/opm/遗传病-交大附属儿医/2018/'+'/'.join(dirname)+'/vcf; done\n')
ofile.write('cp '+infile.split(seperator)[0]+'交付信息表*xlsx'+' /DATA/BPshare/opm/遗传病-交大附属儿医/2018/'+'/'.join(dirname)+'\n')
ofile.write('cp annotation/*_v2.8.xlsx /DATA/BPshare/opm/遗传病-交大附属儿医/2018/'+'/'.join(dirname)+'\n')
ofile.write('cp -r /DATA/BPshare/opm/遗传病-交大附属儿医/2018/'+'/'.join(dirname)+' .\n')
ofile.write('zip -r '+dirname[-1]+'.zip '+dirname[-1]+'\n')
ofile.write('rm -r '+dirname[-1]+'\n')
ofile.write('mkdir bam\n')
ofile.write('cp list bam/list_'+dirname[-1]+'\n')
ofile.write('for i in `cat list`; do mv $i/2_mapping/$i.sort.mkdup.bam* bam; done\n')
ofile.write('for i in `cat list`; do rm $i/2_mapping/$i.sort.bam*; done\n')
ofile.write('scp bam/* klyang@122.112.248.194:/media/bpdata/data/exon_temp2/\n')
ofile.close()
with open('step5_backup.sh') as f:
	for l in f:
		print l.strip('\n')
truth = raw_input('Is the shell script right? If not, please enter no: ')
if truth != '':
	sys.exit()
else:
	ofile.system('bash backup.sh')
