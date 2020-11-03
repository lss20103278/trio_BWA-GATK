#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
version:2.6
by:lss
"""

import numpy as np
import pandas as pd
import sys
reload(sys)
sys.setdefaultencoding("utf8")
import os

if '-redo' not in sys.argv:
	print 'miss -redo ...'
	sys.exit()


def generate_cmd(path):
	os.system('cp /DATA/sslyu/trio_BWA-GATK_2.7/cmd.sh '+path)
	
def generate_sub():
	os.system('cp /DATA/sslyu/trio_BWA-GATK_2.7/sub.sh .')
generate_sub()

strategy = raw_input('strategy(WES PANEL Agilent IDT Agilent_wes brain blood_disease, if not exist, please enter none):')
def generate_dbevn():
	dbevn = open('dbevn.sh', 'w')
	if strategy != 'none':
		with open('/DATA/sslyu/trio_BWA-GATK_2.7/dbevn.sh') as f:
			for l in f:
				if l.startswith('#panel'):
					if strategy == 'Agilent':
						if strategy in l and 'Agilent_wes' not in l:
							dbevn.write(l[1:])
						else:
							dbevn.write(l)
					else:
						if strategy in l:
							dbevn.write(l[1:])
						else:
							dbevn.write(l)
				else:
					dbevn.write(l)
	else:
		with open('/DATA/sslyu/trio_BWA-GATK_2.7/dbevn.sh') as f:
			for l in f:
				dbevn.write(l)
		panel = raw_input('please enter the absolute path of the bed:')
		dbevn.write('panel=\"'+panel+'\"\n')
	dbevn.close()
generate_dbevn()
def generate_run2_sentieon(path,str):
	run2_sentieon = open(path+'/run2_sentieon.sh', 'w')
	n = 0
	with open('/DATA/sslyu/trio_BWA-GATK_2.7/run2_sentieon.sh') as f:
		for l in f:
			n = n+1
			if n>0 and n<= 23:
				run2_sentieon.write(l)
			elif n>= 30:
				run2_sentieon.write(l)
			elif n==24:
				run2_sentieon.write('cmd="'+str+'"\n')
			else:
				continue
	run2_sentieon.close()
import shutil
shutil.copyfile('/DATA/sslyu/trio_BWA-GATK_2.7/sub.sh', 'sub.sh')

kwargs_raw = sys.argv[1:]
kwargs={'-excel':'', '-md5':'', '-redo':''}
for i in range(len(kwargs_raw)):
	for j in kwargs.keys():
		if(kwargs_raw[i]==j):
			kwargs.update({j:kwargs_raw[i+1]})

#read the excel
if(kwargs.get('-excel')!=''):
	infile = kwargs.get('-excel')
else:
	filelist = os.listdir('.')
	infile = ''
	for i in filelist:
		if '分析任务单' in i:
			infile = i
			break
		if '产品信息表' in i:
			infile = i
			break
d1 = pd.read_excel(infile)

if(kwargs.get('-md5')!=''):
	md5 = kwargs.get('-md5')

if(kwargs.get('-redo')!=''):
	redo = kwargs.get('-redo')
	if redo == 'no':
		if '-md5' not in sys.argv:
			print 'miss -md5 ....'
			sys.exit()
#this block has a bug
#ID = []
#for i in range(d1.shape[0]):
#	if np.isnan(d1.iloc[i][u'原始样本ID']):
#		ID.append(d1.iloc[i][u'样本编号'])
#	else:
#		ID.append(d1.iloc[i][u'原始样本ID'])
#d1['sample'] = ID

ID = []
for i in range(len(d1[u'原始样本ID'].isnull())):
	if d1[u'原始样本ID'].isnull()[i]:
		ID.append(d1.iloc[i][u'样本编号'])
	else:
		ID.append(d1.iloc[i][u'原始样本ID'])
ofile = open('step1_prepare.sh', 'w')
ofile.write('mkdir raw\n')
ofile.write('a=('+' '.join(d1[u'样本编号'])+')\n')
#ofile.write('b=('+' '.join(ID)+')\n')
if md5:
	rawdata_dir='/'.join(md5.split('/')[:-1])
	print rawdata_dir
	truth = raw_input('Is the absolute path of the raw data correct? If not, please enter the absolute path of the raw data:')
	if truth != '':
		rawdata_dir=truth
else:
	rawdata_dir=raw_input('please enter the absolute path of the raw data:')
print 'a=('+' '.join(d1[u'样本编号'])+')'
#print 'b=('+' '.join(ID)+')'
print os.listdir(rawdata_dir)
origin = raw_input('please enter the pattern of R1 and R2:')
#ofile.write('c=('+origin+')\n')
#ofile.write('d=(R1 R2)\n')
n = str(len(ID))
if redo == 'yes':
#	os.system('sed -i \'s/^/#&/\' step1_prepare.sh')
	for i in range(len(d1[u'样本编号'])):
		for j in range(2):
			ofile.write('cp '+rawdata_dir+'/'+d1[u'样本编号'][i]+'/*'+origin.split(' ')[j]+'* raw/'+ID[i]+'_R'+str(j+1)+'.fastq.gz;')
		ofile.write('\n')
else:
	ofile.write('for i in `seq '+n+'`; do mkidr -p raw/bk/${a[$i-1]}; awk \'{if(match($2,\"\'${a[$i-1]}\'\")){print $0}}\' '+md5+' > raw/bk/${a[$i-1]}/${a[$i-1]}.md5; for j in `seq 2`; do cp '+rawdata_dir+'/*${a[$i-1]}*${c[$j-1]}* raw/bk/${a[$i-1]}; done; done\n')
	for i in range(len(d1[u'样本编号'])):
		for j in range(2):
			ofile.write('cp '+rawdata_dir+'/*'+d1[u'样本编号'][i]+'*'+origin.split(' ')[j]+'* raw/'+ID[i]+'_R'+str(j+1)+'.fastq.gz;')
		ofile.write('\n')
#	os.system('sed -i \'s/^/#&/\' step1_prepare.sh')
############################################################################################################	
#	for i in range(len(d1[u'样本编号'])):
#		for file in os.listdir(rawdata_dir+str(d1[u'样本编号'][i])):
#			if '_R1.fastq.gz' not in file:
#				print 'the name of fastq file is not appropriate'
#				print os.listdir(rawdata_dir+str(d1[u'样本编号'][i]))
#				cmd1 = raw_input('please enter the name of the R1 fastq:')
#				cmd2 = raw_input('please enter the name of the R2 fastq:')
#				ofile.write('cp '+rawdata_dir+'$a{i}/'+cmd1+' raw/'+'$b{i}_R1.fastq.gz\n')
#				ofile.write('cp '+rawdata_dir+'$a{i}/'+cmd2+' raw/'+'$b{i}_R2.fastq.gz\n')
#				break
#			else:
#				for j in ['_R1.fastq.gz', '_R2.fastq.gz']:
#					if j in file:
#						ofile.write('cp '+rawdata_dir+'$a{i}/'+file+' raw/'+'$b{i}_'+j+'\n')
#	ofile.close()
#else:
#	rawdata_dir=raw_input('please enter the absolute path of the raw data:')
#	
#	flag = 'on'
#	if not os.path.exists('raw/bk'):
#		flag = 'off'
#	if flag == 'off':
#		ofile = open('prepare.sh', 'a')
#		ofile.write('for i in '+' '.join(d1[u'样本编号'])+'; do mkdir -p raw/bk/$i; cp '+rawdata_dir+'/*$i*gz raw/bk/$i; awk \'{if(match($2,\"\'$i\'\")){print $0}}\' '+md5+' > raw/bk/$i/$i.md5; done\n')
#		ofile.close()
#		os.system('bash prepare.sh')
#		os.system('awk \'{print \"#\",$0}\' prepare.sh') #can't execute, don't know why?
#		ofile = open('prepare.sh', 'a')
#		raw_l = os.listdir('raw/bk')
#		for i in range(len(raw_l)):
#			if d1[u'原始样本ID'].isnull()[i]:
#				for file in os.listdir('raw/bk/'+raw_l[i]):
#					if '_R1.fastq.gz' not in file:
#						print 'the name of fastq file is not appropriate'
#						print os.listdir('raw/bk/'+raw_l[i])
#						cmd1 = raw_input('please enter the name of the R1 fastq:')
#						cmd2 = raw_input('please enter the name of the R2 fastq:')
#						ofile.write('cp raw/bk/'+raw_l[i]+'/'+cmd1+' raw/'+raw_l[i]+'_R1.fastq.gz\n')
#						ofile.write('cp raw/bk/'+raw_l[i]+'/'+cmd2+' raw/'+raw_l[i]+'_R2.fastq.gz\n')
#						break
#					else:
#						for j in ['_R1.fastq.gz', '_R2.fastq.gz']:
#							if j in file:
#								ofile.write('cp raw/bk/'+raw_l[i]+'/'+file+' raw/'+raw_l[i]+j+'\n')
#			else:
#				for file in os.listdir('raw/bk/'+raw_l[i]):
#					if '_R1.fastq.gz' not in file:
#						print 'the name of fastq file is not appropriate'
#						print os.listdir('raw/bk/'+raw_l[i])
#						cmd1 = raw_input('please enter the name of the R1 fastq:')
#						cmd2 = raw_input('please enter the name of the R2 fastq:')
#						ofile.write('cp raw/bk/'+raw_l[i]+'/'+cmd1+' raw/'+str(d1.iloc[i][u'原始样本ID'])+'_R1.fastq.gz\n')
#						ofile.write('cp raw/bk/'+raw_l[i]+'/'+cmd2+' raw/'+str(d1.iloc[i][u'原始样本ID'])+'_R2.fastq.gz\n')
#						break
#					else:
#						for j in ['_R1.fastq.gz', '_R2.fastq.gz']:
#							if j in file:
#								ofile.write('cp raw/bk/'+raw_l[i]+'/'+file+' raw/'+str(d1.iloc[i][u'原始样本ID'])+j+'\n')
#		ofile.close()
#	os.system('bash prepare.sh')
#	os.system('awk \'{print \"#\",$0}\' prepare.sh')
	
#	for i in 
#			if d1[u'原始样本ID'].isnull()[i]:
#				if d1.iloc[i][u'样本编号'] in file:
#					if 'R1.fastq' not in file:
#						if '_1.fq' in file:
#							os.rename('raw/'+file,'raw/'+str(d1.iloc[i][u'样本编号'])+'_R1.fastq.gz')
#						if '_2.fq' in file:
#							os.rename('raw/'+file,'raw/'+str(d1.iloc[i][u'样本编号'])+'_R2.fastq.gz')
#					else:
#						os.rename('raw/'+file,'raw/'+str(d1.iloc[i][u'样本编号'])+'_'+'_'.join(file.split('_')[1:]))
#			else:
#				if d1.iloc[i][u'样本编号'] in file:
#					if '
#					os.rename('raw/'+file,'raw/'+str(d1.iloc[i][u'原始样本ID'])+'_'+'_'.join(file.split('_')[1:]))
###############################################################################################################					
d1['sample'] = ID
d1.index = d1['sample']
d1['gender'] = d1[u'性别'].apply(lambda x:'1' if u'男' in x else '2')
pedigree = {}
for i in d1[u'姓名']:
	pedigree[i] = []
for k in pedigree:
	for i in d1[u'姓名']:
		if k in i:
			pedigree[k].append(i)
trio_m = []
trio = []
trio_p = {}
for k in pedigree:
	if len(pedigree[k]) > 1:
		trio.append(k)
		trio_m.extend(pedigree[k])
		trio_p[k] = pedigree[k]
single = [i for i in d1[u'姓名'] if i not in trio_m]
d1.index = d1[u'姓名']
single_ID = [d1.loc[i]['sample'] for i in single]
trio_m_ID = [d1.loc[i]['sample'] for i in trio_m]

if len(trio) > 0:
	ofile.write('cp -r /DATA/sslyu/trio_BWA-GATK_2.7/trio/ .\n')
ofile.write('mkdir annotation\n')
ofile.write('ln -s /DATA/sslyu/trio_BWA-GATK_2.7/*py .\n')
ofile.close()

child = {}
for i in trio:
	child[i] = d1.loc[i]['sample']

for i in d1['sample']:
	if not os.path.exists(str(i)):
		os.mkdir(str(i))
n = 0
with open('/DATA/sslyu/trio_BWA-GATK_2.7/run2_sentieon.sh') as f:
	for l in f:
		n = n+1
		if n>22 and n<=30:
			print l
if len(single) > 0 and len(trio) > 0:
	cmd_single = raw_input('select cmds from following for single (qc clean mapping sorting precalling calling gvcf gtgvcf vqsr hardfilt annotation depth abase qcsum qcvcf): ')
	cmd_trio = raw_input('select cmds from following for trio (qc clean mapping sorting precalling calling gvcf gtgvcf vqsr hardfilt annotation depth abase qcsum qcvcf): ')
	for i in single_ID:
		generate_cmd(str(i))
		generate_run2_sentieon(str(i),cmd_single)
	for i in trio_m_ID:
		generate_cmd(str(i))
		generate_run2_sentieon(str(i),cmd_trio)
elif len(single) == 0:
	cmd_trio = raw_input('select cmds from following for trio (qc clean mapping sorting precalling calling gvcf gtgvcf vqsr hardfilt annotation depth abase qcsum qcvcf): ')
	for i in d1['sample']:
		generate_cmd(str(i))
		generate_run2_sentieon(str(i),cmd_trio)
else:
	cmd_single = raw_input('select cmds from following for single (qc clean mapping sorting precalling calling gvcf gtgvcf vqsr hardfilt annotation depth abase qcsum qcvcf): ')
	for i in d1['sample']:
		generate_cmd(str(i))
		generate_run2_sentieon(str(i),cmd_single)
d1['sample'].to_csv('list', header=None, index=None)

with open('step1_prepare.sh') as f:
	for l in f:
		print l.strip('\n')
truth = raw_input('Is the shell script correct? If not, please enter no: ')
if truth != '':
	sys.exit()
else:
	os.system('bash step1_prepare.sh')
	os.system('sh sub.sh')

