#!/usr/bin/env python
# -*- coding: utf-8 -*-
#########################################################################
# File Name: /DATA/sslyu/trio_BWA-GATK_2.7/annotation.py
# Author: sslyu
# mail: 1573077420@qq.com
# Created Time:Tue 07 Aug 2018 03:12:48 PM CST
#########################################################################

import numpy as np
import pandas as pd
import sys
reload(sys)
sys.setdefaultencoding("utf8")
import os
import re

if len(sys.argv) < 2:
	print 'miss the value of thisminedge...'
	sys.exit()


filelist = os.listdir('.')
infile = []
for i in filelist:
	if '分析任务单' in i:
		infile.append(i)
	if '产品信息表' in i:
		infile.append(i)

if len(infile) == 0:
	print 'please offer an information excel, the name of the excel should include\"分析任务单\", the columns should be 样本编号 原始样本ID 姓名 性别'
	sys.exit()
	
for i in infile:	
  d1 = pd.read_excel(i, dtype=str)
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
  d1['sample'] = ID
  d1.index = d1['sample']
  d1['gender'] = d1[u'性别'].apply(lambda x:'1' if u'男' in x else '2')
  if 'relation' in d1.columns:
  	pedigree = {}
  	pedigree_ID = []
  	for i in d1['pedigree']:
  		if i not in pedigree_ID:
  			pedigree_ID.append(i)
  	for i in pedigree_ID:
  		pedigree[i] = {}
  	for i in d1['sample']:
  		try:
  			if np.isnan(d1['relationship']):
  				pedigree[d1.loc[i]['pedigree']].update({d1.loc[i]['姓名']:i})
  		except:
  			pedigree[d1.loc[i]['pedigree']].update({d1.loc[i]['relationship']:i})
  	
  	for i in range(len(d1['relation'].isnull())):
  		if d1['relation'].isnull()[i]:
  			
  			for k in pedigree:
  				if d1.
  		for k in pedigree:
  		if d1.loc[i]['pedigree'] == k:
  			if d1.loc[i]['relation']
  			pedigree[k].update({d1.loc[i]['relation']:i})
  		
  	
  	
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
  trio_ID = [d1.loc[i]['sample'] for i in trio]
  trio_m_ID = [d1.loc[i]['sample'] for i in trio_m]
  
  child = {}
  ped_ID = {}
  for i in trio:
  	child[i] = d1.loc[i]['sample']
  	for j in d1[u'姓名']:
  		if i in j:
  			ped_ID.setdefault(i,[]).append(j)
  
  if not os.path.exists('annotation'):
  	os.makedirs('annotation')
  def generate_annotation():
  	if len(trio) > 0 and len(single_ID) > 0:
  		ofile1 = open('annotation/list', 'a')
  		ofile2 = open('annotation/list_n_c_f_m', 'a')
  		for k in ped_ID:
  			a = ' '.join(ped_ID[k])
  			if k+u'父' in a and k+u'母' in a:
  				ofile1.write(d1.loc[k]['sample']+'\n')
  			else:
  				ofile2.write(d1.loc[k]['sample']+'\n')
  		ofile1.close()
  		ofile2.close()
  		ofile = open('annotation/list', 'a')
  		for i in single_ID:
  			ofile.write(i+'\n')
  		ofile.close()
  	elif len(single_ID) == 0:
  		ofile1 = open('annotation/list', 'a')
  		ofile2 = open('annotation/list_n_c_f_m', 'a')
  		for k in ped_ID:
  			a = ' '.join(ped_ID[k])
  			if k+u'父' in a and k+u'母' in a:
  				ofile1.write(d1.loc[k]['sample']+'\n')
  			else:
  				ofile2.write(d1.loc[k]['sample']+'\n')
  		ofile1.close()
  		ofile2.close()
  	else:
  		ofile = open('annotation/list', 'a')
  		for i in single_ID:
  			ofile.write(i+'\n')
  		ofile.close()
  	if os.path.getsize('annotation/list') == 0:
  		os.system('rm annotation/list')
  	try:
  		if os.path.getsize('annotation/list_n_c_f_m') == 0:
  			os.system('rm annotation/list_n_c_f_m')
  	except:
  		print 'no not_child_father_mother'
  generate_annotation()

  ofile = open('step3_annotation.sh', 'a')
  if len(single_ID) > 0 and len(trio_ID) > 0:
  	for i in single_ID:
  		ofile.write('mv '+i+'/2_mapping/'+i+'.depth.sample_gene_summary annotation\n')
  		ofile.write('mv '+i+'/3_variants/'+i+'.{ann*,CADD,link,maf} annotation\n')
  	for i in trio_ID:
  		ofile.write('mv '+i+'/2_mapping/'+i+'.depth.sample_gene_summary annotation\n')
  		ofile.write('mv trio/'+i+'/segtrio/'+i+'.{ann*,CADD,link,maf} annotation\n')
  elif len(trio_ID) == 0:
  	for i in single_ID:
  		ofile.write('mv '+i+'/2_mapping/'+i+'.depth.sample_gene_summary annotation\n')
  		ofile.write('mv '+i+'/3_variants/'+i+'.{ann*,CADD,link,maf} annotation\n')
  else:
  	for i in trio_ID:
  		ofile.write('mv '+i+'/2_mapping/'+i+'.depth.sample_gene_summary annotation\n')
  		ofile.write('mv trio/'+i+'/segtrio/'+i+'.{ann*,CADD,link,maf} annotation\n')
  #ofile.write('cp -r /DATA/sslyu/trio_BWA-GATK_2.7/annotation/db/ annotation\n')
  ofile.write('cd annotation\n')
  if os.path.exists('annotation/list_n_c_f_m'):
  	ofile.write('cat list* > list_all\n')
  ofile.write('python /DATA/sslyu/trio_BWA-GATK_2.7/annotation/score_re.py\n')
  if os.path.exists('annotation/list'):
  	ofile.write('nohup python /DATA/sslyu/trio_BWA-GATK_2.7/annotation/annotation_filt_ver2.8.py -l list -thisminedge '+sys.argv[1]+'&\n')
  if os.path.exists('annotation/list_n_c_f_m'):
  	ofile.write('nohup python /DATA/sslyu/trio_BWA-GATK_2.7/annotation/annotation_filt_ver2.8.py -l list_n_c_f_m -c_f_m no -thisminedge '+sys.argv[1]+'-c_f_m no &\n')
  ofile.write('cd ..\n')
  ofile.close()

with open('step3_annotation.sh') as f:
	for l in f:
		print l.strip('\n')
truth = raw_input('Is the shell script correct? If not, please enter no: ')
if truth == 'no':
	sys.exit()
else:
	os.system('sh step3_annotation.sh')
