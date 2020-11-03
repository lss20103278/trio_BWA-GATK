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
import re

father_ID = {}
mother_ID = {}
phenotype2 = {}

def append_sample_gender_excel(infile):
    dataframe = pd.read_excel(infile, dtype=str)
    ID = []
    for i in range(len(dataframe[u'原始样本ID'].isnull())):
        if d1[u'原始样本ID'].isnull()[i]:
            ID.append(dataframe.iloc[i][u'样本编号'])
        else:
            ID.append(dataframe.iloc[i][u'原始样本ID'])
    dataframe['sample'] = ID
    dataframe.index = dataframe['sample']
    dataframe['gender'] = d1[u'性别'].apply(lambda x:'1' if u'男' in x else '2')
    return d1
def append_trio_dataframe(dataframe):
    dataframe.index = dataframe[u'姓名']
    if 'pedigree' not in dataframe.columns:
        pedigree = {}
        for i in dataframe[u'姓名']:
            pedigree[i] = []
        for k in pedigree:
            for i in dataframe[u'姓名']:
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
        pedigree = []
        phenotype = []
        for i in dataframe.index:
            for j in trio:
                if j in i:
                    pedigree.append(j)
                    if j == i:
                        phenotype.append('2')
                    else:
                        phenotype.append('1')
        dataframe['pedigree'] = pedigree
        dataframe['phenotype'] = phenotype

def append_father_mother(dataframe):


def generate_ped_ID(d1):
    ped_ID = {}
    if 'relationship' in d1.columns:
        trio = []
        for i in d1['pedigree']:
            if i not in trio:
                trio.append(i)
        for i in trio:
            ped_ID[i] = {}
        for k in ped_ID:
            for i in d1.index:
                if d1.loc[i]['pedigree'] == k:
                    if d1.loc[i]['relation'] == u'子':
                        ped_ID[k].update({u'子':i})
                    elif d1.loc[i]['relation'] == u'父':
                        ped_ID[i].update({u'父':i})
                    elif d1.loc[i]['relation'] == u'母':
                        ped_ID[i].update({u'母':i})
                    else:
                        ped_ID[i].update({d1.loc[i][u'姓名']:i})
    else:
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
        for i in trio:
            ped_ID[i] = {}
        for k in ped_ID:
            ped_ID[k] = find_father_mother(k)
    return ped_ID

def find_father_mother(d1,k):
    family = {}
    d1.index = d1[u'姓名']
    ped_ID_k = []
    for i in d1.index:
        a = ''
        if k in i:
            a = a+i
            ped_ID_k.append(i)
    if k+u'父' in a and k+u'母' in a:
        key = d1.loc[k]['sample']
        for i in ped_ID_k:
            if re.match(k+u'父',i):
    	        father_name = i
            if re.match(k+u'母',i):
    	        mother_name = i
        family[u'父'] = d1.loc[father_name]['sample']
        family[u'母'] = d1.loc[mother_name]['sample']
    else:
        try:
            father = raw_input('father name of '+i+'(if not exist, please enter 0):')
        father_ID[key] = d1.loc[father_name]['sample']
        mother_ID[key] = d1.loc[father_name]['sample']
        phenotype2[key] = u'2'
        if len(ped_ID_k) == 3:
            for i in [u'父', u'母']:
                for j in ped_ID_k:
                    if re.match(k+i, j):
                        name = j
                key = d1.loc[k]['sample']
                print key
                father_ID[key] = u'0'
                mother_ID[key] = u'0'
                phenotype2[key] = u'1'
        else:
            n_c = [i for i in ped_ID_k if i != k]
            for i in n_c:
                key = d1.loc[i]['sample']
                father = raw_input('father name of '+i+'(if not exist, please enter 0):')
                try:
                    father_ID[key] = d1.loc[father]['sample']
                except:
                    father_ID[key] = '0'
                mother = raw_input('mother name of '+i+'(if not exist, please enter 0):')
                try: 
                    mother_ID[key] = d1.loc[mother]['sample']
                except:
                    mother_ID[key] = '0'
                phenotype2[key] = '1'
    else:
        for i in ped_ID_k:
            key = d1.loc[i]['sample']
            father = raw_input('father name of '+i+'(if not exist, please enter 0):')
            try:
                father_ID[key] = d1.loc[father]['sample']
            except:
                father_ID[key] = u'0'
            mother = raw_input('mother name of '+i+'(if not exist, please enter 0):')
            try:
                mother_ID[key] = d1.loc[mother]['sample']
            except:
                mother_ID[key] = u'0'
            if i == k:
                phenotype2[key] = u'2'
            else:
                phenotype2[key] = u'1'

        if k in i:
            if k == i:
                family[u'子'] = d1.loc[i]['sample']
            elif k in i and u'父' in i:
                family[u'父'] = d1.loc[i]['sample']
            elif k in i and u'母' in i:
                family[u'母'] = d1.loc[i]['sample']
            else:
                family[i] = d1.loc[i]['sample']
    return family

def generate_ped(k1,k2):
    if 'relation' in d1.columns:
        father_ID = {}
        mother_ID = {}
        phenotype2 = {}
        trio = []
        for i in d1['pedigree']:
            if i not in trio:
                trio.append(i)
        for i in trio:
            
    else:
        a = ' '.join(ped_ID[k])
        print a
    d1.index = d1['sample']
    d1[u'father'] = pd.Series(father_ID)
    d1[u'mother'] = pd.Series(mother_ID)
    d1[u'phenotype1'] = u'0'
    d1[u'phenotype2'] = pd.Series(phenotype2)
    peddy = open('trio/ped/peddy.ped', 'a')
    peddy_ID = '0'
    
    ped_f = d1[['sample', u'father', u'mother', u'gender', u'phenotype1']]
    ped_mendel = d1[['sample', u'father', u'mother', u'gender', u'phenotype2']]
    for k in ped_ID:
        f1name = 'trio/ped/'+child[k]+'.ped'
        f2name = 'trio/ped/'+child[k]+'.mendel.ped'
        ped1 = open(f1name, 'w')
        ped2 = open(f2name, 'w')
        a = ' '.join(ped_ID[k])
        key1 = d1['sample'][list(d1[u'姓名']).index(k)]
        ped1.write(key1+' '+' '.join(ped_f.loc[key1])+'\n')
        ped2.write(key1+' '+' '.join(ped_mendel.loc[key1])+'\n')
        peddy.write(peddy_ID+' '+' '.join(ped_mendel.loc[key1])+'\n')
        if k+u'父' and k+u'母' in a:
        	c_f_m = [k+u'父', k+u'母']
        	for i in c_f_m:
        	    for j in ped_ID[k]:
        	        if re.match(i, j):
        		        key2 = d1['sample'][list(d1[u'姓名']).index(j)]
        		        ped1.write(key1+' '+' '.join(ped_f.loc[key2])+'\n')
        		        ped2.write(key1+' '+' '.join(ped_mendel.loc[key2])+'\n')
        		        peddy.write(peddy_ID+' '+' '.join(ped_mendel.loc[key2])+'\n')
        	n_c_f_m = [i for i in ped_ID[k] if i not in c_f_m and i != k]
        	for i in n_c_f_m:
        	    key2 = d1['sample'][list(d1[u'姓名']).index(i)]
        	    ped1.write(key1+' '+' '.join(ped_f.loc[key2])+'\n')
        	    ped2.write(key1+' '+' '.join(ped_mendel.loc[key2])+'\n')
        	    peddy.write(peddy_ID+' '+' '.join(ped_mendel.loc[key2])+'\n')
        else:
        	n_c = [i for i in ped_ID[k] if i != k]
        	for i in n_c:
        	    key2 = d1['sample'][list(d1[u'姓名']).index(i)]
        	    ped1.write(key1+' '+' '.join(ped_f.loc[key2])+'\n')
        	    ped2.write(key1+' '+' '.join(ped_mendel.loc[key2])+'\n')
        	    peddy.write(peddy_ID+' '+' '.join(ped_mendel.loc[key2])+'\n')
        ped1.close()
        ped2.close()
    peddy.close()
if len(sys.argv) < 2:
    print 'missing the value of filtmode...'
    sys.exit()

filelist = os.listdir('.')
infile = []
for i in filelist:
    if '分析任务单' in i:
        infile.append(i)
    if '产品信息表' in i:
        infile.append(i)

if len(infile) != 0:
    for sn in infile:
    	print sn
    	d1 = pd.read_excel(sn)
        ID = []
        for i in range(len(d1[u'原始样本ID'].isnull())):
       	    if d1[u'原始样本ID'].isnull()[i]:
                ID.append(d1.iloc[i][u'样本编号'])
            else:
                ID.append(d1.iloc[i][u'原始样本ID'])
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
        print trio_m_ID
        judgement = raw_input('Do you need to copy gvcf? If not, please enter no:')
        if judgement != 'no':
            ofile = open('cp.sh', 'w')
            for i in trio_m_ID:
             	ofile.write('cp '+i+'/3_variants/gvcf/* trio/gvcf\n')
            ofile.close()
            with open('cp.sh') as f:
             	for l in f:
                    print l.strip('\n')
            judgement1 = raw_input('Is the shell script correct? If not, please enter no:')
            if judgement1 == 'no':
             	cmd = raw_input('Please enter the copy command to copy gvcf:')
             	os.system(cmd)
            else:
                os.system('bash cp.sh')
                os.system('rm cp.sh')
        
        child = {}
        ped_ID = {}
        for i in trio:
            child[i] = d1.loc[i]['sample']
            for j in d1[u'姓名']:
                if i in j:
            	    ped_ID.setdefault(i,[]).append(j)
        
        generate_ped()
        ofile = open('trio/list', 'a')
        for i in child:
        	ofile.write(child[i]+'\n')
        ofile.write('peddy\n')
        ofile.close()
else:
    print 'please offer an information excel, the name of the excel should include \"分析任务单\", the columns should be 样本编号 原始样本ID 姓名 性别'
    sys.exit()


os.chdir('trio')
os.system('sort ped/peddy.ped > tmp |mv tmp ped/peddy.ped')
def generate_cmd():
	ofile = open('cmd.sh', 'w')
	n = 0
	with open('/DATA/sslyu/trio_BWA-GATK_2.7/trio/cmd.sh') as f:
		for l in f:
			n = n+1
			if n == 12:
				ofile.write('filtmode=\''+sys.argv[1]+'\'\n')
			else:
				ofile.write(l)
	ofile.close()
generate_cmd()
with open('ped/peddy.ped') as f:
	for l in f:
		print l.strip('\n')
with open('cmd.sh') as f:
	for l in f:
		print l.strip('\n')
truth = raw_input('Is the peddy.ped correct? If not, please enter no: ')
if truth != '':
	os.system('rm ped/*')
	os.chdir('..')
	sys.exit()
else:
	os.system('sh sub.sh')
	os.chdir('..')
