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

filelist = os.listdir('.')
infile = ''
for i in filelist:
	if '分析任务单' in i:
		seperator = '分析任务单'
		break
	if '产品信息表' in i:
		seperator = '产品信息表'
		break
infile = i
d1 = pd.read_excel(infile)
ID = []
# this block has a problem, wait to be fixed
#for i in range(d1.shape[0]):
#	if np.isnan(d1.iloc[i][u'原始样本ID']):
#		ID.append(d1.iloc[i][u'样本编号'])
#	else:
#		ID.append(d1.iloc[i][u'原始样本ID'])
for i in range(len(d1[u'原始样本ID'].isnull())):
    if d1[u'原始样本ID'].isnull()[i]:
        ID.append(d1.iloc[i][u'样本编号'])
    else:
        ID.append(d1.iloc[i][u'原始样本ID'])
d1['sample'] = ID
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
ped_v = []
for i in d1[u'姓名']:
	if i in single:
		ped_v.append(0)
	else:
		ped_v.append(1)
d1['pedigree'] = ped_v
d1[u'分析流程'] = d1['pedigree'].map(lambda x:'single' if x==0 else 'trio')

relationship = []
for i in d1[u'姓名'].tolist():
	if i in single:
		relationship.append(np.nan)
	else:
		for k in trio:
			if k in i:
				relationship.append('fam'+d1.loc[i]['sample'])

ped_relation = []
if len(trio) > 0 and len(single) > 0: 
	for i in d1[u'姓名'].tolist():
		if i in single:
			ped_relation.append(np.nan)
		else:
			for k in trio:
				if k == i:
					ped_relation.append(u'子')
				elif k != i and k in i:
					v = i[len(k):]
					ped_relation.append(v)
elif len(single) == 0:
	for i in d1[u'姓名'].tolist():
		for k in trio:
			if k == i:
				ped_relation.append(u'子')
			elif k != i and k in i:
				v = i[len(k):]
				ped_relation.append(v)
else:
	for i in d1[u'姓名']:
		ped_relation.append(np.nan)

d2 = pd.DataFrame(columns=['sample', '1Xcov', '20Xcov', 'AvgDepth', 'duplication'])
for i in d1['sample']:
	fn = i+'/2_mapping/'+i+'.qcsum'
	qcsum = pd.read_table(fn)
	d2 = pd.concat([d2,qcsum], ignore_index=True)
d2.index = d2['sample']
d2.to_csv('all.qcsum', sep="\t", index=None)

d3 = pd.DataFrame(columns=['sample', 'indel_genomic', u'snp_genomic', u'ts/tv_genomic',u'snp_exonic', u'ts/tv_exonic'])
qcvcf_list = []
for i in d1[u'姓名'].tolist():
	if i in single:
		qcvcf_list.append(d1.loc[i]['sample'])
	elif i in trio:
		qcvcf_list.append(d1.loc[i]['sample'])
	else:
		continue
d1.index = d1['sample']
for i in qcvcf_list:
	if d1.loc[i][u'姓名'] in single:
		fn = i+'/3_variants/'+i+'.qcvcf'
	else:
		fn = 'trio/'+i+'/triotmp/'+i+'.qcvcf'
	qcvcf = pd.read_table(fn)
	d3 = pd.concat([d3,qcvcf], ignore_index=True)
d3.index = d3['sample']

for i in [u'reads', u'bases(Mb)', u'Q30']:
	if d1[i].isnull().any():
		print i+' is missing, fqstat...'
		os.system('python fqstat.py')
		break
if os.path.exists('all_fqstat'):
	d4 = pd.read_table('all_fqstat', index_col=0)
	d1[[u'reads', u'bases(Mb)', u'Q30']] = d4
	
d5 = pd.merge(d1[[u'批次', u'样本编号', u'sample', u'姓名', u'性别', u'年龄', u'科室', u'分析流程', u'reads', u'bases(Mb)', u'Q30', u'浓度(ng/ul)', u'体积(ul)', u'质量(ug)', u'OD260/OD280', u'OD260/OD230']], d2.drop(u'duplication', axis=1),on="sample")
strategy = os.getcwd().split('/')[-1].split('_')[1]
print 'strategy: '+strategy
truth = raw_input('Is the strategy right? If not, please enter it (WES PANEL Agilent_wes ... ): ')
if truth != '':
	strategy = truth
l1 = ['sample', u'年龄', u'分析流程', 'vcf']
l2 = [u'策略', u'样本间关系', 'vcf', 'QC']
l3 = [strategy, relationship, 'available', u'合格']
for i in range(len(l1)):
	d5.insert(d5.columns.tolist().index(l1[i])+1,l2[i],l3[i])
col1 = d5.columns.tolist()[:12]+d5.columns.tolist()[20:]+d5.columns.tolist()[12:20]
d6 = d5[col1]
colname2 = d6.columns.tolist()
colname2[2] = u'样本ID'
colname2[8] = u'申请科室'
colname2[12:15] = [u'1X_coverage(%)', u'20X_coverage(%)', u'ana_avg_depth(X)']
d6.columns = colname2
ofile = infile.split(seperator)[0]+'交付信息表'+infile.split(seperator)[1]
wb=pd.ExcelWriter(ofile,engine='openpyxl')
d6.to_excel(wb,index=False)
wb.save()

d6 = pd.read_excel('/DATA/sslyu/trio_BWA-GATK_2.7/navicat.xlsx')
infile2 = infile.split(seperator)[0]+'交付信息表'+infile.split(seperator)[1]
d5 = pd.read_excel(infile2)
colname3 = d6.columns
sample_n = d1.shape[0]
start = raw_input('please enter the start number of this batch of samples in the database:')
d6[u'order'] = range(int(start), int(start)+sample_n)
d6[[u'批次', u'样本编号', u'样本ID', u'策略', u'姓名', u'性别', u'年龄',u'样本间关系', u'申请科室',u'exp_reads',u'exp_bases(Mb)', u'exp_Q30', u'exp_浓度(ng/ul)', u'exp_体积(ul)',u'exp_质量(ug)', u'exp_OD260/OD280', u'exp_OD260/OD230', u'exp_QC结论',u'ana_qc_comment',u'ana_pepline', u'ana_vcf']] = d5[[u'批次', u'样本编号', u'样本ID', u'策略', u'姓名', u'性别', u'年龄',u'样本间关系', u'申请科室',u'reads',u'bases(Mb)',u'Q30',u'浓度(ng/ul)',u'体积(ul)',u'质量(ug)',u'OD260/OD280',u'OD260/OD230',u'QC',u'QC',u'分析流程',u'vcf']]
d6.index = d1.index
d6[[u'process_送样日期', u'process_上机日期', u'process_下机日期']] = d1[[u'收样日期', u'上机日期',u'下机日期']]
d6[[u'ana_coverage1X(%)',u'ana_coverage20X(%)', u'ana_avg_depth(X)', u'ana_duplication(%)']] = d2[[u'1Xcov', u'20Xcov', u'AvgDepth', u'duplication']]
d6[[u'ana_indel(genomic)', u'ana_snp(genomic)', u'ana_ts/tv(genomic)',u'ana_snp(exonic)', u'ana_ts/tv(exonic)']] = d3[[u'indel_genomic', u'snp_genomic', u'ts/tv_genomic',u'snp_exonic', u'ts/tv_exonic']]
d6[u'ana_annotation'] = u'v2.8'
d6[u'ana_reporter'] = u'吕珊珊'
d6[u'bk_rawdata'] = u'/anaData/anaData004/children_hos_genetic/rawdata/2018/WES/'
dirname = '/'.join(os.getcwd().split('/')[-2:])
d6[u'bk_reports'] = u'/DATA/BPshare/opm/遗传病-交大附属儿医/2018/'+dirname

d6[u'家系关系'] = ped_relation

ofile = infile.split(seperator)[0]+'navicat_'+infile.split(seperator)[1]
wb=pd.ExcelWriter(ofile, engine='openpyxl')
d6.to_excel(wb,index=False)
wb.save()
