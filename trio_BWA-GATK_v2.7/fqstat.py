#!/usr/bin/env python

import sys
import os
import pandas as pd

ReadNum = []
BaseNum = []
Q30 = []
sn = pd.read_table('list', header=None)
ofile = open('fqstat.sh', 'w')
for i in sn[0]:
	ofile.write('/home/ana005/anaconda2/bin/iTools Fqtools stat -InFq raw/'+i+'*R1* -InFq raw/'+i+'*R2* -OutStat raw/'+i+'.info\n')
ofile.close()
os.system('bash fqstat.sh')
os.system('rm fqstat.sh')
for i in sn[0]:
	with open('raw/'+i+'.info') as f:
		Q30_base = 0.0
		for l in f:
			if 'ReadNum' in l:   
				readnum = l.split(': ')[1].split('\t')[0]
				if 'BaseNum' in l:
					basenum = int(l.split(': ')[2].split('\t')[0])*2/1000000
			if 'BaseQ' in l and 'Q30' in l:
				Q30_base = Q30_base + float(l.strip('\n').split('Q30: ')[1][:-1])
		ReadNum.append(readnum)
		BaseNum.append(basenum)
		Q30.append('{:.2f}'.format(Q30_base/2)+'%')
df = pd.DataFrame(columns=['reads', 'bases(Mb)', 'Q30'])
df['reads'] = ReadNum
df['bases(Mb)'] = BaseNum
df['Q30'] = Q30
df.index = sn[0]
df.to_csv('all_fqstat', sep="\t")
