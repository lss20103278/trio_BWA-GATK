#!/usr/bin/env python
# -*- coding: utf-8 -*-
#########################################################################
# File Name: merge_interval.py
# Author: sslyu
# mail: 1573077420@qq.com
# Created Time:Thu 11 Apr 2019 03:50:55 PM CST
#########################################################################

import numpy as np
import pandas as pd
from collections import *
#this script is to merge different intervals of the same gene into a comprehensive interval

dct = {}
with open('ucscgene_genesymbol.bed') as f:
    for l in f:
        l = l.strip('\n').split('\t')
        dct[l[3]] = {'chr':[],'start':[],'end':[]}

with open('ucscgene_genesymbol.bed') as f:
    for l in f:
        l = l.strip('\n').split('\t')
        dct[l[3]]['chr'].append(l[0])
        dct[l[3]]['start'].append(int(l[1]))
        dct[l[3]]['end'].append(int(l[2]))

ofile = open('ucscgene_genesymbol.uniq.bed', 'w')
for k in dct:
    chr = Counter(dct[k]['chr']) # Counter({'chr1': 3})
    for k_chr in chr:
        ix = []
        for i in range(len(dct[k]['chr'])):
            if dct[k]['chr'][i] == k_chr:
                ix.append(i)
        start_chr = []
        end_chr = []
        for i in ix:
            start_chr.append(dct[k]['start'][i])
            end_chr.append(dct[k]['end'][i])
        start = min(start_chr)
        end = max(end_chr)
        chrom = k_chr
        start = str(start)
        end = str(end)
        ofile.write(chrom+'\t'+start+'\t'+end+'\t'+k+'\n')
ofile.close()
