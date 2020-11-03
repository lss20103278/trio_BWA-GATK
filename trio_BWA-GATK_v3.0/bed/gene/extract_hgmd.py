#!/usr/bin/env python
# -*- coding: utf-8 -*-
#########################################################################
# File Name: extract_hgmd.py
# Author: sslyu
# mail: 1573077420@qq.com
# Created Time:Tue 14 May 2019 04:17:59 PM CST
#########################################################################

import numpy as np
import pandas as pd

hgmd = pd.read_table('/DATA/sslyu/soft/annovar/humandb/hg19_hgmd_2018_spring.txt')
for i in hgmd.index:
    if 'DNA=' not in hgmd.loc[i,'hgmd']:
        hgmd = hgmd.drop(i)
hgmd = hgmd['hgmd'].str.split(';', expand=True)
hgmd[3] = hgmd[3].str.split('=', expand=True)[1]
hgmd[5] = hgmd[5].str.split(':', expand=True)[0]
hgmd[5] = hgmd[5].str.split('=', expand=True)[1]
hgmd = hgmd[[3,5]]
hgmd.drop_duplicates(inplace = True)
hgmd.to_csv('gene_hgvc.csv', index=None, header=None)
