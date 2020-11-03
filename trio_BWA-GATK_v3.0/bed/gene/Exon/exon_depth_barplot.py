#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
version:3.0_20190124
by:lss
update:20190124
"""
ThisVer='v3.0_20190124'
Thisreleasetime='2019/1/24'

import os
import sys
import pandas as pd
import numpy as np

m_path = os.path.split(os.path.abspath(sys.argv[0]))[0]

if len(sys.argv) == 1:
    print """
    Note: there should be two dirs depth and tsv in the current working directory, the depth file of the to be processing samples should be under depth with the format names of sample.depth
    Usage: python %s/exon_depth_barplot.py -gene genesymbol -sn sample
           python %s/exon_depth_barplot.py -gene genesymbol -l list
    """ % (m_path, m_path)
    sys.exit()

##read parameters of cmd line
kwargs_raw = sys.argv[1:]
kwargs={'-gene':'','-l':'list','-sn':''}
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
    list=pd.read_table(kwargs.get('-l'),header=None)[0].tolist()
    print "\nNote:list samples mode"

gene = kwargs.get('-gene')

def exon_plot(gene, sn):
    d = pd.DataFrame()
    av_depth = []
    exon = []
    transcript = []
    chromosome = []
    exonStart = []
    exonEnd = []
    ofile = open(sn+'_'+gene+'_exon.depth.csv', 'w')
    ofile.write('sample,gene,exon ID,hgmd,chromosome,exonStart,exonEnd,depth\n')
    with open('/DATA/sslyu/trio_BWA-GATK_v3.0/bed/gene/Exon/'+gene+'.bed') as f:
        n = 0
        for l in f:
            l = l.strip('\n').split('\t')
            hgmd = l[-1]
            chrom = l[0]
            start = l[1]
            end = l[2]
            cdslen = int(l[2])-int(l[1])
            count = 0
            n = n+1
            with open(sn+'_'+gene+'.bed.depth') as f1:
                for l1 in f1:
                    l1 = l1.strip('\n').split('\t')
                    if l1[0] == chrom:
                        if l1[1] >= start and l1[1] <= end:
                            count = count+int(l1[2])
            depth = float(count)/float(cdslen)
            av_depth.append('{:.2f}'.format(depth))
            #exon.append(sn+'_'+gene+'_exon'+str(n))
            exon.append(l[4])
            transcript.append(l[3])
            chromosome.append(chrom)
            exonStart.append(start)
            exonEnd.append(end)

    for i in range(len(exon)):
        ofile.write(sn+','+transcript[i]+','+exon[i]+','+hgmd+','+chromosome[i]+','+exonStart[i]+','+exonEnd[i]+','+av_depth[i]+'\n')
    ofile.close()

for sn in list:
    sn = str(sn)
    exon_plot(gene, sn)

