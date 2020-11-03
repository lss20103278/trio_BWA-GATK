#!/usr/bin/env python
# -*- coding: utf-8 -*-
#########################################################################
# File Name: exon_depth_barplot_v5.py
# Author: sslyu
# mail: 1573077420@qq.com
# Created Time:Fri 15 Feb 2019 11:03:28 AM CST
#########################################################################

import numpy as np
import pandas as pd
import sys
import os

ThisVer='v5'
Thisreleasetime='2019/2/15'

m_path = os.path.split(os.path.abspath(sys.argv[0]))[0]

if len(sys.argv) == 1:
    print """
    Note: prepare the depth file of samples
    Usage: python %s/exon_depth_barplot.py -gene genename -sn sampleID
           python %s/exon_depth_barplot.py -gene genename -l list
    """ %(m_path, m_path)
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

if not os.path.exists(gene+'.loci'):
    print """
    grep -w %s /DATA/sslyu/soft/annovar/humandb/refseq.sorted.txt |cut -f 3,10-11 |awk '{l=split($2,a,","); split($3,b,","); for (i=1;i<l;i++){OFS="\\t"; print $1,a[i],b[i]}}' |sort |uniq > %s.bed
    awk '{for (i=$2;i<$3;i++){OFS="\\t"; print $1,i,"exon"NR}}' %s.bed > %s.loci
    """ %(gene, gene, gene, gene)
    sys.exit()

def exon_plot(sample, gene):
    loci = pd.read_table(gene+'.loci', header=None)
    loci[1] = loci[1].apply(str)
    depth = pd.read_table(sample+'.depth')
    loci['Locus'] = loci[0]+':'+loci[1]
    d = pd.merge(loci, depth, how='left', on='Locus')
    count = d.groupby(2)['sample_Average_Depth'].sum()
    count = count.fillna(0)
    cdslen = loci.groupby(2)['Locus'].count()
    mean_depth = count.divide(cdslen)
    mean_depth = mean_depth.apply(lambda x:round(x,2))
    mean_depth.to_csv(sample+'_'+gene+'_exon.depth.tsv', header=None, sep="\t")

for sn in list:
    exon_plot(sn, gene)
    #os.system('Rscript /DATA/sslyu/trio_BWA-GATK_3.0/src/exon_depth_barplot.R '+sn+'_'+gene+'_exon.depth.tsv '+gene+' '+sn)
