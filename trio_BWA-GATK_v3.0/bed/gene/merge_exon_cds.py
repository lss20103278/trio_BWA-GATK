#!/usr/bin/env python
# -*- coding: utf-8 -*-
#########################################################################
# File Name: merge_exon_cds.py
# Author: sslyu
# mail: 1573077420@qq.com
# Created Time:Mon 13 May 2019 02:09:24 PM CST
#########################################################################

import numpy as np
import pandas as pd
import sys
import os

def merge_exon_cds(gene):
    def cds_range():
        CDS_start_min = 0
        CDS_end_max = 0
        CDS = pd.read_csv('/DATA/sslyu/trio_BWA-GATK_v3.0/bed/gene/CDS/'+gene+'.bed', sep="\t", header=None)
        CDS_start_min = min(CDS[1])
        CDS_end_max = max(CDS[2])
        chrom = str(CDS[0][0])
        return [chrom,CDS_start_min,CDS_end_max]
    def exon_range():
        exon_start_min = 0
        exon_end_max = 0
        exon = pd.read_csv('/DATA/sslyu/trio_BWA-GATK_v3.0/bed/gene/Exon/'+gene+'.bed', sep="\t", header=None)
        exon_start_min = min(exon[1])
        exon_end_max = max(exon[2])
        chrom = str(exon[0][0])
        return [chrom,exon_start_min,exon_end_max]
    def generate_bed(gene,chrom,start,end):
        ofile = open('/DATA/sslyu/trio_BWA-GATK_v3.0/bed/gene/'+gene+'.bed', 'w')
        ofile.write(chrom+"\t"+start+"\t"+end+"\n")
        ofile.close()
    if os.path.exists('/DATA/sslyu/trio_BWA-GATK_v3.0/bed/gene/CDS/'+gene+'.bed') and os.path.exists('/DATA/sslyu/trio_BWA-GATK_v3.0/bed/gene/Exon/'+gene+'.bed'):
        cds = cds_range()
        exon = exon_range()
        chrom = cds[0]
        start = str(min([cds[1],exon[1]]))
        end = str(max([cds[2],exon[2]]))
        generate_bed(gene,chrom,start,end)
    elif os.path.exists('/DATA/sslyu/trio_BWA-GATK_v3.0/bed/gene/CDS/'+gene+'.bed'):
        cds = cds_range()
        chrom = cds[0]
        start = str(cds[1])
        end = str(cds[2])
        generate_bed(gene,chrom,start,end)
    elif os.path.exists('/DATA/sslyu/trio_BWA-GATK_v3.0/bed/gene/Exon/'+gene+'.bed'):
        exon = exon_range()
        chrom = exon[0]
        start = str(exon[1])
        end = str(exon[2])
        generate_bed(gene,chrom,start,end)
    else:
        print "%s both CDS and exon are lost" % gene

gene = sys.argv[1]
if gene.endswith('.bed'):
    gene = gene[:-4]
merge_exon_cds(gene)
