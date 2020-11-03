#!/usr/bin/env python
# -*- coding: utf-8 -*-
#########################################################################
# File Name: extract_exon_region.py
# Author: sslyu
# mail: 1573077420@qq.com
# Created Time:Fri 10 May 2019 02:48:54 PM CST
#########################################################################

import numpy as np
import pandas as pd
import subprocess
import sys

def extract_exon_region(gene):
    #p1 = subprocess.Popen(['grep', transcript_ID, 'ensemble_GRCh37'],stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    #info = subprocess.Popen(['sort'], stdin = p1.stdout, stdout=subprocess.PIPE).communicate()[0].split('\n')[:-1]
    info = pd.read_csv('/DATA/sslyu/trio_BWA-GATK_v3.0/bed/gene/Exon/ensemble_GRCh37', sep="\t")
    info = info[info['Gene name'] == gene]
    if info.shape[0] > 0:
        max_transcript = max(info['Transcript length (including UTRs and CDS)'])
        info = info[info['Transcript length (including UTRs and CDS)'] == max_transcript]
    
        #gene = info[0].split('\t')[0]
        #print gene
        #hgmd = subprocess.Popen(['grep', gene, '/DATA/ana005/annovar/HGMD/HGMD_2018.spring/HGMD_Advanced_Micro_Lesions.csv'],stdout=subprocess.PIPE, stderr=subprocess.STDOUT).communicate()[0].split('\n')[0].split(',')[4].split(':')[0]
        #if not hgmd.startswith('NM'):
        #    hgmd = subprocess.Popen(['grep', gene, '/DATA/ana005/annovar/HGMD/HGMD_2018.spring/HGMD_Advanced_Micro_Lesions.csv'],stdout=subprocess.PIPE, stderr=subprocess.STDOUT).communicate()[0].split('\n')[1].split(',')[4].split(':')[0]
        hgmd = pd.read_csv('/DATA/sslyu/trio_BWA-GATK_v3.0/bed/gene/gene_hgvc.csv', header=None)
        hgmd = hgmd[hgmd[0] == gene]
        hgvs = ";".join(hgmd[1].unique().tolist())
        info['hgvs'] = hgvs
        info['Chromosome/scaffold name'] = info['Chromosome/scaffold name'].apply(lambda x:"chr"+str(x))
        chroms=[]
        for i in range(22):
            chroms.append("chr"+str(i+1))
        
        chroms.append('chrX')
        chroms.append('chrY')
        for i in info.index:
            if info.loc[i,'Chromosome/scaffold name'] not in chroms:
                info = info.drop(i)
        info = info.sort_values(by=['Exon region start (bp)'])
        output_header = ['Chromosome/scaffold name', 'Exon region start (bp)', 'Exon region end (bp)', 'Gene name', 'Exon stable ID', 'hgvs']
        info[output_header].to_csv('/DATA/sslyu/trio_BWA-GATK_v3.0/bed/gene/Exon/'+gene+'.bed', sep="\t", header=None, index=None)
        print gene
        print info
        print
        #ofile = open(gene+'.bed','w')
        #for i in info:
        #    i = i.split('\t')
        #    chrom = i[3]
        #    exon = i[5]
        #    start = i[1]
        #    end = i[2]
        #    ofile.write("chr"+chrom+"\t"+start+"\t"+end+"\t"+gene+"\t"+exon+"\t"+",".join(hgvs)+"\n")
        #ofile.close()
    else:
        print "%s not in /DATA/sslyu/trio_BWA-GATK_v3.0/bed/gene/Exon/ensemble_GRCh37" % gene

gene = sys.argv[1]
if gene.endswith('.bed'):
    gene = gene[:-4]
extract_exon_region(gene)
