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

def extract_cds_region(cds_ID):
    p1 = subprocess.Popen(['grep', "-w", cds_ID, 'CCDS.20110907.hg19.txt'],stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    info = subprocess.Popen(['sort'], stdin=p1.stdout, stdout=subprocess.PIPE).communicate()[0].strip('\n').split('\t')
    chrom = info[0]
    gene = info[2]
    cds = info[4]
    start = []
    end = []
    for i in info[9].split(', '):
        if i.startswith('['):
            i = i[1:]
        elif i.endswith(']'):
            i = i[:-1]
        else:
            i = i
        start.append(i.split('-')[0])
        end.append(i.split('-')[1])
    hgmd = subprocess.Popen(['grep', gene, '/DATA/ana005/annovar/HGMD/HGMD_2018.spring/HGMD_Advanced_Micro_Lesions.csv'],stdout=subprocess.PIPE, stderr=subprocess.STDOUT).communicate()[0].split('\n')[0].split(',')[4].split(':')[0]
    ofile = open(gene+'.bed','w')
    for i in range(len(start)):
        ofile.write("chr"+chrom+"\t"+start[i]+"\t"+end[i]+"\t"+gene+"\t"+cds+"\t"+hgmd+"\n")
    ofile.close()


def extract_start_end(j):
    if j.startswith('[') and j.endswith(']'):
        j = j[1:-1]
    elif j.startswith('['):
        j = j[1:]
    elif j.endswith(']'):
        j = j[:-1]
    else:
        j = j
    return j

def extract_cds_region(gene):
    p1 = subprocess.Popen(['grep', "-w", gene, '/DATA/sslyu/trio_BWA-GATK_v3.0/bed/gene/CDS/CCDS.20110907.hg19.txt'],stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    info = subprocess.Popen(['sort'], stdin=p1.stdout, stdout=subprocess.PIPE).communicate()[0].split('\n')
    if info[0] != "":
        info = info[:-1] # the last element is ""
        start = []
        end = []
        cds = []
        if len(info) == 1:
            i = info[0]
            i = i.split('\t')
            chrom = i[0]
            for j in i[9].split(', '):
                cds.append(i[4])
                j = extract_start_end(j)
                start.append(j.split('-')[0])
                end.append(j.split('-')[1])
        else:
            cds_region = {}
            for i in info:
                i = i.split('\t')
                if len(i) > 1:
                    for j in i[9].split(', '):
                        j = extract_start_end(j)
                        chrom = i[0]
                        #gene = i[2]
                        cds = i[4]
                        if cds_region.has_key(j):
                            cds_region[j].append(cds)
                        else:
                            cds_region.setdefault(j,[]).append(cds)
            cds = []
            for k in cds_region:
                start.append(k.split('-')[0])
                end.append(k.split('-')[1])
                cds.append(";".join(cds_region[k]))

            #unique_region = {}
            #for k in cds_region:
            #    if len(cds_region[k]) > 1:
            #        start.append(k.split('-')[0])
            #        end.append(k.split('-')[1])
            #        cds.append(cds_region[k][0])
            #    else:
            #        unique_region[k] = cds_region[k]
            #for k in unique_region:
            #    print k+"\t"+unique_region[k][0]
            #start_input = raw_input('please enter start site, seperated by ",":')
            #end_input = raw_input('please enter end site, seperated by ",":')
            #cds_input = raw_input('please enter cds ID, seperated by ",":')
            #start.extend(start_input.split(','))
            #end.extend(end_input.split(','))
            #cds.extend(cds_input.split(','))
                
        #hgmd = subprocess.Popen(['grep', gene, '/DATA/ana005/annovar/HGMD/HGMD_2018.spring/HGMD_Advanced_Micro_Lesions.csv'],stdout=subprocess.PIPE, stderr=subprocess.STDOUT).communicate()[0].split('\n')[0].split(',')[4].split(':')[0]
        #if not hgmd.startswith('NM'):
        #    hgmd = subprocess.Popen(['grep', gene, '/DATA/ana005/annovar/HGMD/HGMD_2018.spring/HGMD_Advanced_Micro_Lesions.csv'],stdout=subprocess.PIPE, stderr=subprocess.STDOUT).communicate()[0].split('\n')[1].split(',')[4].split(':')[0]
        hgmd = pd.read_csv('/DATA/sslyu/trio_BWA-GATK_v3.0/bed/gene/gene_hgvc.csv', header = None)
        hgmd = hgmd[hgmd[0] == gene]
        hgvs = hgmd[1].unique().tolist() # expand=True, return seperate columns
        print hgvs
        index = sorted(range(len(start)), key=lambda k: start[k])
        ofile = open('/DATA/sslyu/trio_BWA-GATK_v3.0/bed/gene/CDS/'+gene+'.bed','w')
        for i in index:
            ofile.write("chr"+chrom+"\t"+start[i]+"\t"+end[i]+"\t"+gene+"\t"+cds[i]+"\t"+";".join(hgvs)+"\n")
        ofile.close()

gene = sys.argv[1]
if gene.endswith('.bed'):
    gene = gene[:-4]
extract_cds_region(gene)
