#!/usr/bin/env python
# -*- coding: utf-8 -*-
#########################################################################
# File Name: peddy.py
# Author: sslyu
# mail: 1573077420@qq.com
# Created Time:Tue 16 Apr 2019 07:00:32 PM CST
#########################################################################

import numpy as np
import pandas as pd
import sys
reload(sys)
sys.setdefaultencoding("utf8")
import os
import subprocess
import shutil
import time

def peddy_check():
    ofile = open('peddy_check.log', 'w')
    sex = pd.read_csv('peddy/pedtrio/results/peddy.sex_check.csv', dtype=str)
    for i in sex.index:
        if sex.loc[i]['error'] == 'True':
            ofile.write(sex.loc[i]['sample_id']+u'分析单上性别为'+sex.loc[i]['ped_sex']+u'，预测性别为'+sex.loc[i]['predicted_sex']+'\n')
    ped = pd.read_csv('peddy/pedtrio/results/peddy.ped_check.csv', dtype=str)
    for i in ped.index:
        if ped.loc[i]['parent_error'] == 'True':
            ofile.write(ped.loc[i]['sample_a']+' and '+ped.loc[i]['sample_b']+' are not child-parent\n')
    ofile.close()
    if os.path.getsize('peddy_check.log') == 0:
        os.remove('peddy_check.log')
    else:
        os.system('cat peddy_check.log')

def peddy():    
    if not os.path.exists('peddy'):
        os.mkdir('peddy')
    mendel = 0
    for i in os.listdir('ped'):
        if 'mendel' in i:
            mendel = mendel+1
    if mendel == 0:
        print 'no mendel.ped ...'
        sys.exit()
    if not os.path.exists('ped/peddy.ped'):
        for i in os.listdir('ped'):
            if 'mendel' in i:
                os.system('cat ped/'+i+' >> peddy_tmp')
        os.system('sort peddy_tmp |uniq |awk \'{OFS="\t"; print \'0\',$2,$3,$4,$5,$6}\' > ped/peddy.ped; rm peddy_tmp')
    shutil.copy('/DATA/sslyu/trio_BWA-GATK_v3.0/src/trio_cmd.sh', 'peddy/trio_cmd.sh')
    #shutil.copy('/DATA/sslyu/trio_BWA-GATK_v3.0/src/trio_cmd_gatk4.sh', 'peddy/trio_cmd.sh')
    os.chdir('peddy')
    p1 = subprocess.Popen(["sbatch", 'trio_cmd.sh'], stdout=subprocess.PIPE)
    jobID = p1.communicate()[0].strip('\n').split()[-1]
    os.chdir('..')
    while True:
        time.sleep(60)
        p1 = subprocess.Popen(["sacct", "-u", "sslyu"], stdout=subprocess.PIPE)
        p2 = subprocess.Popen(["grep",jobID], stdin=p1.stdout, stdout=subprocess.PIPE)
        p3 = subprocess.Popen(["grep","RUNNING"], stdin=p2.stdout, stdout=subprocess.PIPE)
        sacct = p3.communicate()[0]
        if sacct == '':
            os.system('grep -ci error peddy/slurm-'+jobID+'.out')
            break

if not os.path.exists('peddy/pedtrio/results/peddy.sex_check.csv') and not os.path.exists('peddy/pedtrio/results/peddy.ped_check.csv'):
    peddy()        
else:
    peddy_check()
