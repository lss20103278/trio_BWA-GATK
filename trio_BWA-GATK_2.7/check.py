#!/usr/bin/env python
# -*- coding: utf-8 -*-
#########################################################################
# File Name: /DATA/sslyu/trio_BWA-GATK_2.7/check.py
# Author: sslyu
# mail: 1573077420@qq.com
# Created Time:Wed 08 Aug 2018 01:46:07 PM CST
#########################################################################

import os

ofile = open('check.sh', 'w')
ofile.write('for i in `find . -name \"slurm*\"`; do a=`grep -ci error $i`; if [ $a -gt 0 ]; then echo -e $i has errors; fi; done\n')
ofile.close()
os.system('bash check.sh')
os.system('rm check.sh')
os.system('mv trio/peddy/pedtrio/results/peddy.{sex_check.csv,ped_check.csv,html} .')
if os.path.exists('all.qcsum'):
	os.system('cat all.qcsum')
