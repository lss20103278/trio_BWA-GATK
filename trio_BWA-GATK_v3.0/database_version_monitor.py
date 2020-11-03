#!/usr/bin/env python
# -*- coding: utf-8 -*-
#########################################################################
# File Name: database_version_monitor.py
# Author: sslyu
# mail: 1573077420@qq.com
# Created Time:Wed 15 May 2019 03:15:19 PM CST
#########################################################################

import numpy as np
import pandas as pd

os.system('annotate_variation.pl --downdb avdblist --buildver hg19 --webfrom annovar .')
database = ["avsnp", "1000g", "gnomad", "exac03", "esp6500", "dbscsnv", "dbnsfp", "revel", "intervar", "clinvar", "ensembl", "hgmd"]

