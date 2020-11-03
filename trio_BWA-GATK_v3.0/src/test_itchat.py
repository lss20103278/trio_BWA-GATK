#!/usr/bin/env python
# -*- coding: utf-8 -*-
#########################################################################
# File Name: test_itchat.py
# Author: sslyu
# mail: 1573077420@qq.com
# Created Time:Fri 15 Feb 2019 10:03:59 AM CST
#########################################################################

import numpy as np
import pandas as pd

import itchat
itchat.auto_login()
itchat.send('Hello, filehelper', toUserName='filehelper')
