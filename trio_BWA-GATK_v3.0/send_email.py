#!/usr/bin/env python
# -*- coding: utf-8 -*-
#########################################################################
# File Name: /DATA/sslyu/trio_BWA-GATK_v3.0/send_email.py
# Author: sslyu
# mail: 1573077420@qq.com
# Created Time:Mon 22 Apr 2019 08:00:01 PM CST
#########################################################################

import numpy as np
import pandas as pd
import sys
reload(sys)
sys.setdefaultencoding("utf8")
import os
import smtplib                                                                                                                  
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from email.mime.base import MIMEBase
from email import encoders

if not os.path.exists('config'):
    print """
    Examples: 
    prepare config
    format: (delimiter \\t)
    -sn excelfile 
    -panel   IDT-PANEL 
    -filtmode   hard 
    -rawdir    absolute_path_of_rawdata 
    -md5    absolute_path_of_rawdata_of_md5_file 
    -R1 suffix_of_raw_R1 
    -R2 suffix_of_raw_R2
    -outdir absolute_path_of_backup
    -analysis_data the dirname of the analysis data of the rawdata
    -redo   yes
    Note:
    make sure excel colnames contain 批次   收样日期    样本编号    原始样本ID    项目    性别    年龄    科室    浓度(ng/ul)    体积(ul)   总量(ug)    OD260/OD280 OD260/OD230 质检结果    上机日期    下机日期    Reads（M）   bases(Gb)  Q30
    the rawdata and md5 fileare in the same single dir 
    """
    sys.exit()                                                                                                                  

kwargs={'-sn':'', '-panel':'', '-filtmode':'', '-rawdir':'', '-md5':'', '-R1':'', '-R2':'', '-outdir':'', '-analysis_data':'', '-redo':''}

with open('config') as f:
    for l in f:
        l = l.strip('\n').split('\t')
        if len(l) != 0:
            for j in kwargs.keys():
                if l[0] == j:
                    kwargs.update({j:l[1]})
        
sn=kwargs.get('-sn')
panel=kwargs.get('-panel')
filtmode=kwargs.get('-filtmode')
rawdir = kwargs.get('-rawdir')
md5 = kwargs.get('-md5')
R1 = kwargs.get('-R1')
R2 = kwargs.get('-R2')
outdir = kwargs.get('-outdir')
analysis_data = kwargs.get('-analysis_data')
redo = kwargs.get('-redo')

def send_email(choose):
    dirname = '/'.join(os.getcwd().split('/')[-2:])
    for f in os.listdir('.'):
        if 'zip' in f:
            filename = f
    email_user='sslv@basepair.cn'
    email_send1='lss@sibs.ac.cn'
    email_send2=['knqin@basepair.cn', 'jjia@basepair.cn', 'gli@basepair.cn', 'yjsun@basepair.cn'] # in smtplib module, when multi recipients, need a list
    body = open('email.txt').readlines()
    subject = body[1][9:-1]
    print subject
    #subject = sn+'结果'
    msg = MIMEMultipart()
    msg['From'] = email_user
    if choose == 1:
        msg['To'] = ", ".join(email_send1) # in email module, need a string
    else:
        msg['To'] = ", ".join(email_send2) # in email module, need a string
    msg['Subject'] = subject
    body = ''.join(body)                                                                                                        
    msg.attach(MIMEText(body, 'plain', 'utf-8'))
    attachment = open(filename, 'rb')
    part = MIMEBase('application', 'octet-stream')
    part.set_payload((attachment).read())                                                                                       
    encoders.encode_base64(part)
    part.add_header('Content-Disposition', "attachment; filename="+filename)
    msg.attach(part)
    text = msg.as_string()
    mail = smtplib.SMTP('smtp.basepair.cn', 25)
    mail.ehlo()
    mail.starttls()
    mail.login(email_user, 'azby06078961!')
    if choose == 1:
        mail.sendmail(email_user, email_send1, text)
    else:
        mail.sendmail(email_user, email_send2, text)
    mail.quit()

send_email(2)               

