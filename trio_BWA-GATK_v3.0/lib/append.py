#!/usr/bin/env python
# -*- coding: utf-8 -*-
#########################################################################
# File Name: append.py
# Author: sslyu
# mail: 1573077420@qq.com
# Created Time:Fri 19 Apr 2019 03:43:50 PM CST
#########################################################################

import numpy as np
import pandas as pd
import os                                                                                                                       
import sys
reload(sys)
sys.setdefaultencoding("utf8")
import shutil

m_path = os.path.split(os.path.abspath(sys.argv[0]))[0]
print m_path

def append_sample_excel(df):
    ID = []
    if u'样本编号' in df.columns and u'原始样本ID' in df.columns:
        for i in range(len(df[u'原始样本ID'])):
            if df[u'原始样本ID'][i] == 'nan':
                ID.append(str(df.iloc[i][u'样本编号']))
            else:
                ID.append(str(df.iloc[i][u'原始样本ID']))
    df['sample'] = ID
    df['sample'] = df['sample'].apply(str)
    return df

def append_gender(df):
    df['gender'] = df['性别'].apply(lambda x:'1' if '男' in x else '2' if '女' in x else 'unknown')
    return df                                                                                                                   

def generate_pedigree(df):
    pedigree = {}
    df.index = df[u'姓名']
    for i in df[u'姓名']:
        pedigree[i] = []
    for k in pedigree:
        for i in df[u'姓名']:
            if k in i: # the name of the child should be totally included in other family member's names
                pedigree[k].append(i)
    list_keys = pedigree.keys()
    for k in list_keys:
        if len(pedigree[k]) == 1:
            pedigree.pop(k)
    df.index = df[u'姓名']
    return pedigree

def append_relation(df):
    df.index = df[u'姓名']                                                                                                       
    if 'relationship' not in df.columns:
        df['relationship'] = None
    pedigree_dict = generate_pedigree(df)
    trio = pedigree_dict.keys()
    for i in df.index:
        for j in trio:
            if j in i:
                if j == i:                                                                                                      
                    df.loc[i,'relationship'] = u'子'
                elif u'父' in i:
                    df.loc[i,'relationship'] = u'父'  #there is no '父' in i
                elif u'母' in i:
                    df.loc[i,'relationship'] = u'母'  #there is no '母' in i
                elif u'爷爷' in i:
                    df.loc[i,'relationship'] = u'爷爷'  #there is no '爷爷' in i
                elif u'奶奶' in i:
                    df.loc[i,'relationship'] = u'奶奶'  #there is no '奶奶' in i
                elif u'外公' in i:
                    df.loc[i,'relationship'] = u'外公'  #there is no '外公' in i
                elif u'外婆' in i:
                    df.loc[i,'relationship'] = u'外婆'  #there is no '外婆' in i
                elif u'姐姐' in i:
                    df.loc[i,'relationship'] = u'姐姐'  #there is no '姐姐' in i
                elif u'妹妹' in i:
                    df.loc[i,'relationship'] = u'妹妹'  #there is no '妹妹' in i
                elif u'哥哥' in i:
                    df.loc[i,'relationship'] = u'哥哥'  #there is no '哥哥' in i
                elif u'弟弟' in i:
                    df.loc[i,'relationship'] = u'弟弟'  #there is no '弟弟' in i
                else:
                    df.loc[i,'relationship'] = 'other'
    df[u'家系关系'] = None
    df[u'家系关系'] = df['relationship']
    return df

def append_pedigree(df):                                                                                                        
    pedigree_dict = generate_pedigree(df)
    df.index = df[u'姓名']
    df['familyname'] = df['sample']
    for i in df.index:
        for k in pedigree_dict:
            if k in i:  # the name of the child should be totally included in other family member's names
                df.loc[i,'familyname'] = df.loc[k]['sample']
    df[u'样本间关系'] = None
    for i in df.index:
        for k in pedigree_dict:
            if k in i:
                df.loc[i,u'样本间关系'] = 'fam'+df.loc[k]['sample']
    return df                    

def flatten_list(nested):
    if isinstance(nested, list):
        for sublist in nested:                                                                                                  
            for item in flatten_list(sublist):
                yield item
    else:
        yield nested

def append_phenotype(df):
    pedigree_dict = generate_pedigree(df)
    trio = pedigree_dict.keys()
    try:
        single = [i for i in df[u'姓名'] if i not in flatten_list(pedigree_dict.values())] # dict.values()
    except:
        pass                                                                                                                    
    df['phenotype1'] = '0'
    df['phenotype2'] = None
    if len(trio) > 0 and len(single) > 0:
        for i in df.index: 
            if i in trio:
                df.loc[i,'phenotype2'] = '2'
            elif i in single:
                df.loc[i,'phenotype2'] = '2'
            else:
                df.loc[i,'phenotype2'] = '1'
    elif len(trio) > 0:
        for i in df.index:
            if i in trio:
                df.loc[i,'phenotype2'] = '2'
            else:
                df.loc[i,'phenotype2'] = '1'
    else:
        for i in df.index:
            df.loc[i,'phenotype2'] = '2'                                                                                        
    return df

def generate_single(df):
    pedigree_dict = generate_pedigree(df)
    trio = pedigree_dict.keys()
    trio_m = []
    for i in df[u'姓名']:
        for k in pedigree_dict:
            trio_m.extend(pedigree_dict[k])
    single = [i for i in df[u'姓名'] if i not in trio_m]                                                                        
    return single

def append_father_mother(df):
    df['father'] = '0'
    df['mother'] = '0'
    pedigree_dict = generate_pedigree(df)
    trio = pedigree_dict.keys()
    try:
        single = [i for i in df[u'姓名'] if i not in pedigree_dict.values()[0]]
    except:                                                                                                                     
        single = df[u'姓名'].tolist()
    df.index = df[u'姓名']
    for i in trio: # i 患病孩子
        relation = ''
        for j in pedigree_dict[i]: 
            relation = relation+df.loc[j]['relationship']
        for j in df.index:
            if j in pedigree_dict[i]: 
                if u'父' in relation and u'母' in relation: #父母都有
                    if j == i: #j为患病孩子
                        for k in pedigree_dict[i]:
                            if u'父' in k:
                                father_name = k 
                            if u'母' in k:
                                mother_name = k
                        df.loc[j,'father'] = df.loc[father_name]['sample']
                        df.loc[j,'mother'] = df.loc[mother_name]['sample']
                    else: # j为患病孩子家系中其他人                                                                             
                        if len(pedigree_dict[i]) == 3: #家系结构为子父母，父母的父母为0
                            df.loc[j,'father'] = '0' # j为父或母
                        else: #家系结构为子父母和其他家人
                            if u'姐姐' in j or u'妹妹' in j or u'哥哥' in j or u'弟弟' in j: # 家系结构为 孩子，孩子的兄弟姐妹，父母，j为兄弟姐妹
                                for k in pedigree_dict[i]:
                                    if u'父' in k:
                                        father_name = k
                                    if u'母' in k:
                                        mother_name = k
                                df.loc[j,'father'] = df.loc[father_name]['sample']
                                df.loc[j,'mother'] = df.loc[mother_name]['sample']
                            elif u'爷爷' not in relation and u'奶奶' not in relation and u'外婆' not in relation and u'外公' not in relation: 
                                if u'父' in j:
                                    df.loc[j,'father'] = '0'
                                    df.loc[j,'mother'] = '0'
                                if u'母' in j:
                                    df.loc[j,'father'] = '0'
                                    df.loc[j,'mother'] = '0'
                            else:
                                father = raw_input('father ID of '+j+'(if not exist, please enter 0): ')
                                df.loc[j,'father'] = father
                                mother = raw_input('mother ID of '+j+'(if not exist, please enter 0):')
                                df.loc[j,'mother'] = mother
                else: # 父母一方没有
                    if j == i: # j为患病孩子
                        if u'父' in relation: # 只有父亲
                            for k in pedigree_dict[i]:                                                                          
                                if u'父' in k:
                                    father_name = k
                            df.loc[j,'father'] = df.loc[father_name]['sample']
                        if u'母' in relation: # 只有母亲
                            for k in pedigree_dict[i]:
                                if u'母' in k:
                                    mother_name = k
                            df.loc[j,'mother'] = df.loc[mother_name]['sample']
                    elif u'父' in j: # j为父亲
                        if u'爷爷' not in relation and u'奶奶' not in relation and u'外婆' not in relation and u'外公' not in relation: 
                            df.loc[j,'father'] = '0'
                            df.loc[j,'mother'] = '0'
                    elif u'母' in j: # j为母亲
                        if u'爷爷' not in relation and u'奶奶' not in relation and u'外婆' not in relation and u'外公' not in relation: 
                            df.loc[j,'father'] = '0'
                            df.loc[j,'mother'] = '0'
                    else: # j为除父母外的其他家人
                        father = raw_input('father ID of '+j+'(if not exist, please enter 0):')
                        df.loc[j,'father'] = father                                                                             
                        mother = raw_input('mother ID of '+j+'(if not exist, please enter 0):')
                        df.loc[j,'mother'] = mother
    return df                        

