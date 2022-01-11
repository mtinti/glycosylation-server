# -*- coding: utf-8 -*-
"""
Created on Thu Dec 15 13:10:34 2016

@author: mtinti-x
"""
import pandas as pd
import re

def get_df(in_file ='', window = 0):
    df = pd.read_csv(in_file)
    df = df[df['window_hot_index']==window]
    return df

def filter_pep(in_df ='', extend = 10):
    p = re.compile('N.[ST]')
    temp = []
    for pep in in_df['window_seq']:
        #print pep
        pep = pep[extend:extend+3]
        if p.match(pep) != None:
            temp.append(1)
        else:
            temp.append(0)
    in_df['filter']=temp
    in_df = in_df[in_df['filter']==1]
    del in_df['filter']
    return in_df
        
    
    
def extend_df(in_df='', hot_index=10):    
    list_bonus = [(-5,0.35), (-4,0.45), (-3,0.35), (-2,0.35), (-1,0.7), (1,0.6), (3,0.25), (4,0.15), (5,0.2), (6,0.25), (7,0.1) ]
    for item in list_bonus:
        #print item
        temp = []
        for pep in in_df['window_seq']:
            if item[0]< hot_index:
                if pep[item[0]+hot_index] == 'D' or pep[item[0]+hot_index] == 'E':
                    temp.append(item[1])
                else:
                    temp.append(0.0)
            else:
                temp.append(0.0)
        in_df['bonus_'+str(item[0]+hot_index)]=temp
    
    in_df['bonus_all']=in_df[[n for n in in_df.columns if 'bonus_' in n]].sum(axis=1)
    in_df['bonus_max']=in_df[[n for n in in_df.columns if 'bonus_' in n][:-1]].max(axis=1)
    in_df['bonus_presence_E']=[1 if 'E' in pep[4:-4] else 0 for pep in in_df['window_seq']]
    in_df['bonus_presence_D']=[1 if 'D' in pep[4:-4] else 0 for pep in in_df['window_seq']]
    #in_df['bonus_presence_E']=[1 if 'E' in pep[0:4]+pep[-4:] else 0 for pep in in_df['window_seq']]
    #in_df['bonus_presence_D']=[1 if 'D' in pep[0:4]+pep[-4:] else 0 for pep in in_df['window_seq']]
    headers_left = ['bonus_'+str(n[0]+hot_index) for n in list_bonus[0:5]]
    headers_right = ['bonus_'+str(n[0]+hot_index) for n in list_bonus[5:]]                    
    in_df['bonus_left']=in_df[[n for n in in_df.columns if n in  headers_left ]].sum(axis=1)    
    in_df['bonus_right']=in_df[[n for n in in_df.columns if n in headers_right]].sum(axis=1) 
    return in_df