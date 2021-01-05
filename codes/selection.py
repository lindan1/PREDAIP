# -*- coding: utf-8 -*-
"""
Created on Sat Dec 26 16:51:16 2020

@author: DIY
"""
import pandas as pd


def select(fused_feature):
    optimal_feature_index_Path = 'optimal.csv'
    df = pd.read_csv(optimal_feature_index_Path)
    pos_index = df.values
    pos_feature = []
    for i in range(len(pos_index)):
        pos_feature.append(pos_index[i][0])
    optimal_feature = fused_feature[:,pos_feature]
    
    return optimal_feature





#optimal_feature_index_Path = r'C:\Users\DIY\Desktop\PREDAIP\codes\optimal.csv'