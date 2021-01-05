# -*- coding: utf-8 -*-
"""
Created on Sat Dec 26 17:27:42 2020

@author: DIY
"""
# example
#filePath = r'C:\Users\DIY\Desktop\example.txt'

def ExtractSeq(filePath): 
    fastas = open(filePath)
    lines = fastas.read().splitlines()
    sequence = []
    namelist = []
    for line in lines:
        if line[0] == '>':
            namelist.append(line)
        else:
            sequence.append(line)
        
    return sequence,namelist